/*****************************************************************************
 *
 * This MobilityDB code is provided under The PostgreSQL License.
 * Copyright (c) 2016-2025, Université libre de Bruxelles and MobilityDB
 * contributors
 *
 * MobilityDB includes portions of PostGIS version 3 source code released
 * under the GNU General Public License (GPLv2 or later).
 * Copyright (c) 2001-2025, PostGIS contributors
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose, without fee, and without a written
 * agreement is hereby granted, provided that the above copyright notice and
 * this paragraph and the following two paragraphs appear in all copies.
 *
 * IN NO EVENT SHALL UNIVERSITE LIBRE DE BRUXELLES BE LIABLE TO ANY PARTY FOR
 * DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
 * LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
 * EVEN IF UNIVERSITE LIBRE DE BRUXELLES HAS BEEN ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 * UNIVERSITE LIBRE DE BRUXELLES SPECIFICALLY DISCLAIMS ANY WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON
 * AN "AS IS" BASIS, AND UNIVERSITE LIBRE DE BRUXELLES HAS NO OBLIGATIONS TO
 * PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 *
 *****************************************************************************/

/**
 * @file
 * @brief Example: read AIS-like CSV, build TSEQUENCE per MMSI, clean with EKF.
 *
 * Input CSV columns (comma separated), first row is a header:
 * Timestamp,Type of mobile,MMSI,Latitude,Longitude, ... (other columns ignored)
 * Example timestamp format: "01/03/2024 00:00:00" (DD/MM/YYYY HH:MM:SS)
 *
 * For each MMSI we aggregate geographic positions (Longitude, Latitude)
 * into a TSEQUENCE of type T_TDOUBLE2 (values are [lon, lat] in degrees),
 * and apply a Constant-Velocity EKF to clean the trajectory.
 *
 * Build (requires MEOS built with EKF enabled):
 *  gcc -Wall -O2 -DMEOS_EXPERIMENTAL_ANALYTICS \
 *    -I/usr/local/include -o ais_ekf_clean meos/examples/ais_ekf_clean.c \
 *    -L/usr/local/lib -lmeos
 *
 * Run:
 *  ./ais_ekf_clean input.csv
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#include <meos.h>
#include <meos_internal.h>
/* Bring only the minimal temporal definitions to avoid geo deps */
#include <temporal/doublen.h>
#include <temporal/meos_catalog.h>
/* Provide light macros usually supplied by temporal/temporal.h */
#ifndef DatumGetDouble2P
#define DatumGetDouble2P(X)       ((double2 *) DatumGetPointer(X))
#endif
#ifndef Double2PGetDatum
#define Double2PGetDatum(X)       PointerGetDatum(X)
#endif
/* no need for internal geo helpers in this example */


/* EKF header is experimental; include only if available */
#if defined(MEOS_EXPERIMENTAL_ANALYTICS) && (__has_include(<temporal/tekf.h>))
#include <temporal/tekf.h>
#define HAVE_TEKF 1
#else
#undef HAVE_TEKF
#endif

/* Limits for this small example (adjust if needed) */
#define MAX_SHIPS         1024
#define INITIAL_INSTANTS  64
#define MAX_LINE_LEN      2048

typedef struct
{
  TimestampTz t;
  long long mmsi;
  double lon;
  double lat;
} Rec;

typedef struct
{
  long long mmsi;
  int n_inst;
  int cap_inst;
  TInstant **inst; /* T_TDOUBLE2 instants [lon,lat] */
} Track;

static bool parse_timestamp_eu(const char *s, TimestampTz *out)
{
  /* Expect DD/MM/YYYY HH:MM:SS (no timezone). Interpret as UTC. */
  int d=0,m=0,y=0,hh=0,mm=0,ss=0;
  if (!s || strlen(s) < 10)
    return false;
  /* Allow both ' ' and 'T' as separator */
  int matched = sscanf(s, "%d/%d/%d %d:%d:%d", &d, &m, &y, &hh, &mm, &ss);
  if (matched < 3)
    return false;
  if (matched < 6)
  {
    hh = mm = ss = 0;
  }
  char iso[32];
  snprintf(iso, sizeof(iso), "%04d-%02d-%02d %02d:%02d:%02d", y, m, d, hh, mm, ss);
  TimestampTz ts = pg_timestamptz_in(iso, -1);
  if (ts == DT_NOEND)
    return false;
  *out = ts;
  return true;
}

static bool is_nonempty_token(const char *tok)
{
  return tok && *tok && strcmp(tok, "Unknown") != 0;
}

int main(int argc, char **argv)
{
  if (argc < 2)
  {
    fprintf(stderr, "Usage: %s input.csv [output.csv]\n", argv[0]);
    return EXIT_FAILURE;
  }

  /* Initialize MEOS */
  meos_initialize();
  meos_initialize_timezone("UTC");

  FILE *fp = fopen(argv[1], "r");
  if (!fp)
  {
    fprintf(stderr, "Could not open %s\n", argv[1]);
    meos_finalize();
    return EXIT_FAILURE;
  }

  /* Skip header */
  char line[MAX_LINE_LEN];
  if (!fgets(line, sizeof(line), fp))
  {
    fprintf(stderr, "Empty file\n");
    fclose(fp);
    meos_finalize();
    return EXIT_FAILURE;
  }

  Track tracks[MAX_SHIPS] = {0};
  int n_tracks = 0;
  int n_rows = 0, n_rows_ok = 0;

  while (fgets(line, sizeof(line), fp))
  {
    n_rows++;
    /* Tokenize the minimal needed columns */
    /* Columns:
     * 0: Timestamp
     * 1: Type of mobile
     * 2: MMSI
     * 3: Latitude
     * 4: Longitude
     */
    char *p = line;
    char *tokens[5] = {0};
    int idx = 0;
    for (char *tok = strtok(p, ","); tok && idx < 5; tok = strtok(NULL, ","))
      tokens[idx++] = tok;

    if (idx < 5)
      continue; /* incomplete */

    /* Filter: require MMSI, lat, lon, and a parseable timestamp */
    if (!is_nonempty_token(tokens[0]) || !is_nonempty_token(tokens[2]) ||
        !is_nonempty_token(tokens[3]) || !is_nonempty_token(tokens[4]))
      continue;

    /* Optional: ignore Base Station lines */
    if (tokens[1] && strncmp(tokens[1], "Base Station", 12) == 0)
      continue;

    Rec r = {0};
    if (!parse_timestamp_eu(tokens[0], &r.t))
      continue;
    r.mmsi = strtoll(tokens[2], NULL, 10);
    r.lat = strtod(tokens[3], NULL);
    r.lon = strtod(tokens[4], NULL);

    /* Find or add track for MMSI */
    int j = -1;
    for (int i = 0; i < n_tracks; i++)
    {
      if (tracks[i].mmsi == r.mmsi) { j = i; break; }
    }
    if (j == -1)
    {
      if (n_tracks == MAX_SHIPS)
        continue; /* drop if capacity exceeded */
      j = n_tracks++;
      tracks[j].mmsi = r.mmsi;
      tracks[j].n_inst = 0;
      tracks[j].cap_inst = INITIAL_INSTANTS;
      tracks[j].inst = (TInstant **) calloc((size_t)INITIAL_INSTANTS, sizeof(TInstant *));
      if (!tracks[j].inst)
      {
        fprintf(stderr, "Out of memory\n");
        fclose(fp);
        meos_finalize();
        return EXIT_FAILURE;
      }
    }

    /* Create T_TDOUBLE2 instant [lon, lat] at timestamp */
    double2 d2;
    double2_set(r.lon, r.lat, &d2);
    TInstant *ti = tinstant_make(Double2PGetDatum(&d2), T_TDOUBLE2, r.t);

    /* Append, expanding if needed; drop exact-duplicate timestamp */
    int n = tracks[j].n_inst;
    if (n > 0 && tracks[j].inst[n-1]->t == ti->t)
    {
      pfree(ti);
      continue;
    }
    if (n == tracks[j].cap_inst)
    {
      int newcap = tracks[j].cap_inst * 2;
      TInstant **tmp = (TInstant **) realloc(tracks[j].inst, (size_t)newcap * sizeof(TInstant *));
      if (!tmp)
      {
        fprintf(stderr, "Out of memory (expand)\n");
        pfree(ti);
        break;
      }
      tracks[j].inst = tmp;
      tracks[j].cap_inst = newcap;
    }
    tracks[j].inst[tracks[j].n_inst++] = ti;
    n_rows_ok++;
  }
  fclose(fp);

  printf("Read %d rows, accepted %d points, built %d tracks.\n",
         n_rows, n_rows_ok, n_tracks);

  /* Open output CSV */
  const char *outpath = (argc >= 3 ? argv[2] : "ais_ekf_clean_out.csv");
  FILE *fout = fopen(outpath, "w");
  if (!fout)
  {
    fprintf(stderr, "Could not open output file %s\n", outpath);
    /* continue but only print to stdout */
  }
  else
  {
    fprintf(fout, "Timestamp,MMSI,Longitude_Raw,Latitude_Raw,Longitude_Clean,Latitude_Clean\n");
  }

  /* Build sequences and clean with EKF (CV over lon/lat) */
  for (int i = 0; i < n_tracks; i++)
  {
    if (tracks[i].n_inst == 0)
      continue;

    TSequence *seq = tsequence_make((const TInstant **) tracks[i].inst,
                                    tracks[i].n_inst,
                                    true, true, LINEAR, true);

#ifdef HAVE_TEKF
    TEkfModel model;
    if (!tekf_make_cv_model(2, &model))
    {
      printf("MMSI %lld: EKF model creation failed, skipping.\n", tracks[i].mmsi);
      free(seq);
      continue;
    }
    /* CV context and parameters */
    TEkfCvCtx cvctx = { .D = 2, .q_accel_var = 1e-8, .r_meas_var = 1e-6 };
    double P0_diag[4] = { 1e-4, 1e-2, 1e-4, 1e-2 }; /* pos/vel per axis */
    TEkfParams params = {
      .default_dt = 1.0,
      .gate_sigma = 3.5,
      .fill_estimates = true,
      .P0_diag = P0_diag,
      .x0 = NULL,
      .Q_diag = NULL,
      .R_diag = NULL
    };
    int removed = 0;
    Temporal *clean = temporal_ekf_clean((Temporal *) seq, &model, &params, &cvctx, &removed);

    /* Dump CSV rows for this track */
    const TSequence *cseq = (const TSequence *) clean;
    int nraw = seq->count;
    int ncln = cseq ? cseq->count : 0;
    int nmin = (cseq ? (nraw < ncln ? nraw : ncln) : nraw);
    for (int k = 0; k < nmin; k++)
    {
      const TInstant *ir = TSEQUENCE_INST_N(seq, k);
      const TInstant *ic = cseq ? TSEQUENCE_INST_N(cseq, k) : NULL;
      char *ts = pg_timestamptz_out(ir->t);
      const double2 *vr = DatumGetDouble2P(tinstant_value_p(ir));
      double lon_r = vr->a, lat_r = vr->b;
      double lon_c = lon_r, lat_c = lat_r;
      if (ic)
      {
        const double2 *vc = DatumGetDouble2P(tinstant_value_p(ic));
        lon_c = vc->a; lat_c = vc->b;
      }
      if (fout)
        fprintf(fout, "%s,%lld,%.10f,%.10f,%.10f,%.10f\n", ts, tracks[i].mmsi, lon_r, lat_r, lon_c, lat_c);
      pfree(ts);
    }

    /* Optional: also print a brief summary to stdout */
    printf("MMSI %lld cleaned (%d removed).\n", tracks[i].mmsi, removed);

    pfree(seq);
    if (clean) pfree(clean);
#else
    (void) seq;
    /* Without EKF, still dump raw as both raw and clean */
    if (fout)
    {
      for (int k = 0; k < seq->count; k++)
      {
        const TInstant *ir = TSEQUENCE_INST_N(seq, k);
        const double2 *vr = DatumGetDouble2P(tinstant_value_p(ir));
        char *ts = pg_timestamptz_out(ir->t);
        fprintf(fout, "%s,%lld,%.10f,%.10f,%.10f,%.10f\n", ts, tracks[i].mmsi, vr->a, vr->b, vr->a, vr->b);
        pfree(ts);
      }
    }
    printf("MMSI %lld: built sequence with %d instants. Rebuild MEOS with -DMEOS_EXPERIMENTAL_ANALYTICS to enable EKF.\n",
           tracks[i].mmsi, tracks[i].n_inst);
    pfree(seq);
#endif
  }

  if (fout) fclose(fout);

  /* Free allocated instants */
  for (int i = 0; i < n_tracks; i++)
  {
    for (int j = 0; j < tracks[i].n_inst; j++)
      pfree(tracks[i].inst[j]);
    free(tracks[i].inst);
  }

  meos_finalize();
  return EXIT_SUCCESS;
}
