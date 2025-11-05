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
 * @brief Example: read AIS CSV, build TSEQUENCE per MMSI, clean with EKF.
 *
 * Input CSV columns, first row is a header:
 * Timestamp,Type of mobile,MMSI,Latitude,Longitude, ... (other columns ignored)
 * Example timestamp format: "01/03/2024 00:00:00" (DD/MM/YYYY HH:MM:SS)
 *
 * For each MMSI we aggregate geographic positions (Longitude, Latitude)
 * into a TSEQUENCE of (values are [lon, lat] in degrees),
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
#include <math.h>

#include <meos.h>
/* Bring only the minimal temporal definitions to avoid geo deps */
#include <temporal/doublen.h>
#include <temporal/meos_catalog.h>
/* Provide light macros and alloc fallbacks when not building with MEOS */
#ifndef DatumGetDouble2P
#define DatumGetDouble2P(X)       ((double2 *) DatumGetPointer(X))
#endif
#ifndef Double2PGetDatum
#define Double2PGetDatum(X)       PointerGetDatum(X)
#endif
#ifndef palloc
#include <stdlib.h>
#define palloc(sz)                malloc(sz)
#define palloc0(sz)               calloc(1, (sz))
#define pfree(ptr)                free(ptr)
#endif

/* Minimal internal forward decls we rely on (avoid meos_internal.h deps) */
extern TInstant *tinstant_make(Datum value, meosType temptype, TimestampTz t);
extern Datum tinstant_value_p(const TInstant *inst);
extern TInstant *temporal_instant_n(const Temporal *temp, int n); /* returns a copy */

/* Simple local ENU (East, North) approximation around an origin in meters */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
static inline double deg2rad(double d){ return d * (M_PI / 180.0); }
static inline double rad2deg(double r){ return r * (180.0 / M_PI); }

static void enu_from_lonlat(double lat0_deg, double lon0_deg,
                            double lat_deg, double lon_deg,
                            double *E, double *N)
{
  const double R = 6378137.0; /* WGS84 equatorial radius (m) */
  double lat0 = deg2rad(lat0_deg), lon0 = deg2rad(lon0_deg);
  double lat  = deg2rad(lat_deg),  lon  = deg2rad(lon_deg);
  double dlat = lat - lat0;
  double dlon = lon - lon0;
  double clat = cos(lat0);
  *E = R * clat * dlon;
  *N = R * dlat;
}

static void lonlat_from_enu(double lat0_deg, double lon0_deg,
                            double E, double N,
                            double *lat_deg, double *lon_deg)
{
  const double R = 6378137.0;
  double lat0 = deg2rad(lat0_deg), lon0 = deg2rad(lon0_deg);
  double clat = cos(lat0);
  double dlat = N / R;
  double dlon = (clat > 1e-12 ? E / (R * clat) : 0.0);
  double lat = lat0 + dlat;
  double lon = lon0 + dlon;
  *lat_deg = rad2deg(lat);
  *lon_deg = rad2deg(lon);
}
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

/* ENU adapters for EKF model */
typedef struct { int D; double q_accel_var; double r_meas_var; double lat0, lon0; } CvEnuCtx;

static bool enu_z_from_value(Datum v, meosType tt, double *z, void *ctx)
{
  CvEnuCtx *c = (CvEnuCtx *) ctx;
  meosType bt = temptype_basetype(tt);
  if (bt != T_DOUBLE2) return false;
  const double2 *d = DatumGetDouble2P(v);
  double E, N;
  /* d->a = lon, d->b = lat */
  enu_from_lonlat(c->lat0, c->lon0, d->b, d->a, &E, &N);
  z[0] = E; z[1] = N;
  return true;
}

static bool enu_value_from_state(const double *x, meosType tt, Datum *out, void *ctx)
{
  CvEnuCtx *c = (CvEnuCtx *) ctx;
  meosType bt = temptype_basetype(tt);
  if (bt != T_DOUBLE2) return false;
  /* CV state: [E, vE, N, vN] */
  double E = x[0], N = x[2];
  double lat, lon;
  lonlat_from_enu(c->lat0, c->lon0, E, N, &lat, &lon);
  double2 *d2p = palloc(sizeof(double2));
  double2_set(lon, lat, d2p);
  *out = Double2PGetDatum(d2p);
  return true;
}

static bool parse_timestamp_eu(const char *s, TimestampTz *out)
{
  /* Expect DD/MM/YYYY HH:MM:SS (no timezone). Interpret as UTC. */
  int d=0,m=0,y=0,hh=0,mm=0,ss=0;
  if (!s || strlen(s) < 10)
    return false;
  /* Note: this parser expects a space as separator (" ") */
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
#ifdef DT_NOEND
  if (ts == DT_NOEND)
    return false;
#endif
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
  /* Optional third/extra args:
     - 'fill' (default) or 'drop'
     - numeric gate_sigma
     - q=... (process accel var)
     - r=... (measurement var)
     - 'enu' (filter in meters, default) or 'deg' (filter in degrees)
  */
  bool cli_drop = false;
  double cli_gate_sigma = -1.0;
  double cli_q = -1.0;      /* q_accel_var (units depend: m^2/s^4 if ENU, deg^2/s^4 if degrees) */
  double cli_r = -1.0;      /* r_meas_var (units depend: m^2 if ENU, deg^2 if degrees) */
  bool use_enu = true;
  for (int ai = 3; ai < argc; ai++)
  {
    if (strcmp(argv[ai], "drop") == 0) cli_drop = true;
    else if (strcmp(argv[ai], "fill") == 0) cli_drop = false;
    else if (strncmp(argv[ai], "q=", 2) == 0) { cli_q = strtod(argv[ai]+2, NULL); }
    else if (strncmp(argv[ai], "r=", 2) == 0) { cli_r = strtod(argv[ai]+2, NULL); }
    else if (strcmp(argv[ai], "enu") == 0) { use_enu = true; }
    else if (strcmp(argv[ai], "deg") == 0) { use_enu = false; }
    else
    {
      char *endp = NULL; double v = strtod(argv[ai], &endp);
      if (endp && endp != argv[ai]) cli_gate_sigma = v;
    }
  }
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

  /* Build sequences and clean with EKF (CV over lon/lat degrees, or ENU meters if enabled) */
  for (int i = 0; i < n_tracks; i++)
  {
    if (tracks[i].n_inst == 0)
      continue;

    /* Build raw sequence without normalization to preserve all instants */
    TSequence *seq = tsequence_make((const TInstant **) tracks[i].inst,
                                    tracks[i].n_inst,
                                    true, true, LINEAR, false);

#ifdef HAVE_TEKF
    TEkfModel model;
    if (!tekf_make_cv_model(2, &model))
    {
      printf("MMSI %lld: EKF model creation failed, skipping.\n", tracks[i].mmsi);
      pfree(seq);
      continue;
    }
    /* Option A (deg): use input lon/lat directly.
       Option B (enu, default): run the EKF in ENU meters via model adapters.
       We keep the sequence as-is and adapt measurements/state. */
    double lat0 = 0.0, lon0 = 0.0;
    if (use_enu)
    {
      TInstant *i0 = temporal_instant_n((const Temporal *) seq, 1); /* 1-based */
      const double2 *v0 = DatumGetDouble2P(tinstant_value_p(i0));
      lon0 = v0->a; lat0 = v0->b;
      pfree(i0);
    }
    /* CV context and parameters */
    /* Defaults: if ENU (meters): q in m^2/s^4, r in m^2; if degrees: q, r are in degrees units */
    double def_q = use_enu ? 0.04 /* (0.2 m/s^2)^2 */ : 5e-10;
    double def_r = use_enu ? 400.0 /* (20 m)^2 */ : 4e-6;
    CvEnuCtx cvctx = { .D = 2,
                       .q_accel_var = (cli_q > 0 ? cli_q : def_q),
                       .r_meas_var  = (cli_r > 0 ? cli_r : def_r),
                       .lat0 = lat0, .lon0 = lon0 };
    double P0_diag[4] = { 1e-4, 1e-2, 1e-4, 1e-2 }; /* pos/vel per axis */
    TEkfParams params = {
      .default_dt = 1.0, /* used when two observations have non-positive dt */
      .gate_sigma = (cli_gate_sigma > 0 ? cli_gate_sigma : 8.0), /* Mahalanobis threshold in sigmas (0 disables) */
      .fill_estimates = !cli_drop, /* when gated out, emit estimate instead of dropping */
      .P0_diag = P0_diag,
      .x0 = NULL,
      .Q_diag = NULL,
      .R_diag = NULL
    };
    int removed = 0;
    /* When using ENU, adapt model to convert measurements/state */
    if (use_enu)
    {
      /* Measurement: z = [E,N] from [lon,lat] */
      model.z_from_value = enu_z_from_value;
      /* State -> output value: [E,N] -> [lon,lat] */
      model.value_from_state = enu_value_from_state;
    }
    Temporal *clean = temporal_ekf_clean((Temporal *) seq, &model, &params,
                                         (void*)&cvctx,
                                         &removed);

    /* Dump CSV rows for this track */
    const TSequence *cseq = (const TSequence *) clean;
    int nraw = seq->count;
    int ncln = cseq ? cseq->count : 0;
    if (cseq && !params.fill_estimates)
    {
      /* Two-pointer join on timestamps: cleaned times are subset of raw */
      int iraw = 0, icln = 0;
      while (iraw < nraw && icln < ncln)
      {
        TInstant *ir = temporal_instant_n((const Temporal *) seq, iraw + 1);
        TInstant *ic = temporal_instant_n((const Temporal *) cseq, icln + 1);
        if (ir->t == ic->t)
        {
          char *ts = pg_timestamptz_out(ir->t);
          const double2 *vr = DatumGetDouble2P(tinstant_value_p(ir));
          double lon_r = vr->a, lat_r = vr->b;
          const double2 *vc = DatumGetDouble2P(tinstant_value_p(ic));
          double lon_c = vc->a, lat_c = vc->b;
          if (fout)
            fprintf(fout, "%s,%lld,%.10f,%.10f,%.10f,%.10f\n", ts, tracks[i].mmsi,
                    lon_r, lat_r, lon_c, lat_c);
          pfree(ts);
          pfree(ir); pfree(ic);
          iraw++; icln++;
        }
        else if (ir->t < ic->t)
        { pfree(ic); pfree(ir); iraw++; }
        else
        { pfree(ic); pfree(ir); icln++; }
      }
    }
    else
    {
      int nmin = (cseq ? (nraw < ncln ? nraw : ncln) : nraw);
      for (int k = 0; k < nmin; k++)
      {
        TInstant *ir = temporal_instant_n((const Temporal *) seq, k + 1);
        TInstant *ic = cseq ? temporal_instant_n((const Temporal *) cseq, k + 1) : NULL;
        char *ts = pg_timestamptz_out(ir->t);
        const double2 *vr = DatumGetDouble2P(tinstant_value_p(ir));
        double lon_r = vr->a, lat_r = vr->b;
        double lon_c = lon_r, lat_c = lat_r;
        if (ic) { const double2 *vc = DatumGetDouble2P(tinstant_value_p(ic)); lon_c = vc->a; lat_c = vc->b; }
        if (fout)
          fprintf(fout, "%s,%lld,%.10f,%.10f,%.10f,%.10f\n", ts, tracks[i].mmsi, lon_r, lat_r, lon_c, lat_c);
        pfree(ts);
        if (ic) pfree(ic);
        pfree(ir);
      }
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
        TInstant *ir = temporal_instant_n((const Temporal *) seq, k + 1);
        const double2 *vr = DatumGetDouble2P(tinstant_value_p(ir));
        char *ts = pg_timestamptz_out(ir->t);
        fprintf(fout, "%s,%lld,%.10f,%.10f,%.10f,%.10f\n", ts, tracks[i].mmsi, vr->a, vr->b, vr->a, vr->b);
        pfree(ts);
        pfree(ir);
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
