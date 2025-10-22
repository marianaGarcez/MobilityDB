/*
 * AIS EKF demo: read AIS CSV, derive speed/steering, run tgeompoint_ekf.
 *
 * CSV columns (minimum used):
 *  0 Timestamp (dd/mm/yyyy HH:MM:SS)
 *  1 Type of mobile (use "Class A")
 *  2 MMSI (filter by argument)
 *  3 Latitude
 *  4 Longitude
 *  7 SOG (knots)
 *  8 COG (degrees)
 *  6 ROT (deg/min) [optional]
 *
 * Build (in-tree):
 *  gcc -O2 -Imeos/include -Lbuild -lmeos meos/examples/ais_ekf.c -o ais_ekf
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <time.h>

/* Prefer generated export headers when building in-tree */
#if __has_include("meos_export.h")
#  include "meos_export.h"
#else
#  include <meos.h>
#endif
#if __has_include("meos_geo_export.h")
#  include "meos_geo_export.h"
#else
#  include <meos_geo.h>
#endif

#ifndef pfree
#define pfree free
#endif
#include <meos_internal_geo.h>
#include <geo/tgeo_ekf.h>

#define KNOTS_TO_MS 0.514444

typedef struct {
  TimestampTz t;
  double lat, lon;    /* degrees */
  double sog_kn;      /* knots */
  double cog_deg;     /* degrees [0,360) */
  double rot_deg_min; /* degrees per minute; NAN if missing */
  int has_rot;
} AisRow;

/* ALL-mode row with MMSI field */
typedef struct {
  char mmsi[32];
  TimestampTz t;
  double lat, lon;    /* degrees */
  double sog_kn;      /* knots */
  double cog_deg;     /* degrees [0,360) */
  double rot_deg_min; /* degrees per minute; NAN if missing */
  int has_rot;
} AisRowAll;

static int parse_ddmmyyyy(char *dst, size_t dst_sz, const char *src)
{
  /* src: dd/mm/yyyy HH:MM:SS -> dst: yyyy-mm-dd HH:MM:SS */
  int d,m,y,hh,mm,ss;
  if (sscanf(src, "%d/%d/%d %d:%d:%d", &d,&m,&y,&hh,&mm,&ss) != 6)
    return 0;
  snprintf(dst, dst_sz, "%04d-%02d-%02d %02d:%02d:%02d", y,m,d,hh,mm,ss);
  return 1;
}

static double unwrap_angle(double prev, double cur)
{
  /* both in radians; unwrap cur near prev */
  double x = cur;
  while (x - prev > M_PI)  x -= 2.0*M_PI;
  while (x - prev < -M_PI) x += 2.0*M_PI;
  return x;
}

static void ll_to_xy(double lat0_deg, double lon0_deg, double lat_deg, double lon_deg,
                     double *x, double *y)
{
  /* Equirectangular projection around origin (m) */
  const double R = 6378137.0;
  double lat0 = lat0_deg * M_PI/180.0;
  double dlat = (lat_deg - lat0_deg) * M_PI/180.0;
  double dlon = (lon_deg - lon0_deg) * M_PI/180.0;
  *x = R * dlon * cos(lat0);
  *y = R * dlat;
}

static void xy_to_ll(double lat0_deg, double lon0_deg, double x, double y,
                     double *lat_deg, double *lon_deg)
{
  const double R = 6378137.0;
  double lat0 = lat0_deg * M_PI/180.0;
  double dlat = y / R;
  double dlon = (fabs(cos(lat0)) > 1e-12) ? (x / (R * cos(lat0))) : 0.0;
  *lat_deg = lat0_deg + (dlat * 180.0/M_PI);
  *lon_deg = lon0_deg + (dlon * 180.0/M_PI);
}

static int split_csv(char *line, char **fields, int max_fields)
{
  /* Very simple CSV split by comma; empty fields allowed */
  int n = 0;
  char *p = line;
  while (*p && n < max_fields) {
    fields[n++] = p;
    char *c = strchr(p, ',');
    if (!c) break;
    *c = '\0';
    p = c + 1;
  }
  return n;
}

/* Parse either dd/mm/yyyy HH:MM:SS or ISO yyyy-mm-dd HH:MM:SS into timestamptz */
static int parse_ts_any(TimestampTz *out, const char *src)
{
  if (strchr(src, '/')) {
    char iso[64];
    if (!parse_ddmmyyyy(iso, sizeof(iso), src)) return 0;
    *out = timestamptz_in(iso, -1);
    return 1;
  } else {
    *out = timestamptz_in(src, -1);
    return 1;
  }
}

/* Compare AisRow by timestamp */
static int cmp_rows(const void *a, const void *b)
{
  const AisRow *ra = (const AisRow *)a;
  const AisRow *rb = (const AisRow *)b;
  if (ra->t < rb->t) return -1;
  if (ra->t > rb->t) return 1;
  return 0;
}

/* Compare by MMSI then timestamp */
static int cmp_mmsi_t(const void *a, const void *b)
{
  const AisRowAll *ra = (const AisRowAll *)a;
  const AisRowAll *rb = (const AisRowAll *)b;
  int c = strcmp(ra->mmsi, rb->mmsi);
  if (c) return c;
  if (ra->t < rb->t) return -1;
  if (ra->t > rb->t) return 1;
  return 0;
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <MMSI> [csv_file]\n", argv[0]);
    return EXIT_FAILURE;
  }
  const char *mmsi_filter = argv[1];
  FILE *fin = (argc >= 3) ? fopen(argv[2], "r") : stdin;
  int progress = getenv("AIS_EKF_PROGRESS") != NULL;
  long prog_every_rows = 200000;
  const char *s;
  if ((s = getenv("AIS_EKF_PROGRESS_EVERY"))) { long v = atol(s); if (v > 0) prog_every_rows = v; }
  if (!fin) { perror("fopen"); return EXIT_FAILURE; }

  meos_initialize();
  /* Ensure GSL never aborts on numeric errors (Cholesky fallback works) */
  gsl_set_error_handler_off();

  /* ALL MMSIs mode: ./ais_ekf ALL file.csv */
  if (strcmp(mmsi_filter, "ALL") == 0) {
    char buf[4096];
    /* Skip header if present */
    if (fgets(buf, sizeof(buf), fin) == NULL) { fprintf(stderr, "Empty input\n"); return EXIT_FAILURE; }
    if (strstr(buf, "Timestamp") == NULL) fseek(fin, 0, SEEK_SET);

    int progress = getenv("AIS_EKF_PROGRESS") != NULL;
    long prog_every_rows = 200000; const char *s;
    if ((s = getenv("AIS_EKF_PROGRESS_EVERY"))) { long v = atol(s); if (v > 0) prog_every_rows = v; }

    AisRowAll *rows = NULL; size_t nrows = 0, cap = 0; size_t lines_read = 0;
    while (fgets(buf, sizeof(buf), fin)) {
      lines_read++;
      if (progress && (lines_read % (size_t)prog_every_rows == 0)) {
        fprintf(stderr, "[ais_ekf] read %zu lines...\n", lines_read);
        fflush(stderr);
      }
      char line[4096]; strncpy(line, buf, sizeof(line)); line[sizeof(line)-1] = '\0';
      char *f[32]; int nf = split_csv(line, f, 32);
      if (nf < 10) continue;
      if (strcmp(f[1], "Class A") != 0) continue;
      if (!*f[3] || !*f[4]) continue;
      AisRowAll r = {0};
      strncpy(r.mmsi, f[2], sizeof(r.mmsi)-1);
      if (!parse_ts_any(&r.t, f[0])) continue;
      r.lat = atof(f[3]); r.lon = atof(f[4]);
      r.sog_kn = (*f[7]) ? atof(f[7]) : 0.0;
      r.cog_deg = (*f[8]) ? atof(f[8]) : 0.0;
      r.rot_deg_min = (*f[6]) ? atof(f[6]) : NAN;
      r.has_rot = (*f[6] && *f[6] != '\0');
      if (nrows == cap) { cap = cap ? cap*2 : 8192; rows = realloc(rows, cap * sizeof(*rows)); }
      rows[nrows++] = r;
    }
    if (fin != stdin) fclose(fin);
    if (nrows == 0) { fprintf(stderr, "No AIS rows\n"); meos_finalize(); return EXIT_FAILURE; }

    qsort(rows, nrows, sizeof(AisRowAll), cmp_mmsi_t);
    printf("t,lon_raw,lat_raw,lon_filt,lat_filt,sog_kn,cog_deg,rot_deg_min\n");

    size_t start = 0;
    while (start < nrows) {
      size_t end = start + 1; while (end < nrows && strcmp(rows[end].mmsi, rows[start].mmsi) == 0) end++;
      size_t gn = end - start; if (gn < 2) { start = end; continue; }

      TInstant **speed_inst = malloc(gn * sizeof(*speed_inst));
      TInstant **steer_inst = malloc(gn * sizeof(*steer_inst));
      TInstant **gps_inst   = malloc(gn * sizeof(*gps_inst));

      double lat0 = rows[start].lat, lon0 = rows[start].lon;
      GSERIALIZED *origin = geompoint_make2d(0, 0.0, 0.0);
      double phi_prev = rows[start].cog_deg * M_PI/180.0;
      TGeoEkfParams p = { .L=2.83,.H=0.76,.A=3.78,.B=0.5,.q_x=0.5,.q_y=0.5,.q_phi=0.5,.r_x=0.5,.r_y=5.0,.p0_x=10.0,.p0_y=10.0,.p0_phi=0.5,.epsS=0.0,.hold_last_controls=true };

      size_t w = 0; TimestampTz prev_t = 0; double prev_cog = 0.0; int first = 1;
      for (size_t j = start; j < end; j++) {
        TimestampTz t = rows[j].t; if (!first && t <= prev_t) continue;
        double x, y; ll_to_xy(lat0, lon0, rows[j].lat, rows[j].lon, &x, &y);
        GSERIALIZED *pt = geompoint_make2d(0, x, y); gps_inst[w] = tpointinst_make(pt, t); pfree(pt);
        double v_ms = rows[j].sog_kn * KNOTS_TO_MS; speed_inst[w] = tfloatinst_make(v_ms, t);
        double yaw_rate = 0.0; if (!first) {
          double dt = (double)(t - prev_t) / 1000000.0; if (dt > 0) {
            if (rows[j].has_rot) yaw_rate = rows[j].rot_deg_min * (M_PI/180.0)/60.0;
            else { double cog = rows[j].cog_deg * M_PI/180.0; cog = unwrap_angle(prev_cog, cog); yaw_rate = (cog - prev_cog)/dt; prev_cog = cog; }
          }
        } else { prev_cog = rows[j].cog_deg * M_PI/180.0; }
        double steer = (v_ms > 0.1) ? atan2(yaw_rate * p.L, v_ms) : 0.0; steer_inst[w] = tfloatinst_make(steer, t);
        prev_t = t; first = 0; w++;
      }
      if (w < 2) { pfree(origin); free(speed_inst); free(steer_inst); free(gps_inst); start = end; continue; }

      TSequence *speed = tsequence_make((const TInstant **)speed_inst, (int)w, true, true, STEP, true);
      TSequence *steer = tsequence_make((const TInstant **)steer_inst, (int)w, true, true, STEP, true);
      TSequence *gps   = tsequence_make((const TInstant **)gps_inst,   (int)w, true, true, LINEAR, true);
      for (size_t k=0;k<w;k++){ pfree(speed_inst[k]); pfree(steer_inst[k]); pfree(gps_inst[k]); }
      free(speed_inst); free(steer_inst); free(gps_inst);

      Temporal *filt = tgeompoint_ekf((Temporal *)speed, (Temporal *)steer, (Temporal *)gps, origin, phi_prev, &p);
      TimestampTz prev_out_t = 0; int first_out = 1;
      for (size_t j = start; j < end; j++) {
        TimestampTz t = rows[j].t; if (!first_out && t <= prev_out_t) continue; prev_out_t = t; first_out = 0;
        char *ts = timestamptz_out(t);
        double lonr = rows[j].lon, latr = rows[j].lat; double latf=NAN, lonf=NAN;
        GSERIALIZED *gsf=NULL; if (tgeo_value_at_timestamptz(filt, t, false, &gsf)) {
          const POINT2D *pf = GSERIALIZED_POINT2D_P(gsf); xy_to_ll(lat0, lon0, pf->x, pf->y, &latf, &lonf); pfree(gsf);
        }
        printf("%s,%.6f,%.6f,%.6f,%.6f,%.3f,%.3f,", ts, lonr, latr, lonf, latf, rows[j].sog_kn, rows[j].cog_deg);
        if (rows[j].has_rot) printf("%.3f\n", rows[j].rot_deg_min); else printf("\n");
        pfree(ts);
      }
      pfree(speed); pfree(steer); pfree(gps); pfree(filt); pfree(origin);
      start = end;
    }
    meos_finalize();
    return EXIT_SUCCESS;
  }

  /* Read CSV */
  char buf[4096];
  /* Skip header if present */
  if (fgets(buf, sizeof(buf), fin) == NULL) { fprintf(stderr, "Empty input\n"); return EXIT_FAILURE; }
  if (strstr(buf, "Timestamp") == NULL) {
    /* First line is data; push back by resetting file pointer */
    fseek(fin, 0, SEEK_SET);
  }

  AisRow *rows = NULL; size_t nrows = 0, cap = 0;
  size_t lines_read = 0;
  while (fgets(buf, sizeof(buf), fin)) {
    lines_read++;
    if (progress && (lines_read % (size_t)prog_every_rows == 0)) {
      fprintf(stderr, "[ais_ekf] read %zu lines...\n", lines_read);
      fflush(stderr);
    }
    /* Make a mutable copy */
    char line[4096]; strncpy(line, buf, sizeof(line)); line[sizeof(line)-1] = '\0';
    char *f[32]; int nf = split_csv(line, f, 32);
    if (nf < 10) continue;
    /* Filter by type and MMSI */
    if (strcmp(f[1], "Class A") != 0) continue;
    if (strcmp(f[2], mmsi_filter) != 0) continue;

    TimestampTz t; if (!parse_ts_any(&t, f[0])) continue;
    if (!*f[3] || !*f[4]) continue;
    double lat = atof(f[3]);
    double lon = atof(f[4]);
    double sog_kn = (*f[7]) ? atof(f[7]) : 0.0;
    double cog_deg = (*f[8]) ? atof(f[8]) : 0.0;
    double rot = NAN; int has_rot = 0;
    if (*f[6]) { rot = atof(f[6]); has_rot = 1; }

    if (nrows == cap) { cap = cap ? cap*2 : 1024; rows = realloc(rows, cap * sizeof(*rows)); }
    rows[nrows++] = (AisRow){ .t=t, .lat=lat, .lon=lon, .sog_kn=sog_kn, .cog_deg=cog_deg, .rot_deg_min=rot, .has_rot=has_rot };
  }
  if (fin != stdin) fclose(fin);
  if (nrows < 2) { fprintf(stderr, "Not enough rows for MMSI %s\n", mmsi_filter); return EXIT_FAILURE; }

  /* Sort by time and enforce strictly increasing timestamps */
  qsort(rows, nrows, sizeof(AisRow), cmp_rows);
  size_t u = 0;
  for (size_t i = 0; i < nrows; i++) {
    if (u == 0 || rows[i].t > rows[u-1].t) rows[u++] = rows[i];
    /* else skip duplicates/non-increasing */
  }
  nrows = u;
  if (nrows < 2) { fprintf(stderr, "Not enough distinct timestamps for MMSI %s\n", mmsi_filter); return EXIT_FAILURE; }

  /* Build temporals */
  TInstant **speed_inst = malloc(nrows * sizeof(*speed_inst));
  TInstant **steer_inst = malloc(nrows * sizeof(*steer_inst));
  TInstant **gps_inst   = malloc(nrows * sizeof(*gps_inst));

  /* Origin and initial heading */
  double lat0 = rows[0].lat, lon0 = rows[0].lon;
  double x0 = 0.0, y0 = 0.0;
  GSERIALIZED *origin = geompoint_make2d(0 /*SRID unknown*/, x0, y0);
  double phi_prev = rows[0].cog_deg * M_PI/180.0;

  /* Vehicle geometry (defaults) */
  TGeoEkfParams p = {
    .L=2.83, .H=0.76, .A=3.78, .B=0.5,
    .q_x=0.5, .q_y=0.5, .q_phi=0.5,
    .r_x=0.5, .r_y=5.0,
    .p0_x=10.0, .p0_y=10.0, .p0_phi=0.5,
    .epsS=0.0, .hold_last_controls=true
  };

  TimestampTz prev_t = 0; double prev_cog = 0.0; int first = 1; size_t w = 0;
  for (size_t i=0; i<nrows; i++) {
    TimestampTz t = rows[i].t;
    if (!first && t <= prev_t) continue;
    double x, y; ll_to_xy(lat0, lon0, rows[i].lat, rows[i].lon, &x, &y);
    GSERIALIZED *pt = geompoint_make2d(0 /*SRID unknown*/, x, y);
    gps_inst[w] = tpointinst_make(pt, t); pfree(pt);

    double v_ms = rows[i].sog_kn * KNOTS_TO_MS;
    speed_inst[w] = tfloatinst_make(v_ms, t);

    double yaw_rate = 0.0; if (!first) {
      double dt = (double)(t - prev_t) / 1000000.0; if (dt > 0) {
        if (rows[i].has_rot) yaw_rate = rows[i].rot_deg_min * (M_PI/180.0) / 60.0;
        else { double cog = rows[i].cog_deg * M_PI/180.0; cog = unwrap_angle(prev_cog, cog); yaw_rate = (cog - prev_cog) / dt; prev_cog = cog; }
      }
    } else { prev_cog = rows[i].cog_deg * M_PI/180.0; }
    double steer = (v_ms > 0.1) ? atan2(yaw_rate * p.L, v_ms) : 0.0;
    steer_inst[w] = tfloatinst_make(steer, t);
    prev_t = t; first = 0; w++;
  }

  /* Build sequences (step for controls, linear for GPS) */
  TSequence *speed = tsequence_make((const TInstant **)speed_inst, (int)w, true, true, STEP, true);
  TSequence *steer = tsequence_make((const TInstant **)steer_inst, (int)w, true, true, STEP, true);
  TSequence *gps   = tsequence_make((const TInstant **)gps_inst,   (int)w, true, true, LINEAR, true);

  for (size_t i=0;i<w;i++) { pfree(speed_inst[i]); pfree(steer_inst[i]); pfree(gps_inst[i]); }
  free(speed_inst); free(steer_inst); free(gps_inst);

  /* Run EKF */
  Temporal *filt = tgeompoint_ekf((Temporal *)speed, (Temporal *)steer, (Temporal *)gps, origin, phi_prev, &p);
  /* Output EPSG:4326 (WGS84): lon/lat for raw and filtered */
  printf("t,lon_raw,lat_raw,lon_filt,lat_filt,sog_kn,cog_deg,rot_deg_min\n");

  for (size_t i = 0; i < nrows; i++) {
    char *ts = timestamptz_out(rows[i].t);
    double lonr = rows[i].lon, latr = rows[i].lat;
    GSERIALIZED *gsf = NULL; double latf = NAN, lonf = NAN;
    if (tgeo_value_at_timestamptz(filt, rows[i].t, false, &gsf)) {
      const POINT2D *pf = GSERIALIZED_POINT2D_P(gsf);
      xy_to_ll(lat0, lon0, pf->x, pf->y, &latf, &lonf);
      pfree(gsf);
    }
    printf("%s,%.6f,%.6f,%.6f,%.6f,%.3f,%.3f,", ts, lonr, latr, lonf, latf,
           rows[i].sog_kn, rows[i].cog_deg);
    if (rows[i].has_rot) printf("%.3f\n", rows[i].rot_deg_min); else printf("\n");
    pfree(ts);
  }
  pfree(speed); pfree(steer); pfree(gps); pfree(filt);
  pfree(origin);
  free(rows);

  meos_finalize();
  return EXIT_SUCCESS;
}
