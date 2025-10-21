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

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <MMSI> [csv_file]\n", argv[0]);
    return EXIT_FAILURE;
  }
  const char *mmsi_filter = argv[1];
  FILE *fin = (argc >= 3) ? fopen(argv[2], "r") : stdin;
  if (!fin) { perror("fopen"); return EXIT_FAILURE; }

  meos_initialize();

  /* Read CSV */
  char buf[4096];
  /* Skip header if present */
  if (fgets(buf, sizeof(buf), fin) == NULL) { fprintf(stderr, "Empty input\n"); return EXIT_FAILURE; }
  if (strstr(buf, "Timestamp") == NULL) {
    /* First line is data; push back by resetting file pointer */
    fseek(fin, 0, SEEK_SET);
  }

  AisRow *rows = NULL; size_t nrows = 0, cap = 0;
  while (fgets(buf, sizeof(buf), fin)) {
    /* Make a mutable copy */
    char line[4096]; strncpy(line, buf, sizeof(line)); line[sizeof(line)-1] = '\0';
    char *f[32]; int nf = split_csv(line, f, 32);
    if (nf < 10) continue;
    /* Filter by type and MMSI */
    if (strcmp(f[1], "Class A") != 0) continue;
    if (strcmp(f[2], mmsi_filter) != 0) continue;

    char iso[64]; if (!parse_ddmmyyyy(iso, sizeof(iso), f[0])) continue;
    TimestampTz t = timestamp_in(iso, -1);
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

  for (size_t i=0; i<nrows; i++) {
    double x, y; ll_to_xy(lat0, lon0, rows[i].lat, rows[i].lon, &x, &y);
    GSERIALIZED *pt = geompoint_make2d(0 /*SRID unknown*/, x, y);
    gps_inst[i] = tpointinst_make(pt, rows[i].t);
    pfree(pt);

    double v_ms = rows[i].sog_kn * KNOTS_TO_MS;
    speed_inst[i] = tfloatinst_make(v_ms, rows[i].t);

    /* Derive steering from yaw-rate (ROT if present else diff of COG) */
    double yaw_rate = 0.0; /* rad/s */
    if (i > 0) {
      double dt = (double)(rows[i].t - rows[i-1].t) / 1000000.0;
      if (dt > 0) {
        if (rows[i].has_rot) {
          yaw_rate = rows[i].rot_deg_min * (M_PI/180.0) / 60.0;
        } else {
          double cog = rows[i].cog_deg * M_PI/180.0;
          double cog_prev = rows[i-1].cog_deg * M_PI/180.0;
          cog = unwrap_angle(cog_prev, cog);
          yaw_rate = (cog - cog_prev) / dt;
        }
      }
    }
    double steer = 0.0;
    if (v_ms > 0.1) steer = atan2(yaw_rate * p.L, v_ms);
    steer_inst[i] = tfloatinst_make(steer, rows[i].t);
  }

  /* Build sequences (step for controls, linear for GPS) */
  TSequence *speed = tsequence_make((const TInstant **)speed_inst, (int)nrows, true, true, STEP, true);
  TSequence *steer = tsequence_make((const TInstant **)steer_inst, (int)nrows, true, true, STEP, true);
  TSequence *gps   = tsequence_make((const TInstant **)gps_inst,   (int)nrows, true, true, LINEAR, true);

  for (size_t i=0;i<nrows;i++) { pfree(speed_inst[i]); pfree(steer_inst[i]); pfree(gps_inst[i]); }
  free(speed_inst); free(steer_inst); free(gps_inst);

  /* Run EKF */
  Temporal *filt = tgeompoint_ekf((Temporal *)speed, (Temporal *)steer, (Temporal *)gps, origin, phi_prev, &p);

  char *filt_txt = tspatial_as_text(filt, 3);
  printf("%s\n", filt_txt);

  pfree(filt_txt);
  pfree(speed); pfree(steer); pfree(gps); pfree(filt);
  pfree(origin);
  free(rows);

  meos_finalize();
  return EXIT_SUCCESS;
}
