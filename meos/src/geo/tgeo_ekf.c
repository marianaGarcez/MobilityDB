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
 * @brief Extended Kalman filter (EKF) for planar tgeompoint trajectories
 * using speed and steering as controls and GPS tgeompoint as measurements.
 */

/* C */
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* MEOS */
#include <meos.h>
#include <meos_geo.h>
#include <meos_internal.h>
#include <meos_internal_geo.h>
#include "geo/tgeo_spatialfuncs.h"
#include "geo/tgeo_ekf.h"
#include "temporal/temporal.h"
/* Utilities (pfree_array) */
#include "temporal/type_util.h"

/* GSL */
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

/*****************************************************************************
 * Small helpers (angles, LA, timeline)
 *****************************************************************************/

static inline void
wrap_pi(double *a)
{
  while (*a >  M_PI) *a -= 2.0 * M_PI;
  while (*a <= -M_PI) *a += 2.0 * M_PI;
}

/* Wrap angle to [0, 2*pi) for Flink parity */
static inline void
wrap_mod2pi(double *a)
{
  const double two_pi = 2.0 * M_PI;
  double x = fmod(*a, two_pi);
  if (x < 0.0) x += two_pi;
  *a = x;
}

static inline void
mat33_set_identity(double M[9])
{
  M[0]=1; M[1]=0; M[2]=0;
  M[3]=0; M[4]=1; M[5]=0;
  M[6]=0; M[7]=0; M[8]=1;
}

static inline void
mat33_copy(double dst[9], const double src[9])
{
  memcpy(dst, src, sizeof(double) * 9);
}

static inline void
mat33_mul(const double A[9], const double B[9], double C[9])
{
  gsl_matrix_const_view Av = gsl_matrix_const_view_array(A, 3, 3);
  gsl_matrix_const_view Bv = gsl_matrix_const_view_array(B, 3, 3);
  gsl_matrix_view Cv = gsl_matrix_view_array(C, 3, 3);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Av.matrix, &Bv.matrix, 0.0, &Cv.matrix);
}

static inline void
mat33_add_inplace(double A[9], const double B[9])
{
  for (int i=0;i<9;i++) A[i] += B[i];
}

static inline void
mat33_symmetrize(double A[9])
{
  A[1] = A[3] = 0.5 * (A[1] + A[3]);
  A[2] = A[6] = 0.5 * (A[2] + A[6]);
  A[5] = A[7] = 0.5 * (A[5] + A[7]);
}

static inline bool
mat22_inv_safe(double S00, double S01, double S10, double S11,
               double eps, double *i00, double *i01, double *i10, double *i11)
{
  /* Ensure symmetry for numeric stability */
  double a=S00, b=0.5*(S01+S10), d=S11;
  double det = a*d - b*b;
  if (fabs(det) < eps)
  {
    a += eps; d += eps; det = a*d - b*b;
    if (fabs(det) < eps)
      return false;
  }
  double invdet = 1.0 / det;
  *i00 =  d * invdet;
  *i01 = -b * invdet;
  *i10 = -b * invdet;
  *i11 =  a * invdet;
  return true;
}

static int
cmp_timestamptz(const void *a, const void *b)
{
  const TimestampTz *x = (const TimestampTz *)a;
  const TimestampTz *y = (const TimestampTz *)b;
  if (*x < *y) return -1;
  if (*x > *y) return 1;
  return 0;
}

/*****************************************************************************
 * EKF core
 *****************************************************************************/

static void
ekf_predict(double mu[3], double P[9], double dt, double v_in, double steer,
            const TGeoEkfParams *p)
{
  if (dt <= 0.0)
    return;

  /*
   * Flink EKF variable mapping (aliases kept for readability/parity):
   *  - vehicleL/H/A/B  -> p->L/H/A/B
   *  - x_prev,y_prev,phi_prev -> mu[0],mu[1],mu[2] (pre-predict)
   *  - input_speed -> v_in, steering -> steer, timedif -> dt
   *  - speed (mapped) -> v
   *  - x_inc,y_inc,phi_inc -> motion increments
   *  - estimatedPoseVector -> predicted state
   *  - jacobianMatrixGt -> motion Jacobian G
   *  - movementErrorMatrixRt -> process noise Q (diagonal from p->q_*)
   *  - estimatedSigma -> predicted covariance
   */

  const double vehicleL = p->L, vehicleH = p->H, vehicleA = p->A, vehicleB = p->B;
  const double x_prev = mu[0];
  const double y_prev = mu[1];
  const double phi_prev = mu[2];

  const double t = tan(steer);
  const double v = v_in / (1.0 - t * vehicleH / vehicleL);
  const double c = cos(phi_prev);
  const double s = sin(phi_prev);

  /* Motion increments (x_inc, y_inc, phi_inc) */
  const double x_inc  = dt * ( v*c - (v/vehicleL) * t * ( vehicleA*s + vehicleB*c ) );
  const double y_inc  = dt * ( v*s + (v/vehicleL) * t * ( vehicleA*c - vehicleB*s ) );
  const double phi_inc = dt * ( v/vehicleL ) * t;

  /* Predicted state (estimatedPoseVector) */
  double estimatedPoseVector[3] = { x_prev + x_inc, y_prev + y_inc, phi_prev + phi_inc };
  mu[0] = estimatedPoseVector[0];
  mu[1] = estimatedPoseVector[1];
  mu[2] = estimatedPoseVector[2];
  wrap_mod2pi(&mu[2]);

  /* Jacobian (jacobianMatrixGt) */
  double jacobianMatrixGt[9];
  jacobianMatrixGt[0] = 1.0; jacobianMatrixGt[1] = 0.0; jacobianMatrixGt[2] = -dt*( v*s + (v/vehicleL)*t*( vehicleA*c - vehicleB*s ) );
  jacobianMatrixGt[3] = 0.0; jacobianMatrixGt[4] = 1.0; jacobianMatrixGt[5] =  dt*( v*c - (v/vehicleL)*t*( vehicleA*s + vehicleB*c ) );
  jacobianMatrixGt[6] = 0.0; jacobianMatrixGt[7] = 0.0; jacobianMatrixGt[8] =  1.0;

  /* P = G P G^T + Rt, where Rt (movementErrorMatrixRt) is diagonal */
  double GP[9];
  mat33_mul(jacobianMatrixGt, P, GP);
  double GPGt[9];
  double Gt_[9] = { jacobianMatrixGt[0], jacobianMatrixGt[3], jacobianMatrixGt[6],
                    jacobianMatrixGt[1], jacobianMatrixGt[4], jacobianMatrixGt[7],
                    jacobianMatrixGt[2], jacobianMatrixGt[5], jacobianMatrixGt[8] };
  mat33_mul(GP, Gt_, GPGt);

  /* movementErrorMatrixRt (Rt) */
  double movementErrorMatrixRt[9] = {0};
  movementErrorMatrixRt[0] = p->q_x; movementErrorMatrixRt[4] = p->q_y; movementErrorMatrixRt[8] = p->q_phi;

  for (int i=0;i<9;i++) P[i] = GPGt[i] + movementErrorMatrixRt[i];
  mat33_symmetrize(P);
}

static void
ekf_update_xy(double mu[3], double P[9], double zx, double zy,
              const TGeoEkfParams *p)
{
  /* Flink naming: gpsPosition, observationJacobianMatrix (H), gpsErrorMatrix (R) */
  double gpsPosition[2] = { zx, zy };
  /* observationJacobianMatrix is identity on x,y; implicit in the math below */
  /* Innovation y = z - H mu (positionPoseMatrix in Flink x/y path) */
  double positionPoseMatrix[2] = { gpsPosition[0] - mu[0], gpsPosition[1] - mu[1] };

  /* S = H P H^T + R -> top-left 2x2 of P + R (gpsErrorMatrix) */
  double S00 = P[0] + p->r_x;
  double S01 = P[1];
  double S10 = P[3];
  double S11 = P[4] + p->r_y;

  /* kalmanInverse (S^-1) using GSL; regularize with epsS if provided */
  double Sarr[4] = { S00, 0.5*(S01+S10), 0.5*(S01+S10), S11 };
  if (p->epsS > 0.0) { Sarr[0] += p->epsS; Sarr[3] += p->epsS; }
  double Sinvarr[4];
  gsl_matrix_view Sm = gsl_matrix_view_array(Sarr, 2, 2);
  gsl_matrix_view S_inv = gsl_matrix_view_array(Sinvarr, 2, 2);
  gsl_matrix *tmp = gsl_matrix_alloc(2,2);
  gsl_matrix_memcpy(tmp, &Sm.matrix);
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  int chol_status = gsl_linalg_cholesky_decomp(tmp);
  if (chol_status == 0) {
    gsl_linalg_cholesky_invert(tmp);
    gsl_matrix_memcpy(&S_inv.matrix, tmp);
  } else {
    gsl_permutation *perm = gsl_permutation_alloc(2);
    int signum = 0;
    gsl_matrix_memcpy(tmp, &Sm.matrix);
    gsl_linalg_LU_decomp(tmp, perm, &signum);
    gsl_linalg_LU_invert(tmp, perm, &S_inv.matrix);
    gsl_permutation_free(perm);
  }
  gsl_set_error_handler(old_handler);
  gsl_matrix_free(tmp);
  double i00 = Sinvarr[0], i01 = Sinvarr[1], i10 = Sinvarr[2], i11 = Sinvarr[3];

  /* K (kalmanGain) = P H^T S^-1 -> columns 0 and 1 of P times S^-1 (3x2) */
  double K00 = P[0]*i00 + P[1]*i10; /* row0 col0 */
  double K01 = P[0]*i01 + P[1]*i11; /* row0 col1 */
  double K10 = P[3]*i00 + P[4]*i10; /* row1 col0 */
  double K11 = P[3]*i01 + P[4]*i11; /* row1 col1 */
  double K20 = P[6]*i00 + P[7]*i10; /* row2 col0 */
  double K21 = P[6]*i01 + P[7]*i11; /* row2 col1 */

  /* updatedPose = mu + K * innovation */
  mu[0] += K00 * positionPoseMatrix[0] + K01 * positionPoseMatrix[1];
  mu[1] += K10 * positionPoseMatrix[0] + K11 * positionPoseMatrix[1];
  mu[2] += K20 * positionPoseMatrix[0] + K21 * positionPoseMatrix[1];
  wrap_mod2pi(&mu[2]);

  /* Flink-parity covariance update: P = P - K S K^T (use GSL BLAS) */
  /* Build K (3x2) and S (2x2) matrices */
  double Kmat[6] = { K00, K01,
                     K10, K11,
                     K20, K21 };
  /* Use the previously assembled symmetric Sarr from the inversion step */
  double KS[6];
  double KSKT[9];
  gsl_matrix_const_view Kv = gsl_matrix_const_view_array(Kmat, 3, 2);
  gsl_matrix_const_view Sv = gsl_matrix_const_view_array(Sarr, 2, 2);
  gsl_matrix_view KSv = gsl_matrix_view_array(KS, 3, 2);
  /* KS = K * S */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Kv.matrix, &Sv.matrix, 0.0, &KSv.matrix);
  /* KSKT = KS * K^T */
  gsl_matrix_const_view KSTv = gsl_matrix_const_view_array(KS, 3, 2);
  gsl_matrix_view KSKTv = gsl_matrix_view_array(KSKT, 3, 3);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &KSTv.matrix, &Kv.matrix, 0.0, &KSKTv.matrix);

  P[0] -= KSKT[0]; P[1] -= KSKT[1]; P[2] -= KSKT[2];
  P[3] -= KSKT[3]; P[4] -= KSKT[4]; P[5] -= KSKT[5];
  P[6] -= KSKT[6]; P[7] -= KSKT[7]; P[8] -= KSKT[8];
  mat33_symmetrize(P);
}

/*****************************************************************************
 * Public API
 *****************************************************************************/

Temporal *
tgeompoint_ekf(const Temporal *speed, const Temporal *steer,
  const Temporal *gps, const GSERIALIZED *origin, double phi0,
  const TGeoEkfParams *p)
{
  /* Validate arguments */
  VALIDATE_TFLOAT(speed, NULL);
  VALIDATE_TFLOAT(steer, NULL);
  VALIDATE_NOT_NULL(origin, NULL);
  if (gps)
    VALIDATE_TGEOMPOINT(gps, NULL);

  /* Determine SRID and ensure GPS SRID matches */
  int32_t srid = gserialized_get_srid(origin);
  if (gps)
  {
    int32_t srid_gps = tspatial_srid((Temporal *) gps);
    if (! ensure_same_srid(srid, srid_gps))
      return NULL;
  }

  /* Parameters (defaults if NULL) */
  TGeoEkfParams pd;
  if (p)
    memcpy(&pd, p, sizeof(TGeoEkfParams));
  else
  {
    memset(&pd, 0, sizeof(TGeoEkfParams));
    pd.L = 2.83; pd.H = 0.76; pd.A = 3.78; pd.B = 0.5;
    pd.q_x = 0.5; pd.q_y = 0.5; pd.q_phi = 0.5;
    pd.r_x = 0.5; pd.r_y = 5.0;
    pd.p0_x = 1.0; pd.p0_y = 1.0; pd.p0_phi = 0.1;
    pd.epsS = 0.0; pd.hold_last_controls = true;
  }

  /* Build merged timestamp array */
  int nv=0, nd=0, nz=0;
  TimestampTz *tv = temporal_timestamps(speed, &nv);
  TimestampTz *td = temporal_timestamps(steer, &nd);
  TimestampTz *tz = NULL;
  if (gps)
    tz = temporal_timestamps(gps, &nz);

  int nt = nv + nd + (gps ? nz : 0);
  if (nt == 0)
  { if (tv) pfree(tv); if (td) pfree(td); if (tz) pfree(tz); return NULL; }

  TimestampTz *tarr = palloc(sizeof(TimestampTz) * nt);
  int k=0;
  for (int i=0;i<nv;i++) tarr[k++] = tv[i];
  for (int i=0;i<nd;i++) tarr[k++] = td[i];
  if (gps) for (int i=0;i<nz;i++) tarr[k++] = tz[i];
  qsort(tarr, nt, sizeof(TimestampTz), cmp_timestamptz);
  /* Deduplicate */
  int nuniq = 0;
  for (int i=0;i<nt;i++)
  {
    if (i==0 || tarr[i] != tarr[i-1])
      tarr[nuniq++] = tarr[i];
  }
  /* Shrink logically */

  if (tv) pfree(tv); if (td) pfree(td); if (tz) pfree(tz);

  /* Initialize state from origin */
  const POINT2D *p2 = GSERIALIZED_POINT2D_P(origin);
  double mu[3] = { p2->x, p2->y, phi0 };
  double P[9] = {0}; P[0]=pd.p0_x; P[4]=pd.p0_y; P[8]=pd.p0_phi;

  /* Control cache */
  double v_cur = 0.0, delta_cur = 0.0;
  bool v_set = false, d_set = false;

  /* Output instants */
  TInstant **instants = palloc(sizeof(TInstant*) * nuniq);
  int ninsts = 0;

  TimestampTz tprev = tarr[0];
  for (int i=0;i<nuniq;i++)
  {
    TimestampTz t = tarr[i];
    double dt = (i==0) ? 0.0 : ((double)(t - tprev)) / 1000000.0; /* microseconds -> s */

    /* Read controls at t; if absent, keep last if configured */
    double vtmp; if (tfloat_value_at_timestamptz(speed, t, false, &vtmp)) { v_cur = vtmp; v_set = true; }
    double dtmp; if (tfloat_value_at_timestamptz(steer, t, false, &dtmp)) { delta_cur = dtmp; d_set = true; }
    if (!pd.hold_last_controls) { if (!v_set) v_cur = 0.0; if (!d_set) delta_cur = 0.0; }

    /* Predict */
    if (v_set && d_set)
      ekf_predict(mu, P, dt, v_cur, delta_cur, &pd);

    /* Update if GPS available at t */
    if (gps)
    {
      GSERIALIZED *gs = NULL;
      if (tgeo_value_at_timestamptz(gps, t, false, &gs))
      {
        const POINT2D *gpt = GSERIALIZED_POINT2D_P(gs);
        ekf_update_xy(mu, P, gpt->x, gpt->y, &pd);
        pfree(gs);
      }
    }

    /* Emit filtered point */
    GSERIALIZED *outpt = geopoint_make(mu[0], mu[1], 0.0, false, false, srid);
    TInstant *inst = tpointinst_make(outpt, t);
    instants[ninsts++] = inst;
    pfree(outpt);

    tprev = t;
  }

  /* Build sequence (linear interpolation, inclusive bounds) */
  TSequence *seq = tsequence_make((const TInstant **)instants, ninsts,
                                  true, true, LINEAR, NORMALIZE);
  pfree_array((void **)instants, ninsts);
  pfree(tarr);
  return (Temporal *)seq;
}
