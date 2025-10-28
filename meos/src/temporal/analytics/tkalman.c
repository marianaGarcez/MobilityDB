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
 * @brief Linear Kalman Filter (constant-velocity) helpers for temporal cleaning
 *
 * Model (per axis):
 * - State: [pos, vel]^T
 * - F(dt) = [[1, dt],[0,1]], H = [1,0]
 * - Q(dt) = q * [[dt^3/3, dt^2/2],[dt^2/2, dt]]
 * - R = r
 *
 * Outlier gate:
 * - z-score for D=1: |nu|/sqrt(S) > gate_sigma
 * - Mahalanobis for D>1: sum_i nu_i^2 / S_i > gate_sigma^2 (S is diagonal per axis)
 */

#include "temporal/tkalman.h"

#if defined(MEOS_EXPERIMENTAL_ANALYTICS)

/* C */
#include <math.h>
#include <string.h>
/* MEOS */
#include "temporal/temporal.h"
#include "temporal/doublen.h"
#include "temporal/tsequence.h"
#include "temporal/tsequenceset.h"
#include "temporal/type_util.h"

/*****************************************************************************
 * Small linear Kalman filter with CV model, dimension-agnostic via per-axis
 * independent 2x2 filters. This leverages the block-diagonal structure:
 *   F = [[1, dt], [0, 1]] ⊗ I_D
 *   Q = q * [[dt^3/3, dt^2/2], [dt^2/2, dt]] ⊗ I_D
 *   H = [1, 0] ⊗ I_D
 *   R = r * I_D
 * As a result, each axis is independent and can be filtered with a separate
 * 2x2 filter without cross-axis covariance terms.
 *****************************************************************************/

/* Axis state and covariance (2x1 state, 2x2 covariance) */
typedef struct
{
  double x_pos;  /* position */
  double x_vel;  /* velocity */
  double P00, P01, P10, P11; /* covariance matrix */
} TKF_Axis;

/* Guard reasonable D; double2/double3/double4 suggest D<=4 */
#define TKF_MAX_D 4

/* Effective dt: use observed if positive, else params->default_dt (>= 0) */
static inline double
tkf_effective_dt(double observed_dt, const TKalmanParams *params)
{
  double dt = observed_dt;
  if (!(dt > 0.0))
    dt = (params && params->default_dt > 0.0) ? params->default_dt : 0.0;
  return dt;
}

/* Initialize per-axis filter from first observation: x0=[pos0, 0], diag(P0) */
static inline void
tkf_init_axis(double z0, const TKalmanParams *p, TKF_Axis *axis)
{
  axis->x_pos = z0;
  axis->x_vel = 0.0;
  axis->P00 = (p) ? p->init_pos_var : 1.0;
  axis->P11 = (p) ? p->init_vel_var : 1.0;
  axis->P01 = 0.0;
  axis->P10 = 0.0;
}

/* Predict step for one axis using CV model */
static inline void
tkf_predict_axis(double dt, double q_accel_var, TKF_Axis *a)
{
  /* State transition F = [[1, dt],[0,1]] */
  const double xp = a->x_pos;
  const double xv = a->x_vel;
  a->x_pos = xp + dt * xv;
  /* a->x_vel unchanged */

  /* Covariance prediction: P' = F P F^T + Q */
  const double P00 = a->P00;
  const double P01 = a->P01;
  const double P10 = a->P10;
  const double P11 = a->P11;

  /* F P */
  const double FP00 = P00 + dt * P10;      /* [1 dt] * [P00; P10] */
  const double FP01 = P01 + dt * P11;      /* [1 dt] * [P01; P11] */
  const double FP10 = P10;                 /* [0  1] * [P00; P10] */
  const double FP11 = P11;                 /* [0  1] * [P01; P11] */

  /* P' = (F P) F^T + Q */
  const double Ft00 = 1.0, Ft01 = 0.0;
  const double Ft10 = dt,  Ft11 = 1.0;
  double Pp00 = FP00 * Ft00 + FP01 * Ft10; /* FP00*1 + FP01*dt */
  double Pp01 = FP00 * Ft01 + FP01 * Ft11; /* FP00*0 + FP01*1 */
  double Pp10 = FP10 * Ft00 + FP11 * Ft10; /* FP10*1 + FP11*dt */
  double Pp11 = FP10 * Ft01 + FP11 * Ft11; /* FP10*0 + FP11*1 */

  /* Add process noise Q = q * [[dt^3/3, dt^2/2],[dt^2/2, dt]] */
  const double dt2 = dt * dt;
  const double q11 = q_accel_var * (dt2 * dt / 3.0);
  const double q12 = q_accel_var * (dt2 / 2.0);
  const double q22 = q_accel_var * dt;
  Pp00 += q11;
  Pp01 += q12;
  Pp10 += q12;
  Pp11 += q22;

  a->P00 = Pp00; a->P01 = Pp01; a->P10 = Pp10; a->P11 = Pp11;
}

/* Compute innovation for one axis: nu = z - H x_pred with H=[1 0]; S = HPH^T+R */
static inline void
tkf_innovation_axis(double z, double r_var, const TKF_Axis *a,
                    double *nu_out, double *S_out)
{
  const double nu = z - a->x_pos;      /* z - H x with H=[1 0] */
  const double S  = a->P00 + r_var;    /* H P H^T + R = P00 + R */
  *nu_out = nu;
  *S_out = S;
}

/* Update step for one axis using scalar H=[1 0] */
static inline void
tkf_update_axis(double nu, double S, TKF_Axis *a)
{
  /* Kalman gain K = P H^T S^{-1}; with H=[1 0], K = [P00; P10]/S */
  const double invS = 1.0 / S;
  const double Kpos = a->P00 * invS;
  const double Kvel = a->P10 * invS;

  /* State update */
  a->x_pos += Kpos * nu;
  a->x_vel += Kvel * nu;

  /* Covariance update: P = (I - K H) P */
  const double P00 = a->P00;
  const double P01 = a->P01;
  const double P10 = a->P10;
  const double P11 = a->P11;
  a->P00 = P00 - Kpos * P00;   /* also equals (1-Kpos)*P00 */
  a->P01 = P01 - Kpos * P01;
  a->P10 = P10 - Kvel * P00;
  a->P11 = P11 - Kvel * P01;
}

/* Reject outlier decision: z-score for D==1, Mahalanobis for D>1 */
static bool
tkf_reject_outlier(int D, const double *nu, const double *S_diag,
                   const TKalmanParams *p)
{
  const double gate = (p && p->gate_sigma > 0.0) ? p->gate_sigma : 0.0;
  if (gate <= 0.0)
    return false; /* no gating */

  if (D <= 1)
  {
    const double std = sqrt(S_diag[0]);
    if (std <= 0.0)
      return false;
    const double z = fabs(nu[0]) / std;
    return (z > gate);
  }
  else
  {
    /* Mahalanobis distance with diagonal S: d^2 = Σ nu_i^2 / S_i */
    double d2 = 0.0;
    for (int i = 0; i < D; i++)
    {
      const double Si = S_diag[i];
      if (Si <= 0.0)
        continue;
      const double zi = nu[i];
      d2 += (zi * zi) / Si;
    }
    const double gate2 = gate * gate;
    return (d2 > gate2);
  }
}

/* Default parameters tuned for meter/second scales */
/**
 * @brief Default parameters for meter/second scales
 * @return Defaults: gate_sigma=3.5, q_accel_var=0.25, r_meas_var=1.0,
 *         init_pos_var=25.0, init_vel_var=4.0, default_dt=1.0, fill_estimates=false
 */
TKalmanParams
tkalman_default_params(void)
{
  TKalmanParams d;
  d.gate_sigma   = 3.5;
  d.q_accel_var  = 0.25;   /* 0.5^2 */
  d.r_meas_var   = 1.0;    /* 1.0^2 */
  d.init_pos_var = 25.0;   /* 5m sigma squared */
  d.init_vel_var = 4.0;    /* 2 m/s sigma squared */
  d.default_dt   = 1.0;    /* seconds */
  d.fill_estimates = false;
  return d;
}

/*****************************************************************************
 * Helpers for temporal traversal and packing
 *****************************************************************************/

/* Determine payload dimension D from basetype of the temporal */
static int
tkf_dim_from_basetype(meosType basetype)
{
  switch (basetype)
  {
    case T_FLOAT8: return 1;
    case T_DOUBLE2: return 2;
    case T_DOUBLE3: return 3;
    case T_DOUBLE4: return 4;
    default: return 0;
  }
}

/* Read measurement z[0..D-1] from Datum value for given temptype */
static bool
tkf_read_measurement(meosType temptype, Datum value, int D, double *z)
{
  meosType basetype = temptype_basetype(temptype);
  switch (basetype)
  {
    case T_FLOAT8:
    {
      double v = DatumGetFloat8(value);
      if (isnan(v)) return false;
      z[0] = v; return true;
    }
    case T_DOUBLE2:
    {
      const double2 *d = DatumGetDouble2P(value);
      if (!isfinite(d->a) || !isfinite(d->b)) return false;
      z[0] = d->a; z[1] = d->b; return true;
    }
    case T_DOUBLE3:
    {
      const double3 *d = DatumGetDouble3P(value);
      if (!isfinite(d->a) || !isfinite(d->b) || !isfinite(d->c)) return false;
      z[0] = d->a; z[1] = d->b; z[2] = d->c; return true;
    }
    case T_DOUBLE4:
    {
      const double4 *d = DatumGetDouble4P(value);
      if (!isfinite(d->a) || !isfinite(d->b) || !isfinite(d->c) || !isfinite(d->d)) return false;
      z[0] = d->a; z[1] = d->b; z[2] = d->c; z[3] = d->d; return true;
    }
    default:
      return false;
  }
}

/* Build a Datum from axis positions into a caller-provided stack buffer */
typedef union { double2 d2; double3 d3; double4 d4; } TKF_DatumBuf;

static Datum
tkf_build_datum(meosType temptype, int D, const TKF_Axis *axes, TKF_DatumBuf *buf)
{
  switch (temptype_basetype(temptype))
  {
    case T_FLOAT8:
      return Float8GetDatum(axes[0].x_pos);
    case T_DOUBLE2:
      double2_set(axes[0].x_pos, axes[1].x_pos, &buf->d2);
      return Double2PGetDatum(&buf->d2);
    case T_DOUBLE3:
      double3_set(axes[0].x_pos, axes[1].x_pos, axes[2].x_pos, &buf->d3);
      return Double3PGetDatum(&buf->d3);
    case T_DOUBLE4:
      double4_set(axes[0].x_pos, axes[1].x_pos, axes[2].x_pos, axes[3].x_pos, &buf->d4);
      return Double4PGetDatum(&buf->d4);
    default:
      return (Datum) 0;
  }
}

/* Initialize axes from first measurement vector */
static void
tkf_init_axes(int D, const double *z0, const TKalmanParams *p, TKF_Axis *axes)
{
  for (int i = 0; i < D; i++)
    tkf_init_axis(z0[i], p, &axes[i]);
}

/* Predict all axes with shared dt and noise */
static void
tkf_predict_axes(int D, double dt, const TKalmanParams *p, TKF_Axis *axes)
{
  const double q = (p ? p->q_accel_var : 0.0);
  for (int i = 0; i < D; i++)
    tkf_predict_axis(dt, q, &axes[i]);
}

/* Innovation and S per axis with shared R */
static void
tkf_innovation_axes(int D, const double *z, const TKalmanParams *p,
                    const TKF_Axis *axes, double *nu, double *S)
{
  const double r = (p ? p->r_meas_var : 0.0);
  for (int i = 0; i < D; i++)
    tkf_innovation_axis(z[i], r, &axes[i], &nu[i], &S[i]);
}

/* Update all axes */
static void
tkf_update_axes(int D, const double *nu, const double *S, TKF_Axis *axes)
{
  for (int i = 0; i < D; i++)
    tkf_update_axis(nu[i], S[i], &axes[i]);
}

/*****************************************************************************
 * Public API: temporal traversal and cleaning
 *****************************************************************************/

/* Clean a temporal instant: returns a new instant (always accepted) */
static TInstant *
tkf_clean_instant(const TInstant *inst, int D, TKF_Axis *axes,
                  const TKalmanParams *params, int *removed_count)
{
  (void)removed_count;
  /* Initialize from the observation and emit it */
  double z[TKF_MAX_D] = {0};
  if (!tkf_read_measurement(inst->temptype, tinstant_value_p(inst), D, z))
    return tinstant_copy(inst);
  tkf_init_axes(D, z, params, axes);
  TKF_DatumBuf tmp;
  Datum outv = tkf_build_datum(inst->temptype, D, axes, &tmp);
  TInstant *out = tinstant_make(outv, inst->temptype, inst->t);
  return out;
}

/* Clean a temporal sequence: returns new sequence or NULL if empty */
static TSequence *
tkf_clean_sequence(const TSequence *seq, int D, const TKalmanParams *params,
                   int *removed_count)
{
  const int n = seq->count;
  if (n <= 0)
    return NULL;

  TInstant **out_insts = palloc(sizeof(TInstant *) * n);
  int out_count = 0;

  TKF_Axis axes[TKF_MAX_D];
  bool has_state = false;
  TimestampTz prev_t = 0;

  for (int i = 0; i < n; i++)
  {
    const TInstant *inst = TSEQUENCE_INST_N(seq, i);
    double z[TKF_MAX_D] = {0};
    bool has_meas = tkf_read_measurement(inst->temptype, tinstant_value_p(inst), D, z);
    TKF_DatumBuf tmpbuf; /* reuse per-iteration */

    if (!has_state)
    {
      if (!has_meas)
      {
        /* Cannot initialize without a valid measurement; skip */
        prev_t = inst->t; /* keep time for dt to next */
        continue;
      }
      tkf_init_axes(D, z, params, axes);
      Datum outv = tkf_build_datum(inst->temptype, D, axes, &tmpbuf);
      out_insts[out_count++] = tinstant_make(outv, inst->temptype, inst->t);
      has_state = true;
      prev_t = inst->t;
      continue;
    }

    /* Have state: advance with dt and handle measurement/null/outlier */
    double dt = ((double) (inst->t - prev_t)) / 1000000.0;
    dt = tkf_effective_dt(dt, params);
    tkf_predict_axes(D, dt, params, axes);

    if (!has_meas)
    {
      if (params && params->fill_estimates)
      {
        Datum outv = tkf_build_datum(inst->temptype, D, axes, &tmpbuf);
        out_insts[out_count++] = tinstant_make(outv, inst->temptype, inst->t);
      }
      /* else drop */
      prev_t = inst->t;
      continue;
    }

    /* Measurement available: compute innovation and gate */
    double nu[TKF_MAX_D], S[TKF_MAX_D];
    tkf_innovation_axes(D, z, params, axes, nu, S);
    if (tkf_reject_outlier(D, nu, S, params))
    {
      if (params && params->fill_estimates)
      {
        Datum outv = tkf_build_datum(inst->temptype, D, axes, &tmpbuf);
        out_insts[out_count++] = tinstant_make(outv, inst->temptype, inst->t);
      }
      else
      {
        if (removed_count) (*removed_count)++;
      }
      /* no update on rejection */
      prev_t = inst->t;
      continue;
    }

    /* Accept: update and emit filtered estimate */
    tkf_update_axes(D, nu, S, axes);
    Datum outv = tkf_build_datum(inst->temptype, D, axes, &tmpbuf);
    out_insts[out_count++] = tinstant_make(outv, inst->temptype, inst->t);
    prev_t = inst->t;
  }

  if (out_count == 0)
  {
    pfree(out_insts);
    return NULL;
  }

  /* Build sequence preserving flags and bounds */
  interpType interp = MEOS_FLAGS_GET_INTERP(seq->flags);
  TSequence *out_seq = tsequence_make((const TInstant **) out_insts, out_count,
    seq->period.lower_inc, seq->period.upper_inc, interp, NORMALIZE_NO);
  pfree_array((void **) out_insts, out_count);
  return out_seq;
}

static Temporal *
tkf_clean_dispatch(const Temporal *temp, int D, const TKalmanParams *params,
                   int *removed_count)
{
  switch (temp->subtype)
  {
    case TINSTANT:
    {
      const TInstant *inst = (const TInstant *) temp;
      TInstant *out = tkf_clean_instant(inst, D, (TKF_Axis[TKF_MAX_D]){0}, params, removed_count);
      return (Temporal *) out;
    }
    case TSEQUENCE:
    {
      const TSequence *seq = (const TSequence *) temp;
      return (Temporal *) tkf_clean_sequence(seq, D, params, removed_count);
    }
    case TSEQUENCESET:
    {
      const TSequenceSet *ss = (const TSequenceSet *) temp;
      TSequence **seqs = palloc(sizeof(TSequence *) * ss->count);
      int nseqs = 0;
      for (int i = 0; i < ss->count; i++)
      {
        const TSequence *seq = TSEQUENCESET_SEQ_N(ss, i);
        TSequence *clean = tkf_clean_sequence(seq, D, params, removed_count);
        if (clean)
          seqs[nseqs++] = clean;
      }
      return (Temporal *) tsequenceset_make_free(seqs, nseqs, NORMALIZE_NO);
    }
    default:
      return NULL;
  }
}

Temporal *
temporal_kalman_clean(const Temporal *temp,
                      const TKalmanParams *params,
                      int *removed_count)
{
  VALIDATE_NOT_NULL(temp, NULL);
  if (removed_count) *removed_count = 0;

  /* Validate subtype */
  if (!(temp->subtype == TINSTANT || temp->subtype == TSEQUENCE || temp->subtype == TSEQUENCESET))
  {
    meos_error(ERROR, MEOS_ERR_INVALID_ARG,
      "Unsupported temporal subtype; expected TINSTANT/TSEQUENCE/TSEQUENCESET");
    return NULL;
  }

  /* Determine arity from basetype */
  const meosType basetype = temptype_basetype(temp->temptype);
  const int D = tkf_dim_from_basetype(basetype);
  if (!(D >= 1 && D <= TKF_MAX_D))
  {
    meos_error(ERROR, MEOS_ERR_INVALID_ARG_TYPE,
      "Unsupported temporal basetype; expected float8/double2/double3/double4");
    return NULL;
  }

  /* Use defaults when params == NULL, then validate non-negative variances */
  TKalmanParams defaults = tkalman_default_params();
  const TKalmanParams *P = params ? params : &defaults;
  if (P->q_accel_var < 0.0 || P->r_meas_var < 0.0 ||
      P->init_pos_var < 0.0 || P->init_vel_var < 0.0)
  {
    meos_error(ERROR, MEOS_ERR_INVALID_ARG_VALUE,
      "Variance parameters must be non-negative");
    return NULL;
  }

  return tkf_clean_dispatch(temp, D, P, removed_count);
}

#endif /* MEOS_EXPERIMENTAL_ANALYTICS */
