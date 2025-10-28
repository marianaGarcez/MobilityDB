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
 * @brief Extended Kalman Filter (EKF) API for MEOS temporal cleaning (internal)
 */

#ifndef __TEKF_H__
#define __TEKF_H__

#include <stdbool.h>
#include <meos.h>
#include "temporal/meos_catalog.h" /* for meosType */

#if defined(MEOS_EXPERIMENTAL_ANALYTICS)

typedef struct TEkfModel
{
  int N;  /* state dim */
  int M;  /* measurement dim */

  bool (*f)(const double *x, const double *u, double dt,
            double *fx, double *F, void *ctx);
  bool (*h)(const double *x, double *hx, double *H, void *ctx);
  bool (*Q)(double dt, double *Q, void *ctx); /* optional */
  bool (*R)(double *R, void *ctx);            /* optional */

  bool (*z_from_value)(Datum value, meosType temptype, double *z, void *ctx);
  bool (*value_from_state)(const double *x, meosType temptype, Datum *out_value, void *ctx);
} TEkfModel;

typedef struct TEkfParams
{
  double default_dt;
  double gate_sigma;
  bool   fill_estimates;
  const double *P0_diag; /* len N */
  const double *x0;      /* len N */
  const double *Q_diag;  /* len N (fallback if model->Q is NULL) */
  const double *R_diag;  /* len M (fallback if model->R is NULL) */
} TEkfParams;

/* CV model context and builder */
typedef struct TEkfCvCtx {
  int D;
  double q_accel_var;
  double r_meas_var;
} TEkfCvCtx;

bool tekf_make_cv_model(int D, TEkfModel *model);

/* GPS-like range model context and builder */
typedef struct TEkfGpsCtx
{
  int ndim;           /* 2 or 3 */
  int M;              /* number of anchors */
  const double *anchors; /* length M*ndim */
  bool use_bias;      /* include bias in state */
  double q_accel_var; /* (m/s^2)^2 */
  double q_bias_var;  /* bias RW variance */
  double r_meas_var;  /* measurement variance */
} TEkfGpsCtx;

bool tekf_make_gps_model(const TEkfGpsCtx *ctx, TEkfModel *model);

/* Entrypoint */
Temporal * temporal_ekf_clean(const Temporal *temp, const TEkfModel *model,
                              const TEkfParams *params, void *ctx, int *removed_count);

#endif /* MEOS_EXPERIMENTAL_ANALYTICS */

#endif /* __TEKF_H__ */
