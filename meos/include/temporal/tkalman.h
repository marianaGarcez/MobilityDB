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
 * @brief Temporal Kalman cleaning (experimental) for 1D/2D(/3D) numeric series
 *
 * This header declares a lightweight, internal MEOS API that performs
 * outlier cleaning of temporal numeric payloads using a constant-velocity
 * Linear Kalman Filter (LKF) model applied independently per axis. It is
 * intended for library-internal and MEOS clients; it is not exported to the
 * SQL extension layer. The function prototype is gated by the CMake option
 * `MEOS_EXPERIMENTAL_ANALYTICS`.
 *
 * Outlier definition:
 * - Innovation: \f$\nu = z - H x^-\f$; innovation covariance:
 *   \f$S = H P^- H^T + R\f$.
 * - Gate (D=1): reject when `|\f$\nu\f$| / sqrt(S) > gate_sigma` (z-score).
 * - Gate (D>1): reject when Mahalanobis distance \f$\nu^T S^{-1} \nu > gate_\sigma^2\f$.
 * - If `fill_estimates` is true, outliers are replaced by the predicted/filtered
 *   estimate; otherwise they are dropped from the output.
 */

#ifndef __TKALMAN_H__
#define __TKALMAN_H__

/* MEOS */
#include <meos.h>

/*****************************************************************************
 * Parameters
 *****************************************************************************/

/**
 * @brief Parameters for temporal Kalman cleaning
 *
 * The filter assumes a constant-velocity (CV) process model per axis. When
 * input sampling is irregular, a default time step can be specified via
 * @ref default_dt.
 *
 * Units and meanings:
 * - @ref default_dt: seconds; used when observed Δt ≤ 0.
 * - @ref q_accel_var: process noise variance on acceleration `(m/s^2)^2` per axis.
 * - @ref r_meas_var: measurement noise variance per axis (same units as value).
 * - @ref init_pos_var: initial covariance P0 for position (diagonal).
 * - @ref init_vel_var: initial covariance P0 for velocity (diagonal).
 * - @ref gate_sigma: gating threshold in standard deviations.
 * - @ref fill_estimates: if true, emit predicted/filtered estimate instead of dropping outliers.
 */
typedef struct {
  double default_dt;          /**< Default Δt in seconds when sampling is irregular or Δt=0 */
  double q_accel_var;         /**< Process noise variance (m/s^2)^2 applied per axis */
  double r_meas_var;          /**< Measurement noise variance per axis */
  double init_pos_var;        /**< Initial covariance for position (P0) */
  double init_vel_var;        /**< Initial covariance for velocity (P0) */
  double gate_sigma;          /**< N-sigma gating threshold for innovation rejection (e.g., 3.5) */
  bool   fill_estimates;      /**< If true, emit EKF estimate when rejecting outliers; else drop */
} TKalmanParams;

/*****************************************************************************
 * API
 *****************************************************************************/

/**
 * @brief Clean outliers in a temporal numeric series using a CV Kalman filter
 *
 * Supported inputs (others are currently ignored/pass-through):
 * - Subtypes: `TINSTANT`, `TSEQUENCE`, `TSEQUENCESET`.
 * - Interpolation: `DISCRETE`, `STEP`, `LINEAR`.
 * - Basetypes:
 *   - 1D numeric: `tnumber` with base `double` (aka float8).
 *   - Numeric vectors: `tdouble2`, `tdouble3`, `tdouble4`.
 *     Initial implementation focuses on 1D/2D; 3D support is optional.
 *
 * Semantics:
 * - Applies per-axis constant-velocity filtering with configurable process
 *   and measurement noise. For irregular sampling, uses observed Δt when
 *   positive; otherwise falls back to `params->default_dt`.
 * - Outlier gating uses `params->gate_sigma` on the innovation magnitude.
 * - When an observation is rejected and `fill_estimates` is true, an
 *   estimated value at the same timestamp is emitted; otherwise the instant
 *   is dropped and `removed_count` is incremented.
 * - Unsupported basetypes are returned unchanged and `removed_count` is set
 *   to 0.
 *
 * @param[in] temp Temporal value to clean
 * @param[in] params Kalman filter parameters; when NULL, uses @ref tkalman_default_params
 * @param[out] removed_count Optional pointer to receive number of removed
 *             observations; set to 0 if no removals or on pass-through
 * @return A newly allocated cleaned temporal value; NULL on error
 */
#if defined(MEOS_EXPERIMENTAL_ANALYTICS)
/**
 * @brief Return sensible defaults for meter/second scales
 *
 * Defaults:
 * - gate_sigma = 3.5
 * - q_accel_var = 0.5^2
 * - r_meas_var = 1.0^2
 * - init_pos_var = 25.0
 * - init_vel_var = 4.0
 * - default_dt = 1.0
 * - fill_estimates = false
 */
extern TKalmanParams tkalman_default_params(void);

extern Temporal *
temporal_kalman_clean(const Temporal *temp, const TKalmanParams *params,
                      int *removed_count);
#endif /* MEOS_EXPERIMENTAL_ANALYTICS */

/*****************************************************************************
 * End
 *****************************************************************************/

#endif /* __TKALMAN_H__ */
