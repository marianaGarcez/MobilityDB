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
 * using speed and steering temporal controls and optional GPS tgeompoint
 * measurements. Produces a filtered tgeompoint.
 */

#ifndef __TGEO_EKF_H__
#define __TGEO_EKF_H__

/* MEOS */
#include <meos.h>

/*****************************************************************************
 * Parameters and API
 *****************************************************************************/

/**
 * @brief EKF parameters for a bicycle/ACKERMANN-like planar model.
 * Units: meters, seconds, radians.
 */
typedef struct
{
  /* Vehicle geometry */
  double L;   /* wheelbase (m) */
  double H;   /* rear overhang (m) */
  double A;   /* CoM offset X (m) */
  double B;   /* CoM offset Y (m) */

  /* Process noise diagonal (variance) */
  double q_x;    /* sigma^2 for x */
  double q_y;    /* sigma^2 for y */
  double q_phi;  /* sigma^2 for heading */

  /* Measurement noise diagonal (variance) for GPS x,y */
  double r_x;    /* sigma^2 for x */
  double r_y;    /* sigma^2 for y */

  /* Initial covariance diagonal (variance) */
  double p0_x;   /* initial P for x */
  double p0_y;   /* initial P for y */
  double p0_phi; /* initial P for heading */

  /* Numerical stabilizer for 2x2 innovation inverse */
  double epsS;  /* e.g., 1e-9 */

  /* Controls */
  bool hold_last_controls; /* If true, keep last v,delta until updated */
} TGeoEkfParams;

/**
 * @brief Run an EKF over speed and steering to filter a GPS tgeompoint.
 *
 * @param[in] speed   Temporal float (m/s)
 * @param[in] steer   Temporal float (radians)
 * @param[in] gps     Temporal geometry point (planar, same SRID as origin)
 * @param[in] origin  Initial 2D point (GSERIALIZED) defining x0,y0 and SRID
 * @param[in] phi0    Initial heading in radians
 * @param[in] p       Parameters (see TGeoEkfParams); may be NULL for defaults
 * @return Filtered tgeompoint as a TSequence with linear interpolation
 */
extern Temporal *tgeompoint_ekf(const Temporal *speed, const Temporal *steer,
  const Temporal *gps, const GSERIALIZED *origin, double phi0,
  const TGeoEkfParams *p);

#endif /* __TGEO_EKF_H__ */

