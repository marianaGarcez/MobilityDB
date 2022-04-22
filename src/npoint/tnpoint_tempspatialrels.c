/*****************************************************************************
 *
 * This MobilityDB code is provided under The PostgreSQL License.
 * Copyright (c) 2016-2022, Université libre de Bruxelles and MobilityDB
 * contributors
 *
 * MobilityDB includes portions of PostGIS version 3 source code released
 * under the GNU General Public License (GPLv2 or later).
 * Copyright (c) 2001-2022, PostGIS contributors
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
 * @file tnpoint_tempspatialrels.c
 * @brief Temporal spatial relationships for temporal network points.
 *
 * These relationships are applied at each instant and result in a temporal
 * boolean/text. The following relationships are supported:
 * tcontains, tdisjoint, tintersects, ttouches, and tdwithin
 */

#include "npoint/tnpoint_tempspatialrels.h"

/* MobilityDB */
#include "general/temporal_util.h"
#include "point/tpoint_spatialfuncs.h"
#include "point/tpoint_tempspatialrels.h"
#include "npoint/tnpoint_spatialfuncs.h"

/*****************************************************************************
 * Generic functions
 *****************************************************************************/

/**
 * @ingroup libmeos_temporal_spatial_rel
 * @brief Return the temporal disjoint/intersection relationship between the
 * temporal network point and the network point
 */
Temporal *
tinterrel_tnpoint_npoint(const Temporal *temp, const Npoint *np, bool tinter,
  bool restr, Datum atvalue)
{
  ensure_same_srid(tnpoint_srid(temp), npoint_srid(np));
  Temporal *tempgeom = tnpoint_tgeompoint(temp);
  GSERIALIZED *geo = (GSERIALIZED *) DatumGetPointer(npoint_geom(np));
  /* Result depends on whether we are computing tintersects or tdisjoint */
  Temporal *result = tinterrel_tpoint_geo(tempgeom, geo, tinter,
    restr, atvalue);
  pfree(tempgeom);
  pfree(geo);
  return result;
}

/**
 * @ingroup libmeos_temporal_spatial_rel
 * @brief Return the temporal disjoint/intersection relationship between the
 * temporal network point and the geometry
 */
Temporal *
tinterrel_tnpoint_geo(const Temporal *temp, const GSERIALIZED *geo, bool tinter,
  bool restr, Datum atvalue)
{
  if (gserialized_is_empty(geo))
    return NULL;
  ensure_same_srid(tnpoint_srid(temp), gserialized_get_srid(geo));
  Temporal *tempgeom = tnpoint_tgeompoint(temp);
  /* Result depends on whether we are computing tintersects or tdisjoint */
  Temporal *result = tinterrel_tpoint_geo(tempgeom, geo, tinter, restr,
    atvalue);
  pfree(tempgeom);
  return result;
}

/*****************************************************************************/

/**
 * @ingroup libmeos_temporal_spatial_rel
 * @brief Return the temporal contains relationship between the geometry and
 * the temporal network point
 */
Temporal *
tcontains_geo_tnpoint(GSERIALIZED *geo, Temporal *temp, bool restr,
  Datum atvalue)
{
  if (gserialized_is_empty(geo))
    return NULL;
  Temporal *tempgeom = tnpoint_tgeompoint(temp);
  Temporal *result = tcontains_geo_tpoint(geo, tempgeom, restr, atvalue);
  pfree(tempgeom);
  return result;
}

/**
 * @ingroup libmeos_temporal_spatial_rel
 * @brief Return the temporal touches relationship between the temporal network
 * point and the geometry
 */
Temporal *
ttouches_tnpoint_geo(const Temporal *temp, const GSERIALIZED *geo, bool restr,
  Datum atvalue)
{
  if (gserialized_is_empty(geo))
    return NULL;
  ensure_same_srid(tnpoint_srid(temp), gserialized_get_srid(geo));
  Temporal *tempgeom = tnpoint_tgeompoint(temp);
  /* Result depends on whether we are computing tintersects or tdisjoint */
  Temporal *result = ttouches_tpoint_geo(tempgeom, geo, restr, atvalue);
  pfree(tempgeom);
  return result;
}

/**
 * @ingroup libmeos_temporal_spatial_rel
 * @brief Return the temporal touches relationship between the temporal network
 * point and the geometry
 */
Temporal *
ttouches_geo_tnpoint(const GSERIALIZED *geo, const Temporal *temp, bool restr,
  Datum atvalue)
{
  return ttouches_tnpoint_geo(temp, geo, restr, atvalue);
}

/**
 * @ingroup libmeos_temporal_spatial_rel
 * @brief Return the temporal touches relationship between the temporal network
 * point and the network point
 */
Temporal *
ttouches_tnpoint_npoint(const Temporal *temp, const Npoint *np, bool restr,
  Datum atvalue)
{
  ensure_same_srid(tnpoint_srid(temp), npoint_srid(np));
  Temporal *tempgeom = tnpoint_tgeompoint(temp);
  GSERIALIZED *geo = (GSERIALIZED *) DatumGetPointer(npoint_geom(np));
  /* Result depends on whether we are computing tintersects or tdisjoint */
  Temporal *result = ttouches_tpoint_geo(tempgeom, geo, restr, atvalue);
  pfree(tempgeom);
  pfree(geo);
  return result;
}

/**
 * @ingroup libmeos_temporal_spatial_rel
 * @brief Return the temporal touches relationship between the temporal network
 * point and the network point
 */
Temporal *
ttouches_npoint_tnpoint(const Npoint *np, const Temporal *temp, bool restr,
  Datum atvalue)
{
  return ttouches_tnpoint_npoint(temp, np, restr, atvalue);
}

/**
 * @ingroup libmeos_temporal_spatial_rel
 * @brief Return a temporal Boolean that states whether the geometry and the
 * temporal network point are within the given distance
 */
Temporal *
tdwithin_tnpoint_geo(Temporal *temp, GSERIALIZED *geo, Datum dist, bool restr,
  Datum atvalue)
{
  if (gserialized_is_empty(geo))
    return NULL;
  Temporal *tempgeom = tnpoint_tgeompoint(temp);
  Temporal *result = tdwithin_tpoint_geo(tempgeom, geo, dist, restr, atvalue);
  pfree(tempgeom);
  return result;
}

/**
 * @ingroup libmeos_temporal_spatial_rel
 * @brief Return a temporal Boolean that states whether the geometry and the
 * temporal network point are within the given distance
 */
Temporal *
tdwithin_geo_tnpoint(GSERIALIZED *geo, Temporal *temp, Datum dist, bool restr,
  Datum atvalue)
{
  return tdwithin_tnpoint_geo(temp, geo, dist, restr, atvalue);
}

/**
 * @ingroup libmeos_temporal_spatial_rel
 * @brief Return a temporal Boolean that states whether the network point and
 * the temporal network point are within the given distance
 */
Temporal *
tdwithin_tnpoint_npoint(Temporal *temp, Npoint *np, Datum dist, bool restr,
  Datum atvalue)
{
  Datum geom = npoint_geom(np);
  GSERIALIZED *geo = (GSERIALIZED *) PG_DETOAST_DATUM(geom);
  Temporal *tempgeom = tnpoint_tgeompoint(temp);
  Temporal *result = tdwithin_tpoint_geo(tempgeom, geo, dist, restr, atvalue);
  pfree(geo);
  return result;
}

/**
 * @ingroup libmeos_temporal_spatial_rel
 * @brief Return a temporal Boolean that states whether the network point and
 * the temporal network point are within the given distance
 */
Temporal *
tdwithin_npoint_tnpoint(Npoint *np, Temporal *temp, Datum dist, bool restr,
  Datum atvalue)
{
  return tdwithin_tnpoint_npoint(temp, np, dist, restr, atvalue);
}

/**
 * @ingroup libmeos_temporal_spatial_rel
 * @brief Return a temporal Boolean that states whether the temporal network
 * points are within the given distance
 */
Temporal *
tdwithin_tnpoint_tnpoint(Temporal *temp1, Temporal *temp2, Datum dist,
  bool restr, Datum atvalue)
{
  Temporal *sync1, *sync2;
  /* Return NULL if the temporal points do not intersect in time
   * The operation is synchronization without adding crossings */
  if (! intersection_temporal_temporal(temp1, temp2, SYNCHRONIZE_NOCROSS,
    &sync1, &sync2))
    return NULL;

  Temporal *tempgeom1 = tnpoint_tgeompoint(sync1);
  Temporal *tempgeom2 = tnpoint_tgeompoint(sync2);
  Temporal *result = tdwithin_tpoint_tpoint(tempgeom1, tempgeom2, dist,
    restr, atvalue);
  pfree(tempgeom1); pfree(tempgeom2);
  return result;
}

/*****************************************************************************/
/*****************************************************************************/
/*                        MobilityDB - PostgreSQL                            */
/*****************************************************************************/
/*****************************************************************************/

#ifndef MEOS

/*****************************************************************************
 * Generic functions
 *****************************************************************************/

/**
 * Return the temporal disjoint/intersection relationship between the temporal
 * network point and the network point
 */
static Datum
tinterrel_tnpoint_npoint_ext(FunctionCallInfo fcinfo, bool tinter)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Npoint *np = PG_GETARG_NPOINT_P(1);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 3)
  {
    atvalue = PG_GETARG_DATUM(2);
    restr = true;
  }
  Temporal *result = tinterrel_tnpoint_npoint(temp, np, tinter, restr, atvalue);
  PG_FREE_IF_COPY(temp, 0);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

/**
 * Return the temporal disjoint/intersection relationship between the network
 * point and the temporal network point
 */
static Datum
tinterrel_npoint_tnpoint_ext(FunctionCallInfo fcinfo, bool tinter)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(1);
  Npoint *np = PG_GETARG_NPOINT_P(0);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 3)
  {
    atvalue = PG_GETARG_DATUM(2);
    restr = true;
  }
  Temporal *result = tinterrel_tnpoint_npoint(temp, np, tinter, restr, atvalue);
  PG_FREE_IF_COPY(temp, 1);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

/**
 * Return the temporal disjoint/intersection relationship between the geometry
 * and the temporal network point
 */
static Datum
tinterrel_geo_tnpoint_ext(FunctionCallInfo fcinfo, bool tinter)
{
  GSERIALIZED *geo = PG_GETARG_GSERIALIZED_P(0);
  Temporal *temp = PG_GETARG_TEMPORAL_P(1);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 3)
  {
    atvalue = PG_GETARG_DATUM(2);
    restr = true;
  }
  Temporal *result = tinterrel_tnpoint_geo(temp, geo, tinter, restr, atvalue);
  PG_FREE_IF_COPY(geo, 0);
  PG_FREE_IF_COPY(temp, 1);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

/**
 * Return the temporal disjoint/intersection relationship between the temporal
 * network point and the geometry
 */
static Datum
tinterrel_tnpoint_geo_ext(FunctionCallInfo fcinfo, bool tinter)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  GSERIALIZED *geo = PG_GETARG_GSERIALIZED_P(1);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 3)
  {
    atvalue = PG_GETARG_DATUM(2);
    restr = true;
  }
  Temporal *result = tinterrel_tnpoint_geo(temp, geo, tinter, restr, atvalue);
  PG_FREE_IF_COPY(geo, 1);
  PG_FREE_IF_COPY(temp, 0);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

/*****************************************************************************
 * Temporal contains
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Tcontains_geo_tnpoint);
/**
 * Return the temporal contains relationship between the geometry and the
 * temporal network point
 */
PGDLLEXPORT Datum
Tcontains_geo_tnpoint(PG_FUNCTION_ARGS)
{
  GSERIALIZED *geo = PG_GETARG_GSERIALIZED_P(0);
  Temporal *temp = PG_GETARG_TEMPORAL_P(1);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 3)
  {
    atvalue = PG_GETARG_DATUM(2);
    restr = true;
  }
  Temporal *result = tcontains_geo_tnpoint(geo, temp, restr, atvalue);
  PG_FREE_IF_COPY(geo, 0);
  PG_FREE_IF_COPY(temp, 1);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

/*****************************************************************************
 * Temporal disjoint
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Tdisjoint_geo_tnpoint);
/**
 * Return the temporal disjoint relationship between the geometry and the
 * temporal network point
 */
PGDLLEXPORT Datum
Tdisjoint_geo_tnpoint(PG_FUNCTION_ARGS)
{
  return tinterrel_geo_tnpoint_ext(fcinfo, TDISJOINT);
}

PG_FUNCTION_INFO_V1(Tdisjoint_npoint_tnpoint);
/**
 * Return the temporal disjoint relationship between the network point and the
 * temporal network point
 */
PGDLLEXPORT Datum
Tdisjoint_npoint_tnpoint(PG_FUNCTION_ARGS)
{
  return tinterrel_npoint_tnpoint_ext(fcinfo, TDISJOINT);
}

PG_FUNCTION_INFO_V1(Tdisjoint_tnpoint_geo);
/**
 * Return the temporal disjoint relationship between the temporal network point
 * and the geometry
 */
PGDLLEXPORT Datum
Tdisjoint_tnpoint_geo(PG_FUNCTION_ARGS)
{
  return tinterrel_tnpoint_geo_ext(fcinfo, TDISJOINT);
}

PG_FUNCTION_INFO_V1(Tdisjoint_tnpoint_npoint);
/**
 * Return the temporal disjoint relationship between the temporal network point
 * and the network point
 */
PGDLLEXPORT Datum
Tdisjoint_tnpoint_npoint(PG_FUNCTION_ARGS)
{
  return tinterrel_tnpoint_npoint_ext(fcinfo, TDISJOINT);
}

/*****************************************************************************
 * Temporal intersects
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Tintersects_geo_tnpoint);
/**
 * Return the temporal intersects relationship between the geometry and the
 * temporal network point
 */
PGDLLEXPORT Datum
Tintersects_geo_tnpoint(PG_FUNCTION_ARGS)
{
  return tinterrel_geo_tnpoint_ext(fcinfo, TINTERSECTS);
}

PG_FUNCTION_INFO_V1(Tintersects_npoint_tnpoint);
/**
 * Return the temporal intersects relationship between the network point and
 * the temporal network point
 */
PGDLLEXPORT Datum
Tintersects_npoint_tnpoint(PG_FUNCTION_ARGS)
{
  return tinterrel_npoint_tnpoint_ext(fcinfo, TINTERSECTS);
}

PG_FUNCTION_INFO_V1(Tintersects_tnpoint_geo);
/**
 * Return the temporal intersects relationship between the temporal network
 * point and the geometry
 */
PGDLLEXPORT Datum
Tintersects_tnpoint_geo(PG_FUNCTION_ARGS)
{
  return tinterrel_tnpoint_geo_ext(fcinfo, TINTERSECTS);
}

PG_FUNCTION_INFO_V1(Tintersects_tnpoint_npoint);
/**
 * Return the temporal intersects relationship between the temporal network
 * point and the network point
 */
PGDLLEXPORT Datum
Tintersects_tnpoint_npoint(PG_FUNCTION_ARGS)
{
  return tinterrel_tnpoint_npoint_ext(fcinfo, TINTERSECTS);
}

/*****************************************************************************
 * Temporal touches
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Ttouches_geo_tnpoint);
/**
 * Return the temporal touches relationship between the geometry and
 * the temporal network point
 */
PGDLLEXPORT Datum
Ttouches_geo_tnpoint(PG_FUNCTION_ARGS)
{
  GSERIALIZED *geo = PG_GETARG_GSERIALIZED_P(0);
  Temporal *temp = PG_GETARG_TEMPORAL_P(1);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 3)
  {
    atvalue = PG_GETARG_DATUM(2);
    restr = true;
  }
  Temporal *result = ttouches_geo_tnpoint(geo, temp, restr, atvalue);
  PG_FREE_IF_COPY(geo, 0);
  PG_FREE_IF_COPY(temp, 1);
   if (result == NULL)
    PG_RETURN_NULL();
 PG_RETURN_POINTER(result);
}

PG_FUNCTION_INFO_V1(Ttouches_npoint_tnpoint);
/**
 * Return the temporal touches relationship between the network point and
 * the temporal network point
 */
PGDLLEXPORT Datum
Ttouches_npoint_tnpoint(PG_FUNCTION_ARGS)
{
  Npoint *np = PG_GETARG_NPOINT_P(0);
  Temporal *temp = PG_GETARG_TEMPORAL_P(1);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 3)
  {
    atvalue = PG_GETARG_DATUM(2);
    restr = true;
  }
  Temporal *result = ttouches_npoint_tnpoint(np, temp, restr, atvalue);
  PG_FREE_IF_COPY(temp, 1);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

PG_FUNCTION_INFO_V1(Ttouches_tnpoint_geo);
/**
 * Return the temporal touches relationship between the temporal network point
 * and the geometry
 */
PGDLLEXPORT Datum
Ttouches_tnpoint_geo(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  GSERIALIZED *geo = PG_GETARG_GSERIALIZED_P(1);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 3)
  {
    atvalue = PG_GETARG_DATUM(2);
    restr = true;
  }
  Temporal *result = ttouches_tnpoint_geo(temp, geo, restr, atvalue);
  PG_FREE_IF_COPY(temp, 0);
  PG_FREE_IF_COPY(geo, 1);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

PG_FUNCTION_INFO_V1(Ttouches_tnpoint_npoint);
/**
 * Return the temporal touches relationship between the temporal network point
 * and the network point
 */
PGDLLEXPORT Datum
Ttouches_tnpoint_npoint(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Npoint *np = PG_GETARG_NPOINT_P(1);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 3)
  {
    atvalue = PG_GETARG_DATUM(2);
    restr = true;
  }
  Temporal *result = ttouches_tnpoint_npoint(temp, np, restr, atvalue);
  PG_FREE_IF_COPY(temp, 0);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

/*****************************************************************************
 * Temporal dwithin
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Tdwithin_geo_tnpoint);
/**
 * Return a temporal Boolean that states whether the geometry and the
 * temporal network point are within the given distance
 */
PGDLLEXPORT Datum
Tdwithin_geo_tnpoint(PG_FUNCTION_ARGS)
{
  GSERIALIZED *geo = PG_GETARG_GSERIALIZED_P(0);
  Temporal *temp = PG_GETARG_TEMPORAL_P(1);
  Datum dist = PG_GETARG_DATUM(2);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 4)
  {
    atvalue = PG_GETARG_DATUM(3);
    restr = true;
  }
  Temporal *result = tdwithin_geo_tnpoint(geo, temp, dist, restr, atvalue);
  PG_FREE_IF_COPY(geo, 0);
  PG_FREE_IF_COPY(temp, 1);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

PG_FUNCTION_INFO_V1(Tdwithin_tnpoint_geo);
/**
 * Return a temporal Boolean that states whether the temporal network point
 * and the geometry are within the given distance
 */
PGDLLEXPORT Datum
Tdwithin_tnpoint_geo(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  GSERIALIZED *geo = PG_GETARG_GSERIALIZED_P(1);
  Datum dist = PG_GETARG_DATUM(2);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 4)
  {
    atvalue = PG_GETARG_DATUM(3);
    restr = true;
  }
  Temporal *result = tdwithin_tnpoint_geo(temp, geo, dist, restr, atvalue);
  PG_FREE_IF_COPY(temp, 0);
  PG_FREE_IF_COPY(geo, 1);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

PG_FUNCTION_INFO_V1(Tdwithin_npoint_tnpoint);
/**
 * Return a temporal Boolean that states whether the network point and the
 * temporal network point are within the given distance
 */
PGDLLEXPORT Datum
Tdwithin_npoint_tnpoint(PG_FUNCTION_ARGS)
{
  Npoint *np  = PG_GETARG_NPOINT_P(0);
  Temporal *temp = PG_GETARG_TEMPORAL_P(1);
  Datum dist = PG_GETARG_DATUM(2);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 4)
  {
    atvalue = PG_GETARG_DATUM(3);
    restr = true;
  }
  Temporal *result = tdwithin_npoint_tnpoint(np, temp, dist, restr, atvalue);
  PG_FREE_IF_COPY(temp, 1);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

PG_FUNCTION_INFO_V1(Tdwithin_tnpoint_npoint);
/**
 * Return a temporal Boolean that states whether the temporal network point
 * and the network point are within the given distance
 */
PGDLLEXPORT Datum
Tdwithin_tnpoint_npoint(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Npoint *np  = PG_GETARG_NPOINT_P(1);
  Datum dist = PG_GETARG_DATUM(2);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 4)
  {
    atvalue = PG_GETARG_DATUM(3);
    restr = true;
  }
  Temporal *result = tdwithin_tnpoint_npoint(temp, np, dist, restr, atvalue);
  PG_FREE_IF_COPY(temp, 0);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

PG_FUNCTION_INFO_V1(Tdwithin_tnpoint_tnpoint);
/**
 * Return a temporal Boolean that states whether the temporal network points
 * are within the given distance
 */
PGDLLEXPORT Datum
Tdwithin_tnpoint_tnpoint(PG_FUNCTION_ARGS)
{
  Temporal *temp1 = PG_GETARG_TEMPORAL_P(0);
  Temporal *temp2 = PG_GETARG_TEMPORAL_P(1);
  Datum dist = PG_GETARG_DATUM(2);
  bool restr = false;
  Datum atvalue = (Datum) NULL;
  if (PG_NARGS() == 4)
  {
    atvalue = PG_GETARG_DATUM(3);
    restr = true;
  }
  Temporal *result = tdwithin_tnpoint_tnpoint(temp1, temp2, dist, restr,
    atvalue);
  PG_FREE_IF_COPY(temp1, 0);
  PG_FREE_IF_COPY(temp2, 1);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

#endif /* #ifndef MEOS */

/*****************************************************************************/
