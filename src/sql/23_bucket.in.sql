/*****************************************************************************
 *
 * This MobilityDB code is provided under The PostgreSQL License.
 *
 * Copyright (c) 2016-2021, Université libre de Bruxelles and MobilityDB
 * contributors
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

/*
 * Time bucket functions for temporal types.
 * The time bucket function for timestamps is borrowed from TimescaleDB.
 * This function is licensed under the Apache License 2.0.
 * The time bucket function for temporal types generalizes this idea.
 */

-- bucketing of timestamptz happens at UTC time
CREATE OR REPLACE FUNCTION timeBucket(ts TIMESTAMPTZ, bucket_width interval)
  RETURNS TIMESTAMPTZ
  AS 'MODULE_PATHNAME', 'timestamptz_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;

CREATE OR REPLACE FUNCTION timeBucket(ts TIMESTAMPTZ, bucket_width interval, origin TIMESTAMPTZ)
  RETURNS TIMESTAMPTZ
  AS 'MODULE_PATHNAME', 'timestamptz_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;

-- If an interval is given as the third argument, the bucket alignment is offset by the interval.
CREATE OR REPLACE FUNCTION timeBucket(ts TIMESTAMPTZ, bucket_width interval, "offset" interval)
  RETURNS TIMESTAMPTZ
  LANGUAGE SQL IMMUTABLE PARALLEL SAFE STRICT AS
$BODY$
    SELECT @extschema@.timeBucket(ts-"offset", bucket_width)+"offset";
$BODY$;

CREATE OR REPLACE FUNCTION timeBucket(tbool, bucket_width interval)
  RETURNS tbool[]
  AS 'MODULE_PATHNAME', 'temporal_time_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;
CREATE OR REPLACE FUNCTION timeBucket(tint, bucket_width interval)
  RETURNS tint[]
  AS 'MODULE_PATHNAME', 'temporal_time_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;
CREATE OR REPLACE FUNCTION timeBucket(tfloat, bucket_width interval)
  RETURNS tfloat[]
  AS 'MODULE_PATHNAME', 'temporal_time_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;
CREATE OR REPLACE FUNCTION timeBucket(ttext, bucket_width interval)
  RETURNS ttext[]
  AS 'MODULE_PATHNAME', 'temporal_time_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;

CREATE OR REPLACE FUNCTION timeBucket(tbool, bucket_width interval, origin TIMESTAMPTZ)
  RETURNS tbool[]
  AS 'MODULE_PATHNAME', 'temporal_time_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;
CREATE OR REPLACE FUNCTION timeBucket(tint, bucket_width interval, origin TIMESTAMPTZ)
  RETURNS tint[]
  AS 'MODULE_PATHNAME', 'temporal_time_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;
CREATE OR REPLACE FUNCTION timeBucket(tfloat, bucket_width interval, origin TIMESTAMPTZ)
  RETURNS tfloat[]
  AS 'MODULE_PATHNAME', 'temporal_time_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;
CREATE OR REPLACE FUNCTION timeBucket(ttext, bucket_width interval, origin TIMESTAMPTZ)
  RETURNS ttext[]
  AS 'MODULE_PATHNAME', 'temporal_time_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;

/*****************************************************************************/

CREATE OR REPLACE FUNCTION rangeBucket(tint, bucket_width int)
  RETURNS tint[]
  AS 'MODULE_PATHNAME', 'tint_range_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;
CREATE OR REPLACE FUNCTION rangeBucket(tfloat, bucket_width float)
  RETURNS tfloat[]
  AS 'MODULE_PATHNAME', 'tfloat_range_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;

CREATE OR REPLACE FUNCTION rangeBucket(tint, bucket_width int, origin int)
  RETURNS tint[]
  AS 'MODULE_PATHNAME', 'tint_range_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;
CREATE OR REPLACE FUNCTION rangeBucket(tfloat, bucket_width float, origin float)
  RETURNS tfloat[]
  AS 'MODULE_PATHNAME', 'tfloat_range_bucket'
  LANGUAGE C IMMUTABLE PARALLEL SAFE STRICT;

/*****************************************************************************/
