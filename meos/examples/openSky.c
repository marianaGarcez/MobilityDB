/*****************************************************************************
 *
 * This MobilityDB code is provided under The PostgreSQL License.
 * Copyright (c) 2016-2023, Université libre de Bruxelles and MobilityDB
 * contributors
 *
 * MobilityDB includes portions of PostGIS version 3 source code released
 * under the GNU General Public License (GPLv2 or later).
 * Copyright (c) 2001-2023, PostGIS contributors
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
 * @brief A program that uses the MEOS library for a demonstration of
 * the capabilities of the library.
 * Load planes trajectories into the structure
 * Create regions
 * *  Range Queries
 * 1- List the continents each plane pass
 * 2- List planes that were within a region from Regions during a period from
 * Periods. 
 * 3- List the pairs of planes that were both located within a region
 * from Regions during a period from Periods. 
 * 4- List one airport each plane were close to, if any.
 *  * Temporal Aggregate Queries
 * 5- Count the number of trips that were active during minute in the fisrt hour in May 25,
 * 2020. 
 * Distance Queries 
 * 6- list the distance for each plane 
 * 7- List the minimum distance ever between each plane
 *  * Nearest-Neighbor Query
 * 8- For each plane, list the three planes that are closest to it
 * @code
 * gcc -Wall -g -I/usr/local/include  -o openSky openSky.c -L/usr/local/lib
 * -lmeos
 * @endcode
 */

#include <assert.h>
#include <float.h>
#include <proj.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* MEOS */
#include <meos.h>
#include <meos_internal.h>
#include <time.h>

/* Number of instants to send in batch to the file  */
#define NO_INSTANTS_BATCH 500
/* Number of instants to keep when restarting a sequence */
#define NO_INSTANTS_KEEP 2
/* Maximum length in characters of a header record in the input CSV file */
#define MAX_LENGTH_HEADER 74
/* Maximum length in characters of a header record in the input CSV file */
#define MAX_LENGTH_HEADER_AIRPORTS 19
/* Maximum length in characters of a point in the input data */
#define MAX_LENGTH_POINT 148
/* Number of inserts that are sent in bulk */
#define NO_BULK_INSERT 20
/* Maximum number of trips */
#define MAX_TRIPS 60000
/* Maximum number of planes */
#define MAX_PLANES 300
/* Maximum number of airports */
#define MAX_AIRPORTS 4188
/* Represent self distance */
#define INT_MAX 2147483647

typedef struct {
  GSERIALIZED *geom;
} area_record;

typedef struct {
  Timestamp T;
  long int trj_id;
  float lat;
  float lon;
  float velocity;
  float heading;
  float vertrate;
  float baroaltitude;
  float geoaltitude;
} plane_point_record;

typedef struct {
  long int trj_id; /* Identifier of the trip */
  TSequence *trip; /* Array of Latest observations of the trip, by id */
} trip_record;

void readCreateAirports(FILE *airportsUsaFile,
                        area_record airports[MAX_AIRPORTS]) {
  char text_buffer_airport[MAX_LENGTH_HEADER_AIRPORTS];
  fscanf(airportsUsaFile, "%18s\n", text_buffer_airport);
  int i = 0;
  do {
    float lat, lon;
    int read = fscanf(airportsUsaFile, "%f,%f\n", &lat, &lon);
    if (read == 2) {
      char point_buffer_airport[MAX_LENGTH_POINT];
      sprintf(point_buffer_airport, "SRID=4326;Point(%lf %lf)", lon, lat);
      airports[i].geom = geometry_from_hexewkb(point_buffer_airport);
      i++;
    }
    if (read != 2 && !feof(airportsUsaFile))
      printf("Record with missing values ignored\n");
  } while (!feof(airportsUsaFile));
  printf("%d airports read\n", i);
}

void query1(int num_planes, trip_record planes[MAX_TRIPS],
            area_record continents[7]) {
  /* 1- List the continents each plane pass */
  printf("1- List the continents each plane pass\n");
  for (int i = 0; i < num_planes; i++) {
    if (eintersects_tpoint_geo((const Temporal *)planes[i].trip,
                               continents[0].geom))
      printf("Plane %ld passes at Africa\n", planes[i].trj_id);
    if (eintersects_tpoint_geo((const Temporal *)planes[i].trip,
                               continents[1].geom))
      printf("Plane %ld passes at Antarctica\n", planes[i].trj_id);
    if (eintersects_tpoint_geo((const Temporal *)planes[i].trip,
                               continents[2].geom))
      printf("Plane %ld passes at Asia\n", planes[i].trj_id);
    if (eintersects_tpoint_geo((const Temporal *)planes[i].trip,
                               continents[3].geom))
      printf("Plane %ld passes at Europe\n", planes[i].trj_id);
    if (eintersects_tpoint_geo((const Temporal *)planes[i].trip,
                               continents[4].geom))
      printf("Plane %ld passes at North_America\n", planes[i].trj_id);
    if (eintersects_tpoint_geo((const Temporal *)planes[i].trip,
                               continents[5].geom))
      printf("Plane %ld passes at Oceania\n", planes[i].trj_id);
    if (eintersects_tpoint_geo((const Temporal *)planes[i].trip,
                               continents[6].geom))
      printf("Plane %ld passes at South_America\n", planes[i].trj_id);
  }
}

void query2(int num_planes, trip_record planes[MAX_TRIPS],
            const STBox *continentsSTB[7]) {
  printf("2- List planes that were within a continent in May 25, 2020\n");
  for (int i = 0; i < num_planes; i++) {
    if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[0],
                        true))
      printf("Plane %ld was within Africa\n", planes[i].trj_id);
    else if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[1],
                             true))
      printf("Plane %ld was within Antarctica\n", planes[i].trj_id);
    else if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[2],
                             true))
      printf("Plane %ld was within Asia\n", planes[i].trj_id);
    else if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[3],
                             true))
      printf("Plane %ld was within Europe\n", planes[i].trj_id);
    else if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[4],
                             true))
      printf("Plane %ld was within North_America\n", planes[i].trj_id);
    else if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[5],
                             true))
      printf("Plane %ld was within Oceania\n", planes[i].trj_id);
    else if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[6],
                             true))
      printf("Plane %ld was within South_America\n", planes[i].trj_id);
    else
      printf("Plane %ld was not within any continent\n", planes[i].trj_id);
  }
}

void query3(int num_planes, trip_record planes[MAX_TRIPS],
            const STBox *continentsSTB[7]) {
  printf("3- List the pairs of planes that were both located within a region "
         "from Regions during a period from Periods\n");
  for (int i = 0, j = i + 1; i < num_planes && j < num_planes - 1; i++, j++) {
    if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[0],
                        true) &&
        tpoint_at_stbox((const Temporal *)planes[j].trip, continentsSTB[0],
                        true))
      printf("Planes %ld and %ld were both located within Africa\n",
             planes[i].trj_id, planes[j].trj_id);
    else if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[1],
                             true) &&
             tpoint_at_stbox((const Temporal *)planes[j].trip, continentsSTB[1],
                             true))
      printf("Planes %ld and %ld were both located within Antarctica\n",
             planes[i].trj_id, planes[j].trj_id);
    else if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[2],
                             true) &&
             tpoint_at_stbox((const Temporal *)planes[j].trip, continentsSTB[2],
                             true))
      printf("Planes %ld and %ld were both located within Asia\n",
             planes[i].trj_id, planes[j].trj_id);
    else if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[3],
                             true) &&
             tpoint_at_stbox((const Temporal *)planes[j].trip, continentsSTB[3],
                             true))
      printf("Planes %ld and %ld were both located within Europe\n",
             planes[i].trj_id, planes[j].trj_id);
    else if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[4],
                             true) &&
             tpoint_at_stbox((const Temporal *)planes[j].trip, continentsSTB[4],
                             true))
      printf("Planes %ld and %ld were both located within North_America\n",
             planes[i].trj_id, planes[j].trj_id);
    else if (tpoint_at_stbox((const Temporal *)planes[i].trip, continentsSTB[5],
                             true) &&
             tpoint_at_stbox((const Temporal *)planes[j].trip, continentsSTB[5],
                             true))
      printf("Planes %ld and %ld were both located within Oceania\n",
             planes[i].trj_id, planes[j].trj_id);
  }
}

void query4(int num_planes, trip_record planes[MAX_TRIPS],
            area_record airports[MAX_AIRPORTS]) {
  printf("4- List one airport each plane were close to, if any\n");
  for (int i = 0; i < num_planes; i++) {
    for (int j = 0; j < MAX_AIRPORTS; j++) {
      if (edwithin_tpoint_geo((const Temporal *)planes[i].trip,
                              airports[j].geom, 100)) {
        printf("Plane %ld visited airport %d\n", planes[i].trj_id, j);
        break;
      }
    }
  }
}

void query5(int num_planes, trip_record planes[MAX_TRIPS]) {
  printf("5- Count the number of active trips per minute in the first hour on "
         "May 25, 2020.\n");

  // Array to hold the counts for each minute
  int counts[60] = {0};

  // Loop over each minute of the first hour (00:00 to 01:00)
  for (int minute = 0; minute < 60; minute++) {
    // Create start and end timestamp strings for the current minute
    char start_str[30];
    char end_str[30];
    snprintf(start_str, sizeof(start_str), "2020-05-25 00:%02d:00+00", minute);
    snprintf(end_str, sizeof(end_str), "2020-05-25 00:%02d:00+00", minute + 1);

    // Parse the timestamps
    TimestampTz start_time = pg_timestamptz_in(start_str, 0);
    TimestampTz end_time = pg_timestamptz_in(end_str, 0);

    // Create the time span for the current minute
    Span minute_span;
    span_set(start_time,    // Datum lower
             end_time,      // Datum upper
             true,          // lower_inc
             false,         // upper_inc
             T_TIMESTAMPTZ, // basetype
             T_TSTZSPAN,    // spantype
             &minute_span   // Span *s
    );

    int count = 0;

    // Iterate over each trip
    for (int i = 0; i < num_planes; i++) {
      // Assuming planes[i].trip->period is of type Span
      Span *trip_span = &(planes[i].trip->period);

      // Check if the trip's span overlaps with the current minute's span
      if (overlaps_span_span(&minute_span, trip_span)) {
        count++;
      }
    }

    counts[minute] = count;
    printf("Minute %02d:00 - %02d:00: %d active trips\n", minute, minute + 1,
           count);
  }
}

void query6(int num_planes, trip_record planes[MAX_TRIPS]) {
  printf("6- list the distance for each plane\n");
  for (int i = 0; i < num_planes; i++) {
    printf("Plane %ld distance %lf\n", planes[i].trj_id,
           tpoint_length((const Temporal *)planes[i].trip));
  }
}
void query7(int num_planes, trip_record planes[MAX_TRIPS]) {
  printf("7- List the minimum distance ever between each plane\n");
  for (int i = 0, j = 1; i < num_planes - 1; i++, j++) {
    double dist = nad_tpoint_tpoint((const Temporal *)planes[i].trip,
                                    (const Temporal *)planes[j].trip);
    printf("Planes %ld and %ld minimum distance %lf\n", planes[i].trj_id,
           planes[j].trj_id, dist);
  }
}
void query8(int num_planes, trip_record planes[MAX_TRIPS]) {
  printf("8- For each plane, list the Three Closer Planes\n");
  int **distances = malloc(num_planes * sizeof(int *));
  if (distances == NULL) {
    perror("Failed to allocate memory");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < num_planes; i++) {
    distances[i] = malloc(num_planes * sizeof(int));
    if (distances[i] == NULL) {
      perror("Failed to allocate memory");
      exit(EXIT_FAILURE);
    }
  }
  printf("Computing distances between planes\n");
  // Precompute all distances between planes
  for (int i = 0; i < num_planes; i++) {
    for (int j = i + 1; j < num_planes; j++) {
      if (i != j)
        distances[i][j] = nad_tpoint_tpoint((const Temporal *)planes[i].trip,
                                            (const Temporal *)planes[j].trip);
      else
        distances[i][j] = INT_MAX;
    }
  }
  printf("Distances computed\n");
  // Find and print the three closest planes for each plane
  for (int i = 0; i < num_planes; i++) {
    long int closest_planes[3] = {-1, -1, -1};
    int closest_distances[3] = {INT_MAX, INT_MAX, INT_MAX};
    for (int j = 0; j < num_planes; j++) {
      if (i != j) {
        int dist = distances[i][j];
        if (dist == 0)
          continue; // Skip distance zero

        // Find the appropriate position for the current distance
        for (int k = 0; k < 3; k++) {
          if (dist < closest_distances[k]) {
            // Shift the larger distances down
            for (int l = 2; l > k; l--) {
              closest_distances[l] = closest_distances[l - 1];
              closest_planes[l] = closest_planes[l - 1];
            }
            closest_distances[k] = dist;
            closest_planes[k] = planes[j].trj_id;
            break;
          }
        }
      }
    }
    printf("Plane %ld closest planes: %ld (distance %d), %ld (distance %d), "
           "%ld (distance %d)\n",
           planes[i].trj_id, closest_planes[0], closest_distances[0],
           closest_planes[1], closest_distances[1], closest_planes[2],
           closest_distances[2]);
  }
  // Free the allocated memory
  for (int i = 0; i < num_planes; i++) {
    free(distances[i]);
  }
  free(distances);
}

int main(int argc, char **argv) {
  plane_point_record rec;
  area_record continents[7] = {0};
  area_record airports[MAX_AIRPORTS] = {0};
  const STBox *continentsSTB[7] = {0};
  int no_records = 0;
  int no_nulls = 0;
  // int no_planeBB = 0;
  int no_writes = 0;
  const Interval *maxt = pg_interval_in("1 day", -1);
  char point_buffer[MAX_LENGTH_POINT];
  char text_buffer[MAX_LENGTH_HEADER];
  /* Allocate space to build the planes trip */
  trip_record planes[MAX_TRIPS] = {0};
  /* Number of planes */
  int num_planes = 0;
  /* Iterator variable */
  int i;
  /* Return value */
  int return_value = 0;
  // geoaltitude

  /* Initialize MEOS */
  meos_initialize(NULL, NULL);

  /***************************************************************************
   * Section 1: Create regions and airports
   ***************************************************************************/
  printf("Creating regions and airports\n");
  /* Africa */
  continentsSTB[0] = stbox_in("SRID=4326;STBOX XT(((-17.6, -34.8),(51.5, "
                              "37.5)),[2020-05-24,2020-05-26])");
  /* Antarctica */
  continentsSTB[1] = stbox_in(
      "SRID=4326;STBOX XT(((-180, -90), (180, -60)),[2020-05-24,2020-05-26])");
  /* Asia */
  continentsSTB[2] = stbox_in("SRID=4326;STBOX XT(((26.0, -10.0), (180.0, "
                              "81.0)),[2020-05-24,2020-05-26])");
  /* Europe */
  continentsSTB[3] = stbox_in("SRID=4326;STBOX XT(((-31.5, 34.5), (39.0, "
                              "81.0)),[2020-05-24,2020-05-26])");
  /* North America */
  continentsSTB[4] = stbox_in("SRID=4326;STBOX XT(((-168.0, 7.2), (-25.0, "
                              "83.2)),[2020-05-24,2020-05-26])");
  /* Oceania */
  continentsSTB[5] = stbox_in("SRID=4326;STBOX XT(((110.0, -50.0), (180.0, "
                              "0.0)),[2020-05-24,2020-05-26])");
  /* South America */
  continentsSTB[6] = stbox_in("SRID=4326;STBOX XT(((-81.0, -56.0), (-34.0, "
                              "13.0)),[2020-05-24,2020-05-26])");

  continents[0].geom = stbox_to_geo(continentsSTB[0]);
  continents[1].geom = stbox_to_geo(continentsSTB[1]);
  continents[2].geom = stbox_to_geo(continentsSTB[2]);
  continents[3].geom = stbox_to_geo(continentsSTB[3]);
  continents[4].geom = stbox_to_geo(continentsSTB[4]);
  continents[5].geom = stbox_to_geo(continentsSTB[5]);
  continents[6].geom = stbox_to_geo(continentsSTB[6]);

  /* Read the airports file */
  FILE *airportsUsaFile = fopen("data/airports_only_coord_small.csv", "r");
  if (!airportsUsaFile) {
    printf("Error opening airports file\n");
    return_value = 1;
    goto cleanup;
  }
  /* Read the first line of the file with the headers */
  readCreateAirports(airportsUsaFile, airports);

  /***************************************************************************
   * Section 2: Initialize MEOS, open the input CSV file
   ***************************************************************************/

  /* Get start time */
  clock_t t = clock();
  clock_t t1 = clock();

  /* You may substitute the full file path in the first argument of fopen */
  FILE *fileIn = fopen("data/filtered_states.csv", "r");

  printf("file in opened\n");

  if (!fileIn) {
    printf("Error opening input file\n");
    return_value = 1;
    goto cleanup;
  }
  /***************************************************************************
   * Section 3: Read input file line by line and append each observation as a
   * temporal point in MEOS
   ***************************************************************************/
  printf("Reading and creating trajectories\n");
  /* Read the first line of the file with the headers */
  fscanf(fileIn, "%73s\n", text_buffer);

  /* Continue reading the file */
  do {
    int read =
        fscanf(fileIn, "%19[^,],%ld,%f,%f,%f,%f,%f,%f,%f\n", text_buffer,
               &rec.trj_id, &rec.lat, &rec.lon, &rec.velocity, &rec.heading,
               &rec.vertrate, &rec.baroaltitude, &rec.geoaltitude);
    /* Transform the string representing the timestamp into a timestamp value */
    // printf("traject %ld, lat %f, lon %f, vel %f, heading %f, vertrate %f,
    // baroaltitude %f, geoaltitude %f\n", rec.trj_id, rec.lat, rec.lon,
    // rec.velocity, rec.heading, rec.vertrate, rec.baroaltitude,
    // rec.geoaltitude);
    rec.T = pg_timestamp_in(text_buffer, -1);

    if (read == 9)
      no_records++;
    if (read != 9 && !feof(fileIn)) {
      printf("Record with missing values ignored\n");
      no_nulls++;
    }
    if (ferror(fileIn)) {
      printf("Error reading file\n");
      fclose(fileIn);
      fclose(airportsUsaFile);
    }
    /* Find the place to store the new instant */
    int plane = -1;
    for (i = 0; i < MAX_TRIPS; i++) {
      if (planes[i].trj_id == rec.trj_id) {
        plane = i;
        break;
      }
    }
    if (plane < 0) {
      plane = num_planes++;
      if (plane == MAX_TRIPS) {
        printf(
            "The maximum number of planes in the input file is bigger than %d",
            MAX_TRIPS);
        return_value = 1;
        goto cleanup;
      }

      planes[plane].trj_id = rec.trj_id;
    }
    /*
     * Append the latest observation to the corresponding plane.
     * In the input file it is assumed that
     * - The coordinates are given in the WGS84 geographic coordinate system
     * - The timestamps are given in GMT time zone
     */
    char *t_out = pg_timestamp_out(rec.T);
    sprintf(point_buffer, "SRID=4326;Point(%lf %lf %f)@%s+00", rec.lon, rec.lat,
            rec.geoaltitude, t_out);

    /* Append the last observation */
    TInstant *inst = (TInstant *)tgeompoint_in(point_buffer);
    if (!planes[plane].trip)
      planes[plane].trip =
          tsequence_make_exp((const TInstant **)&inst, 1, NO_INSTANTS_BATCH,
                             true, true, LINEAR, false);
    else
      tsequence_append_tinstant(planes[plane].trip, inst, 1000, maxt, true);
  } while (!feof(fileIn));

  t1 = clock() - t1;
  double time_taken1 = ((double)t1) / CLOCKS_PER_SEC;
  printf("\nThe program took %f seconds to read and create trajectories\n",
         time_taken1);
  printf("%d records read.\n%d incomplete records ignored.\n%d planes\n",
         no_records, no_nulls, num_planes);

  /***************************************************************************
   * Section 4: Range Queries
   ***************************************************************************/

  /* 1- List the continents each plane pass */
  // query1(num_planes, planes, continents);

  /* 2- List planes that were within a region from Regions during a period from
   * Periods */
  // query2(num_planes, planes, continentsSTB);

  /* 3- List the pairs of planes that were both located within a region from
   * Regions during a period from Periods */
  // query3(num_planes, planes, continentsSTB);

  /* 4- List one airport each plane were close to, if any */
  // query4(num_planes, planes, airports);

  /***************************************************************************
   * Section 5: Temporal Aggregate Queries
   ***************************************************************************/

  /* 5- Count the number of trips that were active during each hour in May 25,
   * 2020 */
  query5(num_planes, planes);

  /***************************************************************************
   * Section 6: Distance Queries
   ***************************************************************************/

  /* 6- list the distance for each plane */
  // query6(num_planes, planes);

  /* 7- List the minimum distance ever between each plane */
  // query7(num_planes, planes);

  /***************************************************************************
   * Section 7: Nearest-Neighbor Query
   ***************************************************************************/

  /* 8- For each plane, list the three planes that are closest to it */
  // query8(num_planes, planes);

/* Clean up */
cleanup:
  /* Free memory */
  for (i = 0; i < num_planes; i++)
    free(planes[i].trip);
  /* Finalize MEOS */
  meos_finalize();
  /* Close the connection to the logfile */
  fclose(fileIn);
  fclose(airportsUsaFile);
  return return_value;
}