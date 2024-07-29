/*****************************************************************************
 *
 * This MobilityDB code is provided under The PostgreSQL License.
 * Copyright (c) 2016-2024, Université libre de Bruxelles and MobilityDB
 * contributors
 *
 * MobilityDB includes portions of PostGIS version 3 source code released
 * under the GNU General Public License (GPLv2 or later).
 * Copyright (c) 2001-2024, PostGIS contributors
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
 * @brief Port Arrival and Departure Times by MMSI.
 * This query showcases event detection and Temporal Analysis. 
 * It detects the times vessels arrive at and depart from specified ports. 
 * Based on its AIS points and MEOS' atGeom function, it identifies when a vessel 
 * enters a port area and when it leaves.
 *
 * The program can be build as follows
 * @code
 * linux:
 * gcc -Wall -g -I/usr/local/include -o query3 query3.c -L/usr/local/lib -lmeos
 * for macOS:
 * clang -std=gnu99 query3.c -I/opt/homebrew/include  -L/opt/homebrew/lib/ -lmeos -o 03
 * @endcode
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <meos.h>
/* The expandable functions are in the internal MEOS API */
#include <meos_internal.h>

/* Number of instants to send in batch to the file  */
#define NO_INSTANTS_BATCH 100
/* Number of instants to keep when restarting a sequence */
#define NO_INSTANTS_KEEP 2
/* Maximum length in characters of a header record in the input CSV file */
#define MAX_LENGTH_HEADER 1024
/* Maximum length in characters of a point in the input data */
#define MAX_LENGTH_POINT 64
/* Maximum number of trips */
#define MAX_TRIPS 5

typedef struct
{
  Timestamp T;
  long int MMSI;
  double Latitude;
  double Longitude;
  double SOG;
} AIS_record;

typedef struct
{
  long int MMSI;   /* Identifier of the trip */
  TSequence *trip; /* Array of latest observations of the trip, by MMSI */
  bool in_port;    /* True if the ship is in the port */
  int num;         /* Number of times the ship entered the port */
} trip_record;


int time_in_port(const Temporal * traj, const STBox *port, bool in_port){
  if (traj == NULL) {
      printf("traj is NULL\n");
      return 0;
  }
  if (tpoint_at_stbox(traj,port, true)){
    if (in_port == false){
      in_port = true;
      //printf("Ship just entered the port\n");
      return 1;
    }
  }
  else {
    if (in_port == true){
      in_port = false;
      //printf("Ship just left the port\n");
      return 0;
    }
  }
  return 0;
}

int
main(int argc, char **argv)
{
  AIS_record rec;
  int no_records = 0;
  int no_nulls = 0;
  char point_buffer[MAX_LENGTH_POINT];
  char text_buffer[MAX_LENGTH_HEADER];
  /* Allocate space to build the trips */
  trip_record trips[MAX_TRIPS] = {0};
  /* Number of ships */
  int no_ships = 0;
  /* Number of writes */
  int no_writes = 0;
  /* Iterator variables */
  int i, j;
  /* Exit value initialized to 1 (i.e., error) to quickly exit upon error */
  int exit_value = 1;

  /* Initialize MEOS */
  meos_initialize(NULL, NULL);

  /***************************************************************************
   * Section 0: Define the bounding box and time period
   ***************************************************************************/
  const STBox *bay = stbox_in("SRID=4326;STBOX X((3.3615, 53.964367),(3.505853, 55.24544))");
  printf("Created port A\n");

  const STBox *airport = stbox_in("SRID=3414;STBOX X((43239.013,31804.933),(43275.64002746976,39423.317))");
  printf("Created port B\n");

  /***************************************************************************
   * Section 1: Open the output file
   ***************************************************************************/

  /* Open/create the output file */
  /* You may substitute the full file path in the first argument of fopen */
  FILE *file_out = fopen("../data/out/ais_instantsOutQuery3.csv", "w+");

  if (! file_out)
  {
    printf("Error creating/opening the output file\n");
    goto cleanup;
  }

  /***************************************************************************
   * Section 2: Open the input AIS file
   ***************************************************************************/

  /* You may substitute the full file path in the first argument of fopen */
  FILE *file_in = fopen("../data/ais_instants.csv", "r");

  if (! file_in)
  {
    printf("Error opening input file\n");
    goto cleanup;
  }

  /***************************************************************************
   * Section 3: Read input file line by line and append each observation as a
   * temporal point in MEOS
   ***************************************************************************/

  printf("Accumulating %d instants before sending them to the output file\n"
    "(one '*' marker every output file update)\n", NO_INSTANTS_BATCH);

  /* Read the first line of the file with the headers */
  fscanf(file_in, "%1023s\n", text_buffer);

  /* Continue reading the file */
  do
  {
    int read = fscanf(file_in, "%32[^,],%ld,%lf,%lf,%lf\n",
      text_buffer, &rec.MMSI, &rec.Latitude, &rec.Longitude, &rec.SOG);
    /* Transform the string representing the timestamp into a timestamp value */
    rec.T = pg_timestamp_in(text_buffer, -1);

    if (read == 5)
      no_records++;

    if (read != 5 && ! feof(file_in))
    {
      printf("Record with missing values ignored\n");
      no_nulls++;
    }

    if (ferror(file_in))
    {
      printf("Error reading input file\n");
      fclose(file_in);
      fclose(file_out);
    }

    /* Find the place to store the new instant */
    j = -1;
    for (i = 0; i < no_ships; i++)
    {
      if (trips[i].MMSI == rec.MMSI)
      {
        j = i;
        break;
      }
    }
    if (j < 0)
    {
      j = no_ships++;
      if (j == MAX_TRIPS)
      {
        printf("The maximum number of ships in the input file is bigger than %d",MAX_TRIPS);
        goto cleanup;
      }
      trips[j].MMSI = rec.MMSI;
    }
    /*
     * Append the latest observation to the corresponding ship.
     * In the input file it is assumed that
     * - The coordinates are given in the WGS84 geographic coordinate system
     * - The timestamps are given in GMT time zone
     */
    char *t_out = pg_timestamp_out(rec.T);
    sprintf(point_buffer, "SRID=4326;Point(%lf %lf)@%s+00", rec.Longitude,
      rec.Latitude, t_out);
    free(t_out);

    /* Send to the output file the trip if reached the maximum number of instants */
    if (trips[j].trip && trips[j].trip->count == NO_INSTANTS_BATCH)
    {
      char *temp_out = tsequence_out(trips[j].trip, 15);
      fprintf(file_out, "%ld, %s, %d\n",trips[j].MMSI, temp_out, trips[j].num);
      /* Free memory */
      free(temp_out);
      no_writes++;
      printf("*");
      fflush(stdout);
      /* Restart the sequence by only keeping the last instants */
      tsequence_restart(trips[j].trip, NO_INSTANTS_KEEP);
      trips[j].in_port = false;
      trips[j].num = 0;
    }
  
    /* Append the last observation */
    TInstant *inst = (TInstant *) tgeogpoint_in(point_buffer);
    if (! trips[j].trip){
      trips[j].trip = tsequence_make_exp((const TInstant **) &inst, 1, NO_INSTANTS_BATCH, true, true, LINEAR, false);
      trips[j].in_port = false;
      trips[j].num = 0;
    }
    else{
      tsequence_append_tinstant(trips[j].trip, inst, 0.0, NULL, true);
    }
    TInstant *inst2 = (TInstant *) tgeompoint_in(point_buffer);
    int num = time_in_port((const Temporal *)inst2, bay, trips[j].in_port);
    trips[j].num += num;
    free(inst);
  } while (! feof(file_in));

  printf("\n%d records read\n%d incomplete records ignored\n"
    "%d writes to the output file\n", no_records, no_nulls, no_writes);

  /* Close the file */
  fclose(file_in);

  /* State that the program executed successfully */
  exit_value = 0;
  
/* Clean up */
cleanup:
 /* Free memory */
  for (i = 0; i < no_ships; i++)
    free(trips[i].trip);

  /* Finalize MEOS */
  meos_finalize();

  /* Close the connection to the output file */
  fclose(file_out);

  return exit_value;
}
