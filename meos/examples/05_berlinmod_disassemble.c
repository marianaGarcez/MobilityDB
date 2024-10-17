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
 * @brief A simple program that reads from a CSV file synthetic trip data in
 * Brussels generated by the MobilityDB-BerlinMOD generator,
 * https://github.com/MobilityDB/MobilityDB-BerlinMOD
 * disassembles the trips into individual observations, and write them in a
 * CSV file named "berlinmod_instants.csv" ordered by timestamp.
 *
 * The input file is
 * - `berlinmod_trips.csv`: 55 trips from 5 cars during 4 days obtained from the
 *   generator at scale factor 0.005. The input file has been generated with
 *   the following SQL command on the database containing the generated data
 * @code
 * COPY (SELECT tripid, vehid, day, seqno, ashexewkb(trip) AS trip FROM trips WHERE vehid < 6 ORDER BY tripid) TO '/home/user/src/berlinmod_trips.csv' CSV HEADER;
 * @endcode
 * In the above file, the coordinates are given in the 3857 coordinate system,
 * https://epsg.io/3857
 * and the timestamps are given in the Europe/Brussels time zone.
 * This simple program does not cope with erroneous inputs, such as missing
 * fields or invalid values.
 *
 * The program can be build as follows
 * @code
 * gcc -Wall -g -I/usr/local/include -o 05_berlinmod_disassemble 05_berlinmod_disassemble.c -L/usr/local/lib -lmeos
 * @endcode
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <meos.h>

/* Maximum length in characters of a trip in the input data */
#define MAX_LENGTH_TRIP 170001
/* Maximum length in characters of a header in the input CSV file */
#define MAX_LENGTH_HEADER 1024
/* Maximum length in characters of a date in the input data */
#define MAX_LENGTH_DATE 12
/* Maximum number of trips */
#define MAX_NO_TRIPS 64

typedef struct
{
  int tripid;
  int vehid;
  DateADT day;
  int seq;
  Temporal *trip;
} trip_record;


/* Main program */
int main(void)
{
  /* Variables to read the input CSV file */
  char header_buffer[MAX_LENGTH_HEADER];
  char date_buffer[MAX_LENGTH_DATE];
  char trip_buffer[MAX_LENGTH_TRIP];
  /* Arrays to compute the results */
  trip_record trips[MAX_NO_TRIPS] = {0};
  int curr_inst[MAX_NO_TRIPS];

  /* Get start time */
  clock_t t;
  t = clock();

  /* Initialize MEOS */
  meos_initialize();

  /* You may substitute the full file path in the first argument of fopen */
  FILE *file = fopen("data/berlinmod_trips.csv", "r");

  if (! file)
  {
    printf("Error opening input file\n");
    return 1;
  }

  int i = 0;

  /* Read the first line of the file with the headers */
  fscanf(file, "%1023s\n", header_buffer);

  /* Continue reading the file */
  do
  {
    int tripid, vehid, seq;
    int read = fscanf(file, "%d,%d,%10[^,],%d,%170000[^\n]\n",
      &tripid, &vehid, date_buffer, &seq, trip_buffer);
    /* Transform the string representing the date into a date value */
    DateADT day = date_in(date_buffer);
    /* Transform the string representing the trip into a temporal value */
    Temporal *trip = temporal_from_hexwkb(trip_buffer);

    /* Save the trip record */
    trips[i].vehid = vehid;
    trips[i].tripid = tripid;
    trips[i].seq = seq;
    trips[i].day = day;
    trips[i].trip = trip;

    if (read == 5)
      i++;

    if (read != 5 && !feof(file))
    {
      printf("Trip record with missing values\n");
      fclose(file);
      /* Free memory */
      for (int j = 0; j < i; j++)
        free(trips[j].trip);
      return 1;
    }

    if (ferror(file))
    {
      printf("Error reading input file\n");
      fclose(file);
      /* Free memory */
      for (int j = 0; j < i; j++)
        free(trips[j].trip);
      return 1;
    }
  } while (!feof(file));

  int records_in = i;

  /* Close the input file */
  fclose(file);

  /* Open the output file */
  file = fopen("data/berlinmod_instants.csv", "w+");

  /* Write the header line */
  fprintf(file,"tripid,vehid,day,seqno,geom,t\n");

  /* Initialize the current instant for each trip to the first one */
  for (i = 0; i < MAX_NO_TRIPS; i++)
    curr_inst[i] = 1;

  /* Loop until all trips have been processed */
  int records_out = 0;
  while (true)
  {
    /* Take as minimum instant the first instant of the first remaining trip */
    int first = 0;
    while (first < records_in && curr_inst[first] < 0)
      first++;
    if (first == records_in)
      /* All trips have been processed */
      break;
    TInstant *min_inst = temporal_instant_n(trips[first].trip,
      curr_inst[first]);
    int min_trip = first;

    /* Loop for the minimum instant among all remaining trips */
    for (i = first + 1; i < records_in; i++)
    {
      if (curr_inst[i] < 0)
        continue;
      TInstant *inst = temporal_instant_n(trips[i].trip, curr_inst[i]);
      if (min_inst->t > inst->t)
      {
        free(min_inst);
        min_inst = inst;
        min_trip = i;
      }
      else
        free(inst);
    }

    /* Write line in the CSV file */
    char *date_str = date_out(trips[min_trip].day);
    char *geom_str = geo_as_ewkt((GSERIALIZED *) &min_inst->value, 6);
    char *time_str = timestamptz_out(min_inst->t);
    fprintf(file,"%d,%d,%s,%d,%s,%s\n", trips[min_trip].vehid,
      trips[min_trip].vehid, date_str, trips[min_trip].seq, geom_str, time_str);
    free(date_str); free(geom_str); free(time_str); free(min_inst);
    records_out++;

    /* Advance the current instant of the trip */
    curr_inst[min_trip]++;
    if (curr_inst[min_trip] > temporal_num_instants(trips[min_trip].trip))
      curr_inst[min_trip] = -1;
  }

  printf("%d trip records read from file 'berlimod_trips.csv'.", records_in);
  printf("\n%d observation records written in file 'berlimod_instants.csv'.\n", records_out);

  /* Calculate the elapsed time */
  t = clock() - t;
  double time_taken = ((double) t) / CLOCKS_PER_SEC;
  printf("The program took %f seconds to execute\n", time_taken);

  /* Free memory */
  for (i = 0; i < records_in; i++)
    free(trips[i].trip);

  /* Close the ouput file */
  fclose(file);

  /* Finalize MEOS */
  meos_finalize();

  return 0;
}
