#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <meos.h>
#include <meos_internal.h>
#include "temporal/temporal.h"
#include "temporal/tsequence.h"
#include "temporal/type_util.h"
#include "temporal/tkalman.h"

static double now_us(void)
{
  struct timeval tv; gettimeofday(&tv, NULL);
  return tv.tv_sec * 1e6 + tv.tv_usec;
}

static TSequence *make_seq_1d(int n)
{
  TInstant **insts = palloc(sizeof(TInstant *) * n);
  for (int i = 0; i < n; i++)
  {
    double v = sin(0.01 * i);
    insts[i] = tinstant_make(Float8GetDatum(v), T_TFLOAT, (TimestampTz) (i * 1000));
  }
  TSequence *seq = tsequence_make((const TInstant **) insts, n, true, true, STEP, NORMALIZE_NO);
  pfree_array((void **) insts, n);
  return seq;
}

static TSequence *make_seq_2d(int n)
{
  TInstant **insts = palloc(sizeof(TInstant *) * n);
  for (int i = 0; i < n; i++)
  {
    double2 d; double2_set(sin(0.01 * i), cos(0.01 * i), &d);
    insts[i] = tinstant_make(Double2PGetDatum(&d), T_TDOUBLE2, (TimestampTz) (i * 1000));
  }
  TSequence *seq = tsequence_make((const TInstant **) insts, n, true, true, STEP, NORMALIZE_NO);
  pfree_array((void **) insts, n);
  return seq;
}

static TSequence *make_seq_3d(int n)
{
  TInstant **insts = palloc(sizeof(TInstant *) * n);
  for (int i = 0; i < n; i++)
  {
    double3 d; double3_set(sin(0.01 * i), cos(0.01 * i), 0.1 * i, &d);
    insts[i] = tinstant_make(Double3PGetDatum(&d), T_TDOUBLE3, (TimestampTz) (i * 1000));
  }
  TSequence *seq = tsequence_make((const TInstant **) insts, n, true, true, STEP, NORMALIZE_NO);
  pfree_array((void **) insts, n);
  return seq;
}

static void bench_seq(TSequence *seq, const char *label)
{
  TKalmanParams P = tkalman_default_params();
  P.fill_estimates = false;
  const int n = seq->count;
  double t0 = now_us();
  int removed = 0;
  Temporal *clean = temporal_kalman_clean((Temporal *) seq, &P, &removed);
  double t1 = now_us();
  double usec = (t1 - t0);
  double per = usec / (double) n;
  printf("%s: n=%d time=%.0f us (%.3f us/inst) removed=%d\n", label, n, usec, per, removed);
  pfree(clean);
}

int main(void)
{
  meos_initialize();
  meos_initialize_timezone("UTC");

  const int N = 100000; // 1e5
  TSequence *s1 = make_seq_1d(N);
  TSequence *s2 = make_seq_2d(N);
  TSequence *s3 = make_seq_3d(N);

  bench_seq(s1, "D=1");
  bench_seq(s2, "D=2");
  bench_seq(s3, "D=3");

  pfree(s1); pfree(s2); pfree(s3);
  meos_finalize();
  return 0;
}
