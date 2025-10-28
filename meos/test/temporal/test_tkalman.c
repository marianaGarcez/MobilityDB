#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <meos.h>
#include <meos_internal.h>
#include "temporal/temporal.h"
#include "temporal/tsequence.h"
#include "temporal/tsequenceset.h"
#include "temporal/type_util.h"
#include "temporal/tkalman.h"

static inline double gauss_rand(unsigned *seed)
{
  // Box-Muller transform (deterministic via custom LCG)
  // Simple LCG instead of rand() for reproducibility across platforms
  *seed = (*seed) * 1664525u + 1013904223u;
  double u1 = ((*seed) & 0xFFFFFFu) / (double)0x1000000u;
  *seed = (*seed) * 1664525u + 1013904223u;
  double u2 = ((*seed) & 0xFFFFFFu) / (double)0x1000000u;
  if (u1 <= 1e-12) u1 = 1e-12;
  return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

static TSequence *build_1d_sine_with_spikes(int n, double noise_sigma,
                                             const int *spikes, int nspikes,
                                             TimestampTz t0, int64_t step_us)
{
  TInstant **insts = palloc(sizeof(TInstant *) * n);
  unsigned seed = 12345u;
  for (int i = 0; i < n; i++)
  {
    double base = sin(0.1 * i);
    double noise = noise_sigma * gauss_rand(&seed);
    double val = base + noise;
    for (int k = 0; k < nspikes; k++)
      if (i == spikes[k]) val += 50.0; // large spike
    TimestampTz t = t0 + (TimestampTz) step_us * i;
    insts[i] = tinstant_make(Float8GetDatum(val), T_TFLOAT, t);
  }
  TSequence *seq = tsequence_make((const TInstant **) insts, n,
                                  true, true, STEP, NORMALIZE_NO);
  pfree_array((void **) insts, n);
  return seq;
}

static TSequence *build_2d_randomwalk_with_bursts(int n,
                                                  const int *bursts, int nbursts,
                                                  TimestampTz t0, int64_t step_us)
{
  TInstant **insts = palloc(sizeof(TInstant *) * n);
  unsigned seed = 54321u;
  double x = 0.0, y = 0.0;
  for (int i = 0; i < n; i++)
  {
    x += 0.5 * gauss_rand(&seed);
    y += 0.5 * gauss_rand(&seed);
    for (int k = 0; k < nbursts; k++)
      if (i == bursts[k]) { x += 30.0; y -= 30.0; }
    double2 *d = palloc(sizeof(double2));
    double2_set(x, y, d);
    TimestampTz t = t0 + (TimestampTz) step_us * i;
    insts[i] = tinstant_make_free(Double2PGetDatum(d), T_TDOUBLE2, t);
  }
  TSequence *seq = tsequence_make((const TInstant **) insts, n,
                                  true, true, STEP, NORMALIZE_NO);
  pfree_array((void **) insts, n);
  return seq;
}

static void test_1d_spikes(void)
{
  int spikes[10] = {10,30,50,70,90,110,130,150,170,190};
  TSequence *seq = build_1d_sine_with_spikes(200, 0.05, spikes, 10,
                                             0, 1000000);

  // Defaults
  TKalmanParams P = tkalman_default_params();
  P.fill_estimates = false;
  int removed = 0;
  Temporal *clean = temporal_kalman_clean((Temporal *) seq, &P, &removed);
  assert(clean != NULL);
  const TSequence *cseq = (const TSequence *) clean;
  // Expect approximately 10 removals (deterministic spikes at +50)
  assert(removed >= 9 && removed <= 11);
  assert(cseq->count == seq->count - removed);
  // Ensure spike timestamps were removed
  for (int s = 0; s < 10; s++)
  {
    TimestampTz ts = seq->period.lower + (TimestampTz) (1000000LL * spikes[s]);
    bool found = false;
    for (int i = 0; i < cseq->count; i++)
      if (TSEQUENCE_INST_N(cseq, i)->t == ts) { found = true; break; }
    assert(!found);
  }
  pfree(clean);

  // fill_estimates = true
  P.fill_estimates = true;
  removed = 0;
  clean = temporal_kalman_clean((Temporal *) seq, &P, &removed);
  assert(clean != NULL);
  cseq = (const TSequence *) clean;
  assert(removed == 0);
  assert(cseq->count == seq->count);
  // Values at spike timestamps should be bounded (no large spikes)
  for (int s = 0; s < 10; s++)
  {
    TimestampTz ts = seq->period.lower + (TimestampTz) (1000000LL * spikes[s]);
    bool found = false;
    double val_at_spike = 0.0;
    for (int i = 0; i < cseq->count; i++)
      if (TSEQUENCE_INST_N(cseq, i)->t == ts) {
        found = true;
        val_at_spike = DatumGetFloat8(tinstant_value_p(TSEQUENCE_INST_N(cseq, i)));
        break;
      }
    assert(found);
    assert(fabs(val_at_spike) < 10.0);
  }
  pfree(clean);

  pfree(seq);
}

static void test_2d_randomwalk_bursts(void)
{
  int bursts[8] = {12,25,40,60,80,100,140,180};
  TSequence *seq = build_2d_randomwalk_with_bursts(200, bursts, 8,
                                                   0, 1000000);
  TKalmanParams P = tkalman_default_params();
  P.fill_estimates = false;
  int removed = 0;
  Temporal *clean = temporal_kalman_clean((Temporal *) seq, &P, &removed);
  assert(clean != NULL);
  const TSequence *cseq = (const TSequence *) clean;
  // With large vector bursts, expect most to be rejected
  assert(removed >= 6 && removed <= 8);
  assert(cseq->count == seq->count - removed);
  pfree(clean);

  // fill_estimates keeps length
  P.fill_estimates = true;
  removed = 0;
  clean = temporal_kalman_clean((Temporal *) seq, &P, &removed);
  assert(clean != NULL);
  cseq = (const TSequence *) clean;
  assert(removed == 0);
  assert(cseq->count == seq->count);
  pfree(clean);

  pfree(seq);
}

static void test_sequenceset_boundaries_and_gaps(void)
{
  // Two sequences separated by a large gap; spike at first instant of seq2
  const int n1 = 20, n2 = 20;
  TInstant **insts1 = palloc(sizeof(TInstant *) * n1);
  TInstant **insts2 = palloc(sizeof(TInstant *) * n2);
  for (int i = 0; i < n1; i++)
  {
    double val = sin(0.2 * i);
    insts1[i] = tinstant_make(Float8GetDatum(val), T_TFLOAT, (TimestampTz) (i * 1000000LL));
  }
  for (int i = 0; i < n2; i++)
  {
    double val = (i == 0) ? 100.0 : sin(0.2 * i);
    insts2[i] = tinstant_make(Float8GetDatum(val), T_TFLOAT, (TimestampTz) (3600LL * 1000000LL + i * 1000000LL));
  }
  TSequence *seq1 = tsequence_make((const TInstant **) insts1, n1,
                                         true, true, STEP, NORMALIZE_NO);
  TSequence *seq2 = tsequence_make((const TInstant **) insts2, n2,
                                         true, true, STEP, NORMALIZE_NO);
  const TSequence *seqs_arr[2] = { (const TSequence *) seq1, (const TSequence *) seq2 };
  TSequenceSet *ss = tsequenceset_make(seqs_arr, 2, NORMALIZE_NO);
  pfree_array((void **) insts1, n1);
  pfree_array((void **) insts2, n2);

  TKalmanParams P = tkalman_default_params();
  P.fill_estimates = false;
  int removed = 0;
  Temporal *clean = temporal_kalman_clean((Temporal *) ss, &P, &removed);
  assert(clean != NULL);
  const TSequenceSet *css = (const TSequenceSet *) clean;
  assert(css->count == 2);
  const TSequence *cseq2 = TSEQUENCESET_SEQ_N(css, 1);
  // Ensure first instant of second sequence remains (filter resets)
  assert(TSEQUENCE_INST_N(cseq2, 0)->t == TSEQUENCE_INST_N(seq2, 0)->t);
  pfree(clean);

  pfree((void *) seq1);
  pfree((void *) seq2);
  pfree(ss);
}

static void test_edge_cases(void)
{
  // Single instant
  TInstant *inst = tinstant_make(Float8GetDatum(1.0), T_TFLOAT, 0);
  TKalmanParams P = tkalman_default_params();
  int removed = 0;
  Temporal *t = temporal_kalman_clean((Temporal *) inst, &P, &removed);
  assert(t != NULL);
  assert(((TInstant *) t)->t == inst->t);
  assert(removed == 0);
  pfree(t);
  pfree(inst);

  // All points are outliers (except first used for init)
  const int n = 10;
  TInstant **insts = palloc(sizeof(TInstant *) * n);
  for (int i = 0; i < n; i++)
  {
    double val = (i == 0) ? 0.0 : 100.0;
    insts[i] = tinstant_make(Float8GetDatum(val), T_TFLOAT, (TimestampTz) (i * 1000000LL));
  }
  TSequence *seq = tsequence_make((const TInstant **) insts, n,
                                  true, true, STEP, NORMALIZE_NO);
  pfree_array((void **) insts, n);
  P.fill_estimates = false;
  removed = 0;
  Temporal *clean = temporal_kalman_clean((Temporal *) seq, &P, &removed);
  assert(clean != NULL);
  assert(((TSequence *) clean)->count == 1);
  assert(removed == n - 1);
  pfree(clean);

  // All points identical (constant series)
  insts = palloc(sizeof(TInstant *) * n);
  for (int i = 0; i < n; i++)
    insts[i] = tinstant_make(Float8GetDatum(5.0), T_TFLOAT, (TimestampTz) (i * 1000000LL));
  seq = tsequence_make((const TInstant **) insts, n,
                       true, true, STEP, NORMALIZE_NO);
  pfree_array((void **) insts, n);
  removed = 0;
  clean = temporal_kalman_clean((Temporal *) seq, &P, &removed);
  assert(clean != NULL);
  assert(((TSequence *) clean)->count == n);
  assert(removed == 0);
  pfree(clean);

  // Very irregular sampling
  const int m = 12;
  insts = palloc(sizeof(TInstant *) * m);
  int64_t steps[12] = {1,10,2,20,3,30,4,40,5,50,6,60};
  TimestampTz tcur = 0;
  for (int i = 0; i < m; i++)
  {
    tcur += steps[i] * 1000000LL;
    insts[i] = tinstant_make(Float8GetDatum(sin(0.2 * i)), T_TFLOAT, tcur);
  }
  seq = tsequence_make((const TInstant **) insts, m,
                       true, true, STEP, NORMALIZE_NO);
  pfree_array((void **) insts, m);
  removed = 0;
  clean = temporal_kalman_clean((Temporal *) seq, &P, &removed);
  assert(clean != NULL);
  assert(((TSequence *) clean)->count == m);
  pfree(clean);

  pfree(seq);
}

int main(void)
{
  meos_initialize();
  meos_initialize_timezone("UTC");

  test_1d_spikes();
  test_2d_randomwalk_bursts();
  test_sequenceset_boundaries_and_gaps();
  test_edge_cases();

  meos_finalize();
  printf("test_tkalman: OK\n");
  return 0;
}
