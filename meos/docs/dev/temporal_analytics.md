# Kalman cleaning (internal)

This section shows how to call the internal temporal Kalman cleaner from a C pipeline using MEOS. The cleaner applies a constant‑velocity (CV) linear Kalman Filter per axis and removes or replaces outliers based on an innovation gate.

Example (double2, sequence):

```c
#include <meos.h>
#include <meos_internal.h>
#include "temporal/temporal.h"
#include "temporal/tsequence.h"
#include "temporal/tkalman.h"

void run_kalman_example(void) {
  meos_initialize();
  meos_initialize_timezone("UTC");

  // Build a short 2D sequence
  TInstant *insts[3];
  double2 a; double2_set(0.0, 0.0, &a);
  double2 b; double2_set(1.0, 1.0, &b);
  double2 c; double2_set(50.0, -50.0, &c); // big outlier
  insts[0] = tinstant_make(Double2PGetDatum(&a), T_TDOUBLE2, 0);
  insts[1] = tinstant_make(Double2PGetDatum(&b), T_TDOUBLE2, 1000000);
  insts[2] = tinstant_make(Double2PGetDatum(&c), T_TDOUBLE2, 2000000);
  TSequence *seq = tsequence_make((const TInstant **)insts, 3, true, true, STEP, NORMALIZE_NO);

  // Defaults, drop outliers (set fill_estimates=true to replace instead)
  TKalmanParams P = tkalman_default_params();
  P.fill_estimates = false;
  int removed = 0;
  Temporal *clean = temporal_kalman_clean((Temporal *)seq, &P, &removed);

  // Use `clean` in downstream processing...
  // (e.g., serialize, analyze, etc.)

  pfree(clean);
  pfree(seq);
  meos_finalize();
}
```

Notes:
- Subtypes: supports TINSTANT, TSEQUENCE, TSEQUENCESET (resets per sequence).
- Interpolation: preserves DISCRETE, STEP, LINEAR.
- Types: float8 (1D), double2 (2D), double3 (3D), double4 (4D).
- Outliers: z‑score (D=1) or Mahalanobis (D>1), controlled by `gate_sigma`.
