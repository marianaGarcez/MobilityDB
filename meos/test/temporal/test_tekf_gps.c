#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <meos.h>
#include <meos_internal.h>
#include "temporal/temporal.h"
#include "temporal/tsequence.h"
#include "temporal/tekf.h"
#include "temporal/type_util.h"

static inline double gauss(unsigned *seed){ *seed = (*seed)*1664525u + 1013904223u; double u1 = ((*seed)&0xFFFFFFu)/(double)0x1000000u; *seed = (*seed)*1664525u + 1013904223u; double u2 = ((*seed)&0xFFFFFFu)/(double)0x1000000u; if(u1<=1e-12) u1=1e-12; return sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);} 

static bool z_from_pos_value(Datum value, meosType tt, double *z, void *vctx){ TEkfGpsCtx *ctx=(TEkfGpsCtx*)vctx; int D=ctx->ndim, M=ctx->M; const double *a=ctx->anchors; if(temptype_basetype(tt)==T_DOUBLE3 && D==3){ const double3 *p=DatumGetDouble3P(value); for(int i=0;i<M;i++){ double dx=p->a-a[i*D+0], dy=p->b-a[i*D+1], dz=p->c-a[i*D+2]; z[i]=sqrt(dx*dx+dy*dy+dz*dz);} return true;} if(temptype_basetype(tt)==T_DOUBLE2 && D==2){ const double2 *p=DatumGetDouble2P(value); for(int i=0;i<M;i++){ double dx=p->a-a[i*D+0], dy=p->b-a[i*D+1]; z[i]=sqrt(dx*dx+dy*dy);} return true;} return false; }

static double mse_positions(const TSequence *seq, const double *truth){ double se=0; int n=seq->count; for(int i=0;i<n;i++){ const TInstant*inst=TSEQUENCE_INST_N(seq,i); meosType bt=temptype_basetype(inst->temptype); if(bt==T_DOUBLE3){ const double3*p=DatumGetDouble3P(tinstant_value_p(inst)); double dx=p->a-truth[3*i+0], dy=p->b-truth[3*i+1], dz=p->c-truth[3*i+2]; se+=dx*dx+dy*dy+dz*dz; } else if(bt==T_DOUBLE2){ const double2*p=DatumGetDouble2P(tinstant_value_p(inst)); double dx=p->a-truth[3*i+0], dy=p->b-truth[3*i+1]; se+=dx*dx+dy*dy; } } return se/(double)n; }

int main(void){
  meos_initialize(); meos_initialize_timezone("UTC");
  // Use CV model (positions as measurements). This reduces to a standard KF.
  TEkfModel model; assert(tekf_make_cv_model(3, &model));
  TEkfCvCtx cv = { .D = 3, .q_accel_var = 0.25, .r_meas_var = 4.0 };
  TEkfParams params={0}; params.default_dt=1.0; params.gate_sigma=3.5; params.fill_estimates=true;
  const int n=60; TInstant **insts=palloc(sizeof(TInstant*)*n); double *truth=palloc(sizeof(double)*3*n); unsigned seed=4242u;
  for(int i=0;i<n;i++){ double t=(double)i; double tx=5.0*t, ty=3.0*t, tz=10.0; truth[3*i+0]=tx; truth[3*i+1]=ty; truth[3*i+2]=tz; double nx=tx+2.0*gauss(&seed), ny=ty+2.0*gauss(&seed), nz=tz+2.0*gauss(&seed); if(i%15==0 && i>0){ nx+=60.0; ny-=60.0; } double3 d; double3_set(nx,ny,nz,&d); insts[i]=tinstant_make(Double3PGetDatum(&d), T_TDOUBLE3, (TimestampTz)(i*1000000LL)); }
  TSequence *seq=tsequence_make((const TInstant**)insts, n, true, true, STEP, NORMALIZE_NO); pfree_array((void**)insts,n);
  double mse_in=mse_positions(seq,truth);
  int removed=0; Temporal *clean=temporal_ekf_clean((Temporal*)seq,&model,&params,&cv,&removed);
  assert(clean);
  const TSequence*cseq=(const TSequence*)clean;
  assert(cseq->count==seq->count);
  double mse_out=mse_positions(cseq,truth);
  assert(mse_out < 0.5*mse_in);
  pfree(clean); pfree(seq); pfree(truth);
  meos_finalize(); printf("test_tekf_gps: OK\n"); return 0;
}
