/* See plan: EKF engine + CV + GPS models */
#include "temporal/tekf.h"

#if defined(MEOS_EXPERIMENTAL_ANALYTICS)

#include <math.h>
#include <string.h>
#include <meos.h>
#include <meos_internal.h>
#include "temporal/temporal.h"
#include "temporal/doublen.h"
#include "temporal/tsequence.h"
#include "temporal/tsequenceset.h"
#include "temporal/type_util.h"

static inline void mat_clear(double *A, int m, int n){ memset(A,0,sizeof(double)*(size_t)m*(size_t)n);} 
static inline void mat_copy(double *d,const double*s,int m,int n){ memcpy(d,s,sizeof(double)*(size_t)m*(size_t)n);} 
static inline void mat_add(double *C,const double*A,const double*B,int m,int n){int mn=m*n;for(int i=0;i<mn;i++)C[i]=A[i]+B[i];}
static inline void mat_sub(double *C,const double*A,const double*B,int m,int n){int mn=m*n;for(int i=0;i<mn;i++)C[i]=A[i]-B[i];}
static inline void mat_mul(const double*A,const double*B,double*C,int arows,int acols,int bcols){for(int i=0;i<arows;i++){for(int j=0;j<bcols;j++){double s=0;for(int k=0;k<acols;k++)s+=A[i*acols+k]*B[k*bcols+j];C[i*bcols+j]=s;}}}
static inline void mat_transpose(const double*A,double*At,int m,int n){for(int i=0;i<m;i++)for(int j=0;j<n;j++)At[j*m+i]=A[i*n+j];}

static bool chol_decompose(const double *A, double *L, int m){
  mat_copy(L,A,m,m);
  for(int i=0;i<m;i++){
    for(int j=i;j<m;j++){
      double sum=L[i*m+j];
      for(int k=0;k<i;k++) sum-=L[i*m+k]*L[j*m+k];
      if(i==j){ if(sum<=0.0) return false; L[i*m+i]=sqrt(sum);} else { L[j*m+i]=sum/L[i*m+i]; }
    }
    for(int j=i+1;j<m;j++) L[i*m+j]=0.0;
  }
  return true;
}

static void chol_solve(const double *L, const double *b, double *x, int m){
  double *y=palloc(sizeof(double)*(size_t)m);
  for(int i=0;i<m;i++){ double sum=b[i]; for(int k=0;k<i;k++) sum-=L[i*m+k]*y[k]; y[i]=sum/L[i*m+i]; }
  for(int i=m-1;i>=0;i--){ double sum=y[i]; for(int k=i+1;k<m;k++) sum-=L[k*m+i]*x[k]; x[i]=sum/L[i*m+i]; }
  pfree(y);
}

/* CV model */
typedef struct { int D; double q_accel_var; double r_meas_var; } CV_ModelCtx;

static bool cv_f(const double *x,const double*u,double dt,double*fx,double*F,void*v){(void)u;CV_ModelCtx*c=v;int D=c->D;int N=2*D;for(int i=0;i<D;i++){double p=x[2*i], v2=x[2*i+1]; fx[2*i]=p+dt*v2; fx[2*i+1]=v2;} mat_clear(F,N,N); for(int i=0;i<D;i++){F[(2*i)*N+(2*i)]=1; F[(2*i)*N+(2*i+1)]=dt; F[(2*i+1)*N+(2*i+1)]=1;} return true;}
static bool cv_h(const double *x,double*hx,double*H,void*v){CV_ModelCtx*c=v;int D=c->D;int N=2*D; for(int i=0;i<D;i++) hx[i]=x[2*i]; mat_clear(H,D,N); for(int i=0;i<D;i++) H[i*N+2*i]=1; return true;}
static bool cv_Q(double dt,double*Q,void*v){CV_ModelCtx*c=v;int D=c->D;int N=2*D; mat_clear(Q,N,N); double q=c->q_accel_var, dt2=dt*dt; double q11=q*(dt2*dt/3.0), q12=q*(dt2/2.0), q22=q*dt; for(int i=0;i<D;i++){int p=2*i, ve=2*i+1; Q[p*N+p]+=q11; Q[p*N+ve]+=q12; Q[ve*N+p]+=q12; Q[ve*N+ve]+=q22;} return true;}
static bool cv_R(double*R,void*v){CV_ModelCtx*c=v;int M=c->D; mat_clear(R,M,M); for(int i=0;i<M;i++) R[i*M+i]=c->r_meas_var; return true;}

bool tekf_make_cv_model(int D, TEkfModel *model){ if(D<=0||D>4||!model) return false; model->N=2*D; model->M=D; model->f=cv_f; model->h=cv_h; model->Q=cv_Q; model->R=cv_R; model->z_from_value=NULL; model->value_from_state=NULL; return true; }

/* GPS-like model */
static bool gps_f(const double *x,const double*u,double dt,double*fx,double*F,void*v){(void)u;const TEkfGpsCtx*c=v;int D=c->ndim;int N=2*D+(c->use_bias?1:0);for(int i=0;i<D;i++){double p=x[2*i], ve=x[2*i+1]; fx[2*i]=p+dt*ve; fx[2*i+1]=ve;} if(c->use_bias) fx[2*D]=x[2*D]; mat_clear(F,N,N); for(int i=0;i<D;i++){F[(2*i)*N+(2*i)]=1; F[(2*i)*N+(2*i+1)]=dt; F[(2*i+1)*N+(2*i+1)]=1;} if(c->use_bias) F[(2*D)*N+(2*D)]=1; return true;}
static bool gps_h(const double *x,double*hx,double*H,void*v){const TEkfGpsCtx*c=v;int D=c->ndim, M=c->M, N=2*D+(c->use_bias?1:0), b=(c->use_bias?2*D:-1); const double *a=c->anchors; mat_clear(H,M,N); for(int i=0;i<M;i++){ double s=0, diff[4]={0}; for(int d=0;d<D;d++){ double pd=x[2*d], sd=a[i*D+d], dd=pd-sd; diff[d]=dd; s+=dd*dd; } double r=sqrt(s); if(c->use_bias) r+=x[b]; hx[i]=r; double denom=fmax(r-(c->use_bias?x[b]:0),1e-12); for(int d=0;d<D;d++) H[i*N+2*d]=diff[d]/denom; if(c->use_bias) H[i*N+b]=1.0; } return true;}
static bool gps_Q(double dt,double*Q,void*v){const TEkfGpsCtx*c=v;int D=c->ndim, N=2*D+(c->use_bias?1:0); mat_clear(Q,N,N); double q=c->q_accel_var, dt2=dt*dt; double q11=q*(dt2*dt/3.0), q12=q*(dt2/2.0), q22=q*dt; for(int i=0;i<D;i++){int p=2*i,v2=2*i+1; Q[p*N+p]+=q11; Q[p*N+v2]+=q12; Q[v2*N+p]+=q12; Q[v2*N+v2]+=q22;} if(c->use_bias){int b=2*D; Q[b*N+b]+=c->q_bias_var*fmax(dt,1e-6);} return true;}
static bool gps_R(double*R,void*v){const TEkfGpsCtx*c=v;int M=c->M; mat_clear(R,M,M); for(int i=0;i<M;i++) R[i*M+i]=c->r_meas_var; return true;}
static bool gps_value_from_state(const double *x, meosType tt, Datum *out, void*v){const TEkfGpsCtx*c=v;int D=c->ndim; meosType bt=temptype_basetype(tt); int outD = (bt==T_DOUBLE3?3: bt==T_DOUBLE2?2: bt==T_FLOAT8?1:0); if(outD!=D) return false; double vec[4]={0}; for(int d=0;d<D;d++) vec[d]=x[2*d]; switch(bt){ case T_FLOAT8: *out=Float8GetDatum(vec[0]); return true; case T_DOUBLE2: { double2 *d2=palloc(sizeof(double2)); double2_set(vec[0],vec[1],d2); *out=Double2PGetDatum(d2); return true;} case T_DOUBLE3: { double3 *d3=palloc(sizeof(double3)); double3_set(vec[0],vec[1],vec[2],d3); *out=Double3PGetDatum(d3); return true;} default: return false; }}
bool tekf_make_gps_model(const TEkfGpsCtx *ctx, TEkfModel *model){ if(!ctx||!model) return false; if(!(ctx->ndim==2||ctx->ndim==3)) return false; if(ctx->M<=0||!ctx->anchors) return false; model->N=2*ctx->ndim+(ctx->use_bias?1:0); model->M=ctx->M; model->f=gps_f; model->h=gps_h; model->Q=gps_Q; model->R=gps_R; model->z_from_value=NULL; model->value_from_state=gps_value_from_state; return true; }

static int basetype_dim(meosType t){switch(t){case T_FLOAT8:return 1;case T_DOUBLE2:return 2;case T_DOUBLE3:return 3;case T_DOUBLE4:return 4;default:return 0;}}
static bool read_z_from_value(Datum v, meosType tt, int M, double *z){meosType bt=temptype_basetype(tt); int D=basetype_dim(bt); if(D!=M) return false; switch(bt){case T_FLOAT8: z[0]=DatumGetFloat8(v); return true; case T_DOUBLE2:{ const double2 *d=DatumGetDouble2P(v); z[0]=d->a; z[1]=d->b; return true;} case T_DOUBLE3:{ const double3 *d=DatumGetDouble3P(v); z[0]=d->a; z[1]=d->b; z[2]=d->c; return true;} case T_DOUBLE4:{ const double4 *d=DatumGetDouble4P(v); z[0]=d->a; z[1]=d->b; z[2]=d->c; z[3]=d->d; return true;} default: return false; }}
static Datum pack_value_from_vec(meosType tt, const double *v){meosType bt=temptype_basetype(tt); switch(bt){case T_FLOAT8: return Float8GetDatum(v[0]); case T_DOUBLE2:{ double2 *d=palloc(sizeof(double2)); double2_set(v[0],v[1],d); return Double2PGetDatum(d);} case T_DOUBLE3:{ double3 *d=palloc(sizeof(double3)); double3_set(v[0],v[1],v[2],d); return Double3PGetDatum(d);} case T_DOUBLE4:{ double4 *d=palloc(sizeof(double4)); double4_set(v[0],v[1],v[2],v[3],d); return Double4PGetDatum(d);} default: return (Datum)0; }}

static void build_R(const TEkfModel *m,const TEkfParams *p,double *R,void*ctx){int M=m->M; mat_clear(R,M,M); if(m->R && m->R(R,ctx)) return; if(p&&p->R_diag){for(int i=0;i<M;i++) R[i*M+i]=p->R_diag[i]; return;} for(int i=0;i<M;i++) R[i*M+i]=1e-6;}
static void build_Q(const TEkfModel *m,const TEkfParams *p,double dt,double *Q,void*ctx){int N=m->N; mat_clear(Q,N,N); if(m->Q && m->Q(dt,Q,ctx)) return; if(p&&p->Q_diag){for(int i=0;i<N;i++) Q[i*N+i]=p->Q_diag[i]; return;} for(int i=0;i<N;i++) Q[i*N+i]=1e-8;}
static void build_P0(const TEkfModel *m,const TEkfParams *p,double *P0){int N=m->N; if(p&&p->P0_diag){for(int i=0;i<N;i++) P0[i]=p->P0_diag[i]; return;} for(int i=0;i<N;i++) P0[i]=1.0;}
static void build_x0(const TEkfModel *m,const TEkfParams *p,const double *z,double *x0){int N=m->N, M=m->M; if(p&&p->x0){memcpy(x0,p->x0,sizeof(double)*(size_t)N); return;} for(int i=0;i<N;i++) x0[i]=0.0; int k=(M<N?M:N); for(int i=0;i<k;i++) x0[i]=z?z[i]:0.0;}

typedef struct{int N,M; double *x,*P,*fx,*F,*Q,*FP,*Ft,*hx,*H,*Ht,*y,*HP,*PHt,*S,*L,*K;} TEkfWs;
static TEkfWs *ws_create(int N,int M){TEkfWs*w=palloc0(sizeof(TEkfWs)); w->N=N;w->M=M; w->x=palloc(sizeof(double)*N); w->P=palloc(sizeof(double)*N*N); w->fx=palloc(sizeof(double)*N); w->F=palloc(sizeof(double)*N*N); w->Q=palloc(sizeof(double)*N*N); w->FP=palloc(sizeof(double)*N*N); w->Ft=palloc(sizeof(double)*N*N); w->hx=palloc(sizeof(double)*M); w->H=palloc(sizeof(double)*M*N); w->Ht=palloc(sizeof(double)*N*M); w->y=palloc(sizeof(double)*M); w->HP=palloc(sizeof(double)*M*N); w->PHt=palloc(sizeof(double)*N*M); w->S=palloc(sizeof(double)*M*M); w->L=palloc(sizeof(double)*M*M); w->K=palloc(sizeof(double)*N*M); return w;}
static void ws_free(TEkfWs*w){ if(!w) return; pfree(w->x); pfree(w->P); pfree(w->fx); pfree(w->F); pfree(w->Q); pfree(w->FP); pfree(w->Ft); pfree(w->hx); pfree(w->H); pfree(w->Ht); pfree(w->y); pfree(w->HP); pfree(w->PHt); pfree(w->S); pfree(w->L); pfree(w->K); pfree(w);} 

static inline double effective_dt(double obs, const TEkfParams *p){ double dt=obs; if(!(dt>0.0)) dt=(p&&p->default_dt>0.0)?p->default_dt:0.0; return dt;}

static Temporal * ekf_clean_instant(const TInstant *inst){ return (Temporal *) tinstant_copy(inst);} 

static TSequence * ekf_clean_sequence(const TSequence *seq, const TEkfModel *model, const TEkfParams *params_in, void *ctx, int *removed){
  const int N=model->N, M=model->M; int n=seq->count; if(n<=0) return NULL; TEkfParams defp, *params=(TEkfParams *)&defp; if(params_in) params=(TEkfParams *)params_in; else { memset(&defp,0,sizeof(defp)); defp.default_dt=1.0; defp.gate_sigma=3.5; defp.fill_estimates=false; }
  TEkfWs *w=ws_create(N,M); double *P0=palloc(sizeof(double)*N); build_P0(model,params,P0); for(int i=0;i<N;i++) for(int j=0;j<N;j++) w->P[i*N+j]=(i==j)?P0[i]:0.0;
  TInstant **out=palloc(sizeof(TInstant *)*n); int outc=0; bool has=false; TimestampTz prev=0;
  for(int i=0;i<n;i++){
    const TInstant *inst=TSEQUENCE_INST_N(seq,i); Datum val=tinstant_value_p(inst); double *z=palloc(sizeof(double)*M); bool have=false; if(model->z_from_value) have=model->z_from_value(val,inst->temptype,z,ctx); else have=read_z_from_value(val,inst->temptype,M,z);
    if(!has){ if(!have){ prev=inst->t; pfree(z); continue; } build_x0(model,params,z,w->x); /* emit initial */
      if(model->value_from_state){ Datum outv; if(model->value_from_state(w->x,inst->temptype,&outv,ctx)){ out[outc++]=tinstant_make(outv,inst->temptype,inst->t); DATUM_FREE(outv,temptype_basetype(inst->temptype)); } else { Datum outv2=pack_value_from_vec(inst->temptype,z); out[outc++]=tinstant_make(outv2,inst->temptype,inst->t); DATUM_FREE(outv2,temptype_basetype(inst->temptype)); } }
      else { Datum outv=pack_value_from_vec(inst->temptype,z); out[outc++]=tinstant_make(outv,inst->temptype,inst->t); DATUM_FREE(outv,temptype_basetype(inst->temptype)); }
      has=true; prev=inst->t; pfree(z); continue; }
    double dt=effective_dt(((double)(inst->t - prev))/1000000.0, params);
    build_Q(model,params,dt,w->Q,ctx); model->f(w->x,NULL,dt,w->fx,w->F,ctx); memcpy(w->x,w->fx,sizeof(double)*N); mat_mul(w->F,w->P,w->FP,N,N,N); mat_transpose(w->F,w->Ft,N,N); mat_mul(w->FP,w->Ft,w->P,N,N,N); mat_add(w->P,w->P,w->Q,N,N);
    if(!have){ if(params->fill_estimates){ Datum outv; bool ok=false; if(model->value_from_state) ok = model->value_from_state(w->x, inst->temptype, &outv, ctx); if(!ok){ model->h(w->x,w->hx,w->H,ctx); outv = pack_value_from_vec(inst->temptype,w->hx); } out[outc++]=tinstant_make(outv,inst->temptype,inst->t); DATUM_FREE(outv,temptype_basetype(inst->temptype)); } prev=inst->t; pfree(z); continue; }
    build_R(model,params,w->S,ctx); model->h(w->x,w->hx,w->H,ctx); for(int ii=0;ii<M;ii++) w->y[ii]=z[ii]-w->hx[ii]; mat_mul(w->H,w->P,w->HP,M,N,N); mat_transpose(w->H,w->Ht,M,N); mat_mul(w->HP,w->Ht,w->S,M,N,M); double *Rt=palloc(sizeof(double)*M*M); build_R(model,params,Rt,ctx); mat_add(w->S,w->S,Rt,M,M); if(!chol_decompose(w->S,w->L,M)){ for(int ii=0;ii<M;ii++) w->S[ii*M+ii]+=1e-9; if(!chol_decompose(w->S,w->L,M)){ if(params->fill_estimates){ Datum outv; bool ok=false; if(model->value_from_state) ok = model->value_from_state(w->x, inst->temptype, &outv, ctx); if(!ok){ outv = pack_value_from_vec(inst->temptype, w->hx); } out[outc++]=tinstant_make(outv,inst->temptype,inst->t); DATUM_FREE(outv,temptype_basetype(inst->temptype)); } else if(removed) (*removed)++; prev=inst->t; pfree(Rt); pfree(z); continue; }} double *tmpy=palloc(sizeof(double)*M); chol_solve(w->L,w->y,tmpy,M); double d2=0; for(int ii=0;ii<M;ii++) d2+=w->y[ii]*tmpy[ii]; pfree(tmpy); bool reject=(params->gate_sigma>0.0)&&(d2>params->gate_sigma*params->gate_sigma); if(reject){ if(params->fill_estimates){ Datum outv; bool ok=false; if(model->value_from_state) ok = model->value_from_state(w->x, inst->temptype, &outv, ctx); if(!ok){ outv = pack_value_from_vec(inst->temptype, w->hx); } out[outc++]=tinstant_make(outv,inst->temptype,inst->t); DATUM_FREE(outv,temptype_basetype(inst->temptype)); } else if(removed) (*removed)++; prev=inst->t; pfree(Rt); pfree(z); continue; }
    mat_transpose(w->HP,w->PHt,M,N); for(int j=0;j<N;j++){ double *rhs=palloc(sizeof(double)*M), *col=palloc(sizeof(double)*M); for(int ii=0;ii<M;ii++) rhs[ii]=w->PHt[j*M+ii]; chol_solve(w->L,rhs,col,M); for(int ii=0;ii<M;ii++) w->K[j*M+ii]=col[ii]; pfree(rhs); pfree(col);} for(int ii=0;ii<N;ii++){ double s=0; for(int jj=0;jj<M;jj++) s+=w->K[ii*M+jj]*w->y[jj]; w->x[ii]+=s; } double *KHP=palloc(sizeof(double)*N*N); mat_mul(w->K,w->HP,KHP,N,M,N); mat_sub(w->P,w->P,KHP,N,N); pfree(KHP);
    { Datum outv; bool ok=false; if(model->value_from_state) ok = model->value_from_state(w->x, inst->temptype, &outv, ctx); if(!ok){ model->h(w->x,w->hx,w->H,ctx); outv = pack_value_from_vec(inst->temptype,w->hx); } out[outc++]=tinstant_make(outv,inst->temptype,inst->t); DATUM_FREE(outv,temptype_basetype(inst->temptype)); } prev=inst->t; pfree(Rt); pfree(z);
  }
  TSequence *res=NULL; if(outc>0) res=tsequence_make((const TInstant **)out,outc, seq->period.lower_inc, seq->period.upper_inc, MEOS_FLAGS_GET_INTERP(seq->flags), NORMALIZE_NO); pfree_array((void**)out,outc); pfree(P0); ws_free(w); return res;
}

static Temporal * ekf_dispatch(const Temporal *t, const TEkfModel *m, const TEkfParams *p, void *ctx, int *removed)
{
  switch (t->subtype)
  {
    case TINSTANT:
      return ekf_clean_instant((const TInstant *) t);
    case TSEQUENCE:
      return (Temporal *) ekf_clean_sequence((const TSequence *) t, m, p, ctx, removed);
    case TSEQUENCESET:
    {
      const TSequenceSet *ss = (const TSequenceSet *) t;
      TSequence **arr = palloc(sizeof(TSequence *) * ss->count);
      int k = 0;
      for (int i = 0; i < ss->count; i++)
      {
        TSequence *s = ekf_clean_sequence(TSEQUENCESET_SEQ_N(ss, i), m, p, ctx, removed);
        if (s) arr[k++] = s;
      }
      return (Temporal *) tsequenceset_make_free(arr, k, NORMALIZE_NO);
    }
    default:
      return NULL;
  }
}

Temporal * temporal_ekf_clean(const Temporal *temp, const TEkfModel *model, const TEkfParams *params, void *ctx, int *removed)
{
  VALIDATE_NOT_NULL(temp, NULL);
  VALIDATE_NOT_NULL(model, NULL);
  if (removed) *removed = 0;
  if (!(temp->subtype == TINSTANT || temp->subtype == TSEQUENCE || temp->subtype == TSEQUENCESET))
  {
    meos_error(ERROR, MEOS_ERR_INVALID_ARG, "Unsupported temporal subtype for EKF");
    return NULL;
  }
  if (model->N <= 0 || model->M <= 0 || !model->f || !model->h)
  {
    meos_error(ERROR, MEOS_ERR_INVALID_ARG, "Invalid EKF model (dims or callbacks)");
    return NULL;
  }
  return ekf_dispatch(temp, model, params, ctx, removed);
}

#endif
