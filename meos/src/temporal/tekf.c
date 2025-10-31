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

/* Basic matrix operations */

/* clear matrix */
static inline void mat_clear(double *A, int m, int n){ 
  memset(A,0,sizeof(double)*(size_t)m*(size_t)n);
} 

/* copy matrix */
static inline void mat_copy(double *d,const double*s,int m,int n){ 
  memcpy(d,s,sizeof(double)*(size_t)m*(size_t)n);
}

/* add, subtract, multiply, transpose */
static inline void mat_add(double *C,const double*A,const double*B,int m,int n){ 
  int mn=m*n;
  for(int i=0;i<mn;i++)
    C[i]=A[i]+B[i];
}
static inline void mat_sub(double *C,const double*A,const double*B,int m,int n){ 
  int mn=m*n;
  for(int i=0;i<mn;i++)
    C[i]=A[i]-B[i];
}
static inline void mat_mul(const double*A,const double*B,double*C,int arows,int acols,int bcols){ 
  for(int i=0;i<arows;i++){
    for(int j=0;j<bcols;j++){
      double s=0;
      for(int k=0;k<acols;k++)
        s+=A[i*acols+k]*B[k*bcols+j];
      C[i*bcols+j]=s;
    }
  }
}
static inline void mat_transpose(const double*A,double*At,int m,int n){ 
  for(int i=0;i<m;i++)
    for(int j=0;j<n;j++)
      At[j*m+i]=A[i*n+j];
}

/* Cholesky decomposition and solver for positive-definite matrices */
static bool chol_decompose(const double *A, double *L, int m){
  mat_copy(L,A,m,m);
  for(int i=0;i<m;i++){
    for(int j=i;j<m;j++){
      double sum=L[i*m+j];
      for(int k=0;k<i;k++) 
        sum-=L[i*m+k]*L[j*m+k];
      if(i==j){ 
        if(sum<=0.0) return false; 
        L[i*m+i]=sqrt(sum);
      } 
      else { 
        L[j*m+i]=sum/L[i*m+i]; 
      }
    }
    for(int j=i+1;j<m;j++) 
      L[i*m+j]=0.0;
  }
  return true;
}

/* Solve L*y = b for y */
static void chol_solve(const double *L, const double *b, double *x, int m){
  double *y=palloc(sizeof(double)*(size_t)m);
  for(int i=0;i<m;i++){ 
    double sum=b[i]; 
    for(int k=0;k<i;k++) 
      sum-=L[i*m+k]*y[k]; 
    y[i]=sum/L[i*m+i]; 
  }
  for(int i=m-1;i>=0;i--){ 
    double sum=y[i]; 
    for(int k=i+1;k<m;k++) 
      sum-=L[k*m+i]*x[k]; 
    x[i]=sum/L[i*m+i]; 
  }
  pfree(y);
}


/* f is the state transition function */
static bool cv_f(const double *x,const double*u,double dt,double*fx,double*F,void*v){
  (void)u;
  CV_ModelCtx*c=v;
  int D=c->D;
  int N=2*D;
  for(int i=0;i<D;i++){
    double p=x[2*i], v2=x[2*i+1]; 
    fx[2*i]=p+dt*v2; 
    fx[2*i+1]=v2;
  } 
  mat_clear(F,N,N); 
  for(int i=0;i<D;i++){
    F[(2*i)*N+(2*i)]=1; 
    F[(2*i)*N+(2*i+1)]=dt; 
    F[(2*i+1)*N+(2*i+1)]=1;
  } 
  return true;
}

/* h is the measurement function */
static bool cv_h(const double *x,double*hx,double*H,void*v){
  CV_ModelCtx*c=v;
  int D=c->D;
  int N=2*D;
  for(int i=0;i<D;i++) 
    hx[i]=x[2*i];
  mat_clear(H,D,N);
  for(int i=0;i<D;i++) 
    H[i*N+2*i]=1;
  return true;
}
/* Q is the process noise covariance */
static bool cv_Q(double dt,double*Q,void*v){
  CV_ModelCtx*c=v;
  int D=c->D;
  int N=2*D;
  mat_clear(Q,N,N);
  double q=c->q_accel_var, dt2=dt*dt;
  double q11=q*(dt2*dt/3.0), q12=q*(dt2/2.0), q22=q*dt;
  for(int i=0;i<D;i++){
    int p=2*i, ve=2*i+1;
    Q[p*N+p]+=q11;
    Q[p*N+ve]+=q12;
    Q[ve*N+p]+=q12;
    Q[ve*N+ve]+=q22;
  }
  return true;
}

/* R is the measurement noise covariance */
static bool cv_R(double*R,void*v){
  CV_ModelCtx*c=v;
  int M=c->D;
  mat_clear(R,M,M);
  for(int i=0;i<M;i++)  
    R[i*M+i]=c->r_meas_var;
  return true;
}

/* Builder for CV model */
bool tekf_make_cv_model(int D, TEkfModel *model){
  if(D<=0||D>4||!model) 
    return false;
  model->N=2*D;
  model->M=D;
  model->f=cv_f;
  model->h=cv_h;
  model->Q=cv_Q;
  model->R=cv_R;
  model->z_from_value=NULL;
  model->value_from_state=NULL;
  return true;
}

/* GPS-like model */
/* Define context for GPS-like model */
/* f is the state transition function */
static bool gps_f(const double *x,const double*u,double dt,double*fx,double*F,void*v){
  (void)u;
  const TEkfGpsCtx*c=v;
  int D=c->ndim;
  int N=2*D+(c->use_bias?1:0);
  for(int i=0;i<D;i++){
    double p=x[2*i], ve=x[2*i+1]; 
    fx[2*i]=p+dt*ve; 
    fx[2*i+1]=ve;
  } 
  if(c->use_bias) 
    fx[2*D]=x[2*D]; 
  mat_clear(F,N,N); 
  for(int i=0;i<D;i++){
    F[(2*i)*N+(2*i)]=1; 
    F[(2*i)*N+(2*i+1)]=dt; 
    F[(2*i+1)*N+(2*i+1)]=1;
  } 
  if(c->use_bias) F[(2*D)*N+(2*D)]=1; 
  return true;
}

/* h is the measurement function */
static bool gps_h(const double *x,double*hx,double*H,void*v){
  const TEkfGpsCtx*c=v;
  int D=c->ndim, M=c->M, N=2*D+(c->use_bias?1:0), b=(c->use_bias?2*D:-1); 
  const double *a=c->anchors; 
  mat_clear(H,M,N); 
  for(int i=0;i<M;i++){ 
    double s=0, diff[4]={0}; 
    for(int d=0;d<D;d++){ 
      double pd=x[2*d], sd=a[i*D+d], dd=pd-sd; 
      diff[d]=dd; 
      s+=dd*dd; 
    } 
    double r=sqrt(s); 
    if(c->use_bias) r+=x[b]; 
    hx[i]=r; 
    double denom=fmax(r-(c->use_bias?x[b]:0),1e-12); 
    for(int d=0;d<D;d++) 
      H[i*N+2*d]=diff[d]/denom; 
    if(c->use_bias) 
      H[i*N+b]=1.0; 
  } 
  return true;
}

/* Q is the process noise covariance (CV per axis + optional bias RW) */
bool gps_Q(double dt, double *Q, void *v)
{
  const TEkfGpsCtx *c = v;
  int D = c->ndim;
  int N = 2*D + (c->use_bias ? 1 : 0);
  mat_clear(Q, N, N);

  /* Constant-acceleration model per axis */
  double q = c->q_accel_var;
  double dt2 = dt * dt;
  double q11 = q * (dt2*dt/3.0);
  double q12 = q * (dt2/2.0);
  double q22 = q * dt;
  for (int i = 0; i < D; i++)
  {
    int p = 2*i, vix = 2*i+1;
    Q[p*N + p]     += q11;
    Q[p*N + vix]   += q12;
    Q[vix*N + p]   += q12;
    Q[vix*N + vix] += q22;
  }

  /* Optional bias random walk */
  if (c->use_bias)
  {
    int b = 2*D;
    double qb = c->q_bias_var;
    if (qb > 0.0 && dt > 0.0)
      Q[b*N + b] += qb * dt;
  }
  return true;
}

/* R is the measurement noise covariance */
static bool gps_R(double*R,void*v){
  const TEkfGpsCtx*c=v;
  int M=c->M; 
  mat_clear(R,M,M); 
  for(int i=0;i<M;i++) 
    R[i*M+i]=c->r_meas_var; 
  return true;
}
/* value_from_state: extract position (and bias if applicable) from state */
/* return false if state is invalid as NaN */
static bool gps_value_from_state(const double *x, meosType tt, Datum *out, void*v){
  const TEkfGpsCtx*c=v;
  int D=c->ndim; 
  meosType bt=temptype_basetype(tt); 
  int outD = (bt==T_DOUBLE3?3: bt==T_DOUBLE2?2: bt==T_FLOAT8?1:0); 
  if(outD!=D) 
    return false; 

  double vec[4]={0}; 
  for(int d=0;d<D;d++) 
    vec[d]=x[2*d]; 

  switch(bt){ 
    case T_FLOAT8: *out=Float8GetDatum(vec[0]); 
    return true; 
    case T_DOUBLE2: { 
      double2 *d2=palloc(sizeof(double2)); 
      double2_set(vec[0],vec[1],d2); 
      *out=Double2PGetDatum(d2); 
      return true;
    } case T_DOUBLE3: { 
      double3 *d3=palloc(sizeof(double3)); 
      double3_set(vec[0],vec[1],vec[2],d3); 
      *out=Double3PGetDatum(d3);
      return true;
    } 
    default: return false;
  }
  
}

/* Builder for GPS-like model */
bool tekf_make_gps_model(const TEkfGpsCtx *ctx, TEkfModel *model){ 
  if(!ctx||!model) return false; 
  if(!(ctx->ndim==2||ctx->ndim==3)) 
    return false; 
  if(ctx->M<=0||!ctx->anchors) 
    return false; 
  model->N=2*ctx->ndim+(ctx->use_bias?1:0); 
  model->M=ctx->M; 
  model->f=gps_f; 
  model->h=gps_h; 
  model->Q=gps_Q; 
  model->R=gps_R; 
  model->z_from_value=NULL; 
  model->value_from_state=gps_value_from_state; 
  return true; 
}

/* Helpers to read/write measurement values */
static int basetype_dim(meosType t){
  switch(t){
    case T_FLOAT8:return 1;
    case T_DOUBLE2:return 2;
    case T_DOUBLE3:return 3;
    case T_DOUBLE4:return 4;
    default:return 0;
  }
}
/* z is of length M */
static bool read_z_from_value(Datum v, meosType tt, int M, double *z){
  meosType bt=temptype_basetype(tt);
  int D=basetype_dim(bt);
  if(D!=M) return false;
  switch(bt){
    case T_FLOAT8: z[0]=DatumGetFloat8(v); return true;
    case T_DOUBLE2:{ 
      const double2 *d=DatumGetDouble2P(v); z[0]=d->a; z[1]=d->b; return true;
    }
    case T_DOUBLE3:{ 
      const double3 *d=DatumGetDouble3P(v); z[0]=d->a; z[1]=d->b; z[2]=d->c; return true;
    }
    case T_DOUBLE4:{ 
      const double4 *d=DatumGetDouble4P(v); z[0]=d->a; z[1]=d->b; z[2]=d->c; z[3]=d->d; return true;
    }
    default: return false;
  }
}

/* pack value from vector */
static Datum pack_value_from_vec(meosType tt, const double *v){
  meosType bt=temptype_basetype(tt);
  switch(bt){
    case T_FLOAT8: return Float8GetDatum(v[0]);
    case T_DOUBLE2:{ 
      double2 *d=palloc(sizeof(double2)); double2_set(v[0],v[1],d); return Double2PGetDatum(d);
    }
    case T_DOUBLE3:{ 
      double3 *d=palloc(sizeof(double3)); double3_set(v[0],v[1],v[2],d); return Double3PGetDatum(d);
    }
    case T_DOUBLE4:{ 
      double4 *d=palloc(sizeof(double4)); double4_set(v[0],v[1],v[2],v[3],d); return Double4PGetDatum(d);
    }
    default: return (Datum)0;
  }
}

/* Builders for initial matrices/vectors */
static void build_R(const TEkfModel *m,const TEkfParams *p,double *R,void*ctx){
  int M=m->M;
  mat_clear(R,M,M);
  if(m->R && m->R(R,ctx)) 
    return;
  if(p&&p->R_diag){
    for(int i=0;i<M;i++) R[i*M+i]=p->R_diag[i];
    return;
  }
  for(int i=0;i<M;i++)
    R[i*M+i]=1e-6;
}

/* Q is the process noise covariance */
static void build_Q(const TEkfModel *m,const TEkfParams *p,double dt,double *Q,void*ctx){
  int N=m->N;
  mat_clear(Q,N,N);
  if(m->Q && m->Q(dt,Q,ctx)) 
    return;
  if(p&&p->Q_diag){
    for(int i=0;i<N;i++) Q[i*N+i]=p->Q_diag[i];
    return;
  }
  for(int i=0;i<N;i++) 
    Q[i*N+i]=1e-8;
}

/* P0 is the initial error covariance */
static void build_P0(const TEkfModel *m,const TEkfParams *p,double *P0){
  int N=m->N;
  if(p&&p->P0_diag){
    for(int i=0;i<N;i++) 
      P0[i]=p->P0_diag[i]; 
    return;
  }
  for(int i=0;i<N;i++) 
    P0[i]=1.0;
}

/* x0 is the initial state */
static void build_x0(const TEkfModel *m,const TEkfParams *p,const double *z,double *x0){
  int N=m->N, M=m->M;
  if(p&&p->x0){
    memcpy(x0,p->x0,sizeof(double)*(size_t)N); 
    return;
  }
  for(int i=0;i<N;i++) 
    x0[i]=0.0;
  int k=(M<N?M:N);
  for(int i=0;i<k;i++) 
    x0[i]=z?z[i]:0.0;
}

/* Rough trilateration initializer for GPS ranges (least-squares), sets pos and bias */
static bool gps_init_from_ranges(const TEkfGpsCtx *ctx, const double *z, double *x0)
{
  if (!ctx || !z || !x0) return false;
  const int D = ctx->ndim;
  const int M = ctx->M;
  if (M < D + 1) 
    return false; /* need at least D+1 anchors */
  /* Build linear system A p = b from differences to anchor 0 */
  const double *a0 = &ctx->anchors[0];
  double *A = palloc0(sizeof(double) * (size_t)(M-1) * (size_t)D);
  double *b = palloc0(sizeof(double) * (size_t)(M-1));
  /* Build linear system */
  for (int i = 1; i < M; i++)
  {
    const double *ai = &ctx->anchors[i * D];
    for (int d = 0; d < D; d++)
      A[(i-1)*D + d] = ai[d] - a0[d];
    double ai2 = 0.0, a0p = 0.0;
    for (int d = 0; d < D; d++) { 
      ai2 += ai[d]*ai[d]; 
      a0p += a0[d]*a0[d]; 
    }
    b[i-1] = 0.5 * ( (z[0]*z[0] - z[i]*z[i]) - (a0p - ai2) );
  }
  /* Solve least squares via normal equations: (A^T A) p = A^T b */
  double *At = palloc(sizeof(double) * (size_t)D * (size_t)(M-1));
  double *AtA = palloc(sizeof(double) * (size_t)D * (size_t)D);
  double *Atb = palloc(sizeof(double) * (size_t)D);
  mat_transpose(A, At, M-1, D);
  mat_mul(At, A, AtA, D, M-1, D);
  /* Cholesky solve AtA p = At b */
  for (int d = 0; d < D; d++)
  {
    double sum = 0.0; 
    for (int i = 0; i < M-1; i++) 
      sum += At[d*(M-1) + i] * b[i];
    Atb[d] = sum;
  }
  double *L = palloc(sizeof(double) * (size_t)D * (size_t)D);
  if (!chol_decompose(AtA, L, D)) { /* regularize */
    for (int d = 0; d < D; d++) 
      AtA[d*D + d] += 1e-6;
    if (!chol_decompose(AtA, L, D)) { 
      pfree(A); 
      pfree(b); 
      pfree(At); 
      pfree(AtA); 
      pfree(Atb); 
      pfree(L); 
      return false; 
    }
  }
  double *pvec = palloc(sizeof(double) * (size_t)D);
  /* solution in pvec */
  chol_solve(L, Atb, pvec, D);
  for (int d = 0; d < D; d++) 
    x0[2*d] = pvec[d], x0[2*d+1] = 0.0;
  if (ctx->use_bias)
  {
    /* bias ~ z0 - ||p - a0|| */
    double diff2 = 0.0; 
    for (int d = 0; d < D; d++) { 
      double dd = pvec[d] - a0[d]; diff2 += dd*dd; 
    }
    x0[2*D] = z[0] - sqrt(diff2);
  }
  pfree(A); pfree(b); pfree(At); pfree(AtA); pfree(Atb); pfree(L); pfree(pvec);
  return true;
}


/* Workspace management */
static TEkfWs *ws_create(int N,int M){
  TEkfWs*w=palloc0(sizeof(TEkfWs)); 
  w->N=N;w->M=M; 
  w->x=palloc(sizeof(double)*N); 
  w->P=palloc(sizeof(double)*N*N); 
  w->fx=palloc(sizeof(double)*N); 
  w->F=palloc(sizeof(double)*N*N); 
  w->Q=palloc(sizeof(double)*N*N); 
  w->FP=palloc(sizeof(double)*N*N); 
  w->Ft=palloc(sizeof(double)*N*N); 
  w->hx=palloc(sizeof(double)*M); 
  w->H=palloc(sizeof(double)*M*N); 
  w->Ht=palloc(sizeof(double)*N*M);
  w->y=palloc(sizeof(double)*M); 
  w->HP=palloc(sizeof(double)*M*N); 
  w->PHt=palloc(sizeof(double)*N*M); 
  w->S=palloc(sizeof(double)*M*M); 
  w->L=palloc(sizeof(double)*M*M); 
  w->K=palloc(sizeof(double)*N*M); return w;
}

/* Free workspace */
static void ws_free(TEkfWs*w){ 
  if(!w) return; 
  pfree(w->x); 
  pfree(w->P); 
  pfree(w->fx); 
  pfree(w->F); 
  pfree(w->Q); 
  pfree(w->FP); 
  pfree(w->Ft); 
  pfree(w->hx); 
  pfree(w->H); 
  pfree(w->Ht); 
  pfree(w->y); 
  pfree(w->HP); 
  pfree(w->PHt); 
  pfree(w->S); 
  pfree(w->L); 
  pfree(w->K); 
  pfree(w);} 

/* Determine dt */
static inline double effective_dt(double obs, const TEkfParams *p){ 
  double dt=obs; 
  if(!(dt>0.0)) 
    dt=(p&&p->default_dt>0.0)?p->default_dt:0.0; 
  return dt;
}

/* EKF cleaning of instant: just return a copy */
static Temporal * ekf_clean_instant(const TInstant *inst){ 
  return (Temporal *) tinstant_copy(inst);
}

/* EKF cleaning of sequence */
/* the algorithm takes one temporal sequence 
and runs an Extended Kalman Filter (EKF) over it */
static TSequence * ekf_clean_sequence(const TSequence *seq, const TEkfModel *model, const TEkfParams *params_in, void *ctx, int *removed){
  const int N=model->N, M=model->M; 
  int n=seq->count; 

  /* if no observations */
  if(n <= 0) return NULL; 
  /* For very short sequences, just return a copy */
  if (n < 3) return tsequence_copy(seq);

  /* Prepare parameters */
  TEkfParams defp, *params=(TEkfParams *)&defp; 
  if(params_in) 
    params=(TEkfParams *)params_in; 
  else { 
    memset(&defp,0,sizeof(defp)); 
    defp.default_dt=1.0; 
    defp.gate_sigma=3.5; 
    defp.fill_estimates=false; 
  }
  /* Create workspace and initialize, meaning allocate memory for all matrices */
  TEkfWs *w=ws_create(N,M); 
  double *P0=palloc(sizeof(double)*N); 
  build_P0(model,params,P0); 
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      w->P[i*N+j]=(i==j)?P0[i]:0.0;

  /* Output instants */
  TInstant **out=palloc(sizeof(TInstant *)*n); 
  int outc=0; 
  bool has=false; 
  TimestampTz prev=0;

  /* Main EKF loop */
  for(int i=0;i<n;i++){
    const TInstant *inst=TSEQUENCE_INST_N(seq,i); 
    Datum val=tinstant_value_p(inst); 
    double *z=palloc(sizeof(double)*M); 
    bool have=false; 

    /* Extract measurement from value */
    /* if the model provides a way to extract measurements from the value 
    else use a default method that reads the value into the measurement vector */
    if(model->z_from_value) 
      have=model->z_from_value(val,inst->temptype,z,ctx); 
    else 
      have=read_z_from_value(val,inst->temptype,M,z);

    /* Initialization */
    if(!has){ 
      if(!have){ 
        prev=inst->t; 
        pfree(z); 
        continue; 
      }
      /* GPS-specific init if applicable */
      if (model->h == gps_h) { 
        (void) gps_init_from_ranges((const TEkfGpsCtx *) ctx, z, w->x); 
      }
      else { 
        build_x0(model,params,z,w->x); 
      }
      /* emit initial instant, if there is a valid measurement */
      if (model->value_from_state) {
        Datum outv;
        if (model->value_from_state(w->x, inst->temptype, &outv, ctx)) {
          out[outc++] = tinstant_make(outv, inst->temptype, inst->t);
          DATUM_FREE(outv, temptype_basetype(inst->temptype));
        } else {
          Datum outv2 = pack_value_from_vec(inst->temptype, z);
          out[outc++] = tinstant_make(outv2, inst->temptype, inst->t);
          DATUM_FREE(outv2, temptype_basetype(inst->temptype));
        }
      } /* else just copy measurement */
      else {
        Datum outv = pack_value_from_vec(inst->temptype, z);
        out[outc++] = tinstant_make(outv, inst->temptype, inst->t);
        DATUM_FREE(outv, temptype_basetype(inst->temptype));
      }
      has=true; prev=inst->t; pfree(z); continue; 
    }
    /* Prediction step */
    double dt=effective_dt(((double)(inst->t - prev))/1000000.0, params);
    /* Predict state and error covariance */
    build_Q(model,params,dt,w->Q,ctx); 
    /* x = f(x) */
    model->f(w->x,NULL,dt,w->fx,w->F,ctx); 
    /* Jacobian */
    memcpy(w->x,w->fx,sizeof(double)*N); 
    /* P = F*P*F^T + Q */
    mat_mul(w->F,w->P,w->FP,N,N,N); 
    /* F^T */
    mat_transpose(w->F,w->Ft,N,N); 
    /* P = F*P*F^T + Q */
    mat_mul(w->FP,w->Ft,w->P,N,N,N); 
    /* + Q */
    mat_add(w->P,w->P,w->Q,N,N);

    /* Update step */
    if(!have){ 
      /* if the model provides a way to extract measurements from the value 
      else use a default method that reads the value into the measurement vector */
      if(params->fill_estimates){ 
        Datum outv; 
        bool ok=false; 
        if(model->value_from_state) 
          ok = model->value_from_state(w->x, inst->temptype, &outv, ctx); 
        if(!ok){ 
          model->h(w->x,w->hx,w->H,ctx); 
          outv = pack_value_from_vec(inst->temptype,w->hx); 
        } 
        out[outc++]=tinstant_make(outv,inst->temptype,inst->t); 
        DATUM_FREE(outv,temptype_basetype(inst->temptype)); 
      } prev=inst->t; pfree(z); 
      continue; 
    }
    build_R(model,params,w->S,ctx); 
    model->h(w->x,w->hx,w->H,ctx); 
    for(int ii=0;ii<M;ii++) 
      w->y[ii]=z[ii]-w->hx[ii]; 

    mat_mul(w->H,w->P,w->HP,M,N,N); 
    mat_transpose(w->H,w->Ht,M,N); 
    mat_mul(w->HP,w->Ht,w->S,M,N,M); 

    double *Rt=palloc(sizeof(double)*M*M); 
    build_R(model,params,Rt,ctx); 
    mat_add(w->S,w->S,Rt,M,M); 

    /* gating and outlier rejection */
    if(!chol_decompose(w->S,w->L,M)){ 
      for(int ii=0;ii<M;ii++) 
        w->S[ii*M+ii]+=1e-9; 
      if(!chol_decompose(w->S,w->L,M)){ 
        if(params->fill_estimates){ 
          Datum outv; 
          bool ok=false; 
          if(model->value_from_state) ok = model->value_from_state(w->x, inst->temptype, &outv, ctx); 
          if(!ok){ 
            outv = pack_value_from_vec(inst->temptype, w->hx); 
          } 
          out[outc++]=tinstant_make(outv,inst->temptype,inst->t); 
          DATUM_FREE(outv,temptype_basetype(inst->temptype)); 
        } 
        else if(removed) (*removed)++; 
        prev=inst->t; 
        pfree(Rt); 
        pfree(z); 
        continue; 
      }} 
      /* Mahalanobis distance */
      double *tmpy=palloc(sizeof(double)*M); 
      chol_solve(w->L,w->y,tmpy,M); 
      double d2=0; 
      for(int ii=0;ii<M;ii++) 
        d2+=w->y[ii]*tmpy[ii]; 
      pfree(tmpy); 
      bool reject=(params->gate_sigma>0.0)&&(d2>params->gate_sigma*params->gate_sigma); 
      if(reject){ 
        if(params->fill_estimates){ 
          Datum outv; 
          bool ok=false; 
          if(model->value_from_state) ok = model->value_from_state(w->x, inst->temptype, &outv, ctx); 
          if(!ok){ 
            outv = pack_value_from_vec(inst->temptype, w->hx); 
          } 
          out[outc++]=tinstant_make(outv,inst->temptype,inst->t); 
          DATUM_FREE(outv,temptype_basetype(inst->temptype)); 
        } 
        else if(removed) 
          (*removed)++; 
        prev=inst->t; 
        pfree(Rt); 
        pfree(z); 
        continue; 
      }
    mat_transpose(w->HP,w->PHt,M,N); 
    /* Compute Kalman gain K = P*H^T*S^-1 via Cholesky solve */
    for(int j=0;j<N;j++){ 
      double *rhs=palloc(sizeof(double)*M), *col=palloc(sizeof(double)*M); 
      for(int ii=0;ii<M;ii++) 
        rhs[ii]=w->PHt[j*M+ii]; 
      chol_solve(w->L,rhs,col,M); 
      for(int ii=0;ii<M;ii++) 
      w->K[j*M+ii]=col[ii]; 
      pfree(rhs); pfree(col);
    } 
    /* Update state x = x + K*y */
    for(int ii=0;ii<N;ii++){ 
      double s=0; 
      for(int jj=0;jj<M;jj++) 
        s+=w->K[ii*M+jj]*w->y[jj]; 
        w->x[ii]+=s; 
      } 
    /* Update error covariance P = P - K*H*P */
      double *KHP=palloc(sizeof(double)*N*N); 
      mat_mul(w->K,w->HP,KHP,N,M,N); 
      mat_sub(w->P,w->P,KHP,N,N); 
      pfree(KHP);
      /* Emit cleaned instant */
      Datum outv; 
      bool ok=false; 
      if(model->value_from_state) 
        ok = model->value_from_state(w->x, inst->temptype, &outv, ctx); 
      if(!ok){ 
        model->h(w->x,w->hx,w->H,ctx); 
        outv = pack_value_from_vec(inst->temptype,w->hx); 
      } 
      out[outc++]=tinstant_make(outv,inst->temptype,inst->t); 
      DATUM_FREE(outv,temptype_basetype(inst->temptype)); 
    prev=inst->t; 
    pfree(Rt); 
    pfree(z);
  }
  TSequence *res=NULL; 
  if(outc>0) 
    res=tsequence_make((const TInstant **)out,outc, seq->period.lower_inc, seq->period.upper_inc, MEOS_FLAGS_GET_INTERP(seq->flags), NORMALIZE_NO); 
  pfree_array((void**)out,outc); pfree(P0); ws_free(w); return res;
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
