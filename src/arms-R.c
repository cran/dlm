/* adaptive rejection metropolis sampling                                 */
/* original C code by Wally Gilks <Wally.Gilks@mrc-bsu.cam.ac.uk>         */
/* R port by Giovanni Petris <GPetris@uark.edu>                           */

/* *********************************************************************** */

#include           <R.h>
#include           <Rmath.h>  /* for RNG interface */
#include           <Rdefines.h>
#include           <stdio.h>
#include           <math.h>
#include           <stdlib.h>

/* *********************************************************************** */

typedef struct point {    /* a point in the x,y plane */
  double x,y;             /* x and y coordinates */
  double ey;              /* exp(y-ymax+YCEIL) */
  double cum;             /* integral up to x of rejection envelope */
  int f;                  /* is y an evaluated point of log-density */
  struct point *pl,*pr;   /* envelope points to left and right of x */
} POINT;

/* *********************************************************************** */

typedef struct envelope {  /* attributes of the entire rejection envelope */
  int cpoint;              /* number of POINTs in current envelope */
  int npoint;              /* max number of POINTs allowed in envelope */
  int *neval;              /* number of function evaluations performed */
  double ymax;             /* the maximum y-value in the current envelope */
  POINT *p;                /* start of storage of envelope POINTs */
  double *convex;          /* adjustment for convexity */
} ENVELOPE;

/* *********************************************************************** */

typedef struct metropolis { /* for metropolis step */
  int on;            /* whether metropolis is to be used */
  double xprev;      /* previous Markov chain iterate */
  double yprev;      /* current log density at xprev */
} METROPOLIS;

/* *********************************************************************** */

#define XEPS  0.00001            /* critical relative x-value difference */
#define YEPS  0.1                /* critical y-value difference */
#define EYEPS 0.001              /* critical relative exp(y) difference */
#define YCEIL 50.                /* maximum y avoiding overflow in exp(y) */

/* *********************************************************************** */

/* declarations for functions defined in this file */

SEXP arms(SEXP bounds, SEXP myldens, SEXP yprev, SEXP size, SEXP rho);

int initial (double *xinit, int ninit, double xl, double xr, int npoint,
	     SEXP myldens, ENVELOPE *env, double *convex, int *neval,
             METROPOLIS *metrop, SEXP rho);

void sample(ENVELOPE *env, POINT *p);

void invert(double prob, ENVELOPE *env, POINT *p);

int test(ENVELOPE *env, POINT *p, SEXP myldens, METROPOLIS *metrop, SEXP rho);

int update(ENVELOPE *env, POINT *p, SEXP myldens, METROPOLIS *metrop, SEXP rho);


void cumulate(ENVELOPE *env);

int meet (POINT *q, ENVELOPE *env, METROPOLIS *metrop);

double area(POINT *q);

double expshift(double y, double y0);

double logshift(double y, double y0);

double perfunc(SEXP myldens, ENVELOPE *env, double x, SEXP rho);

void display(FILE *f, ENVELOPE *env);

double u_random();

/* *********************************************************************** */
# define NINIT 4
SEXP arms (SEXP bounds, SEXP myldens, SEXP yprev, SEXP size, SEXP rho) {

/* to perform derivative-free adaptive rejection sampling with metropolis step */

/* bounds      : boundaries of the support of the density */
/* myldens     : R function to evaluate log density */
/* yprev       : previous value from markov chain */
/* size        : number of sampled values to be obtained */
/* rho         : R environment in which the logdensity is evaluated */

    double xl, xr, xinit[NINIT], convex=1.0;
    int i, npoint=100, nsamp, neval, err;
    SEXP ysamp;         /* sampled values */
    ENVELOPE *env;      /* rejection envelope */
    POINT pwork;        /* a working point, not yet incorporated in envelope */
    int msamp=0;        /* the number of x-values currently sampled */
    METROPOLIS *metrop; /* to hold bits for metropolis step */

    nsamp = INTEGER(size)[0];
    xl = REAL(bounds)[0];
    xr = REAL(bounds)[1];
    for (i=0; i<NINIT; i++)
        xinit[i] = xl + (i + 1.0) * (xr - xl)/(NINIT + 1.0);
    PROTECT( ysamp = NEW_NUMERIC(nsamp) );
        

    /* set up space required for envelope */
    /* env = (ENVELOPE *)malloc(sizeof(ENVELOPE)); */
    env = (ENVELOPE *)Calloc(1, ENVELOPE);
    if(env == NULL){
        /* insufficient space */
        error("insufficient space");
    }

    /* start setting up metropolis struct */
    /* metrop = (METROPOLIS *)malloc(sizeof(METROPOLIS)); */
    metrop = (METROPOLIS *)Calloc(1, METROPOLIS);
    if(metrop == NULL){
        /* insufficient space */
         error("insufficient space");
    }
    metrop->on = 1; 

    /* set up initial envelope */
    err = initial(xinit,NINIT,xl,xr,npoint,myldens,env,&convex,
                  &neval,metrop,rho);
    if(err) 
        error("Can set err...");

    /* finish setting up metropolis struct (can only do this after */
    /* setting up env) */
    if(metrop->on){
        if((REAL(yprev)[0] < xl) || (REAL(yprev)[0] > xr)){
            /* previous markov chain iterate out of range */
            error("previous markov chain iterate out of range");
        }
        metrop->xprev = REAL(yprev)[0];
        metrop->yprev = perfunc(myldens,env,REAL(yprev)[0],rho);
    }

    /* now do adaptive rejection */
    do {
        /* sample a new point */
        sample (env,&pwork);

        /* perform rejection (and perhaps metropolis) tests */
        i = test(env,&pwork,myldens,metrop,rho);
        if(i == 1){
            /* point accepted */
            REAL(ysamp)[msamp++] = pwork.x;
        } else if (i != 0) {
            /* envelope error - violation without metropolis */
            error("envelope error - violation without metropolis");
        }  
    } while (msamp < nsamp);

    /* nsamp points now sampled */

    /* free space */
    Free(env->p);
    Free(env);
    Free(metrop);

    UNPROTECT(1);
    return ysamp;
}

/* *********************************************************************** */

int initial (double *xinit, int ninit, double xl, double xr, int npoint,
	     SEXP myldens, ENVELOPE *env, double *convex, int *neval,
             METROPOLIS *metrop, SEXP rho)

/* to set up initial envelope */

/* xinit        : initial x-values */
/* ninit        : number of initial x-values */
/* xl,xr        : lower and upper x-bounds */
/* npoint       : maximum number of POINTs allowed in envelope */
/* myldens      : R function to evaluate log density */
/* *env         : rejection envelope attributes */
/* *convex      : adjustment for convexity */
/* *neval       : current number of function evaluations */
/* *metrop      : for metropolis step */
/* rho          : R environment in which the logdensity is evaluated */

{
  int i,j,k,mpoint;
  POINT *q;

  if(ninit<3){
    /* too few initial points */
    return 1001;
  }

  mpoint = 2*ninit + 1;
  if(npoint < mpoint){
    /* too many initial points */
    return 1002;
  }

  if((xinit[0] <= xl) || (xinit[ninit-1] >= xr)){
    /* initial points do not satisfy bounds */
    return 1003;
  }

  for(i=1; i<ninit; i++){
    if(xinit[i] <= xinit[i-1]){
      /* data not ordered */
      return 1004;
    }
  }

  if(*convex < 0.0){
    /* negative convexity parameter */
    return 1008;
  }

  /* copy convexity address to env */
  env->convex = convex;

  /* copy address for current number of function evaluations */
  env->neval = neval;
  /* initialise current number of function evaluations */
  *(env->neval) = 0;

  /* set up space for envelope POINTs */
  env->npoint = npoint;
  /* env->p = (POINT *)malloc(npoint*sizeof(POINT)); */
  env->p = (POINT *)Calloc(npoint, POINT);
  if(env->p == NULL){
    /* insufficient space */
    return 1006;
  }

  /* set up envelope POINTs */
  q = env->p;
  /* left bound */
  q->x = xl;
  q->f = 0;
  q->pl = NULL;
  q->pr = q+1;
  for(j=1, k=0; j<mpoint-1; j++){
    q++;
    if(j%2){
      /* point on log density */
      q->x = xinit[k++];
      q->y = perfunc(myldens,env,q->x,rho);
      q->f = 1;
    } else {
      /* intersection point */
      q->f = 0;
    }
    q->pl = q-1;
    q->pr = q+1;
  }
  /* right bound */
  q++;
  q->x = xr;
  q->f = 0;
  q->pl = q-1;
  q->pr = NULL;

  /* calculate intersection points */
  q = env->p;
  for (j=0; j<mpoint; j=j+2, q=q+2){
    if(meet(q,env,metrop)){
      /* envelope violation without metropolis */
      return 2000;
    }
  }

  /* exponentiate and integrate envelope */
  cumulate(env);

  /* note number of POINTs currently in envelope */
  env->cpoint = mpoint;

  return 0;
}

/* *********************************************************************** */

void sample(ENVELOPE *env, POINT *p)

/* To sample from piecewise exponential envelope */
/* *env    : envelope attributes */
/* *p      : a working POINT to hold the sampled value */

{
  double prob;

  /* sample a uniform */
  prob = u_random();
  /* get x-value correponding to a cumulative probability prob */
  invert(prob,env,p);

  return;
}

/* *********************************************************************** */

void invert(double prob, ENVELOPE *env, POINT *p)

/* to obtain a point corresponding to a qiven cumulative probability */
/* prob    : cumulative probability under envelope */
/* *env    : envelope attributes */
/* *p      : a working POINT to hold the sampled value */

{
  double u,xl,xr,yl,yr,eyl,eyr,prop,z;
  POINT *q;

  /* find rightmost point in envelope */
  q = env->p;
  while(q->pr != NULL)q = q->pr;

  /* find exponential piece containing point implied by prob */
  u = prob * q->cum;
  while(q->pl->cum > u)q = q->pl;

  /* piece found: set left and right POINTs of p, etc. */
  p->pl = q->pl;
  p->pr = q;
  p->f = 0;
  p->cum = u;

  /* calculate proportion of way through integral within this piece */
  prop = (u - q->pl->cum) / (q->cum - q->pl->cum);

  /* get the required x-value */
  if (q->pl->x == q->x){
    /* interval is of zero length */
    p->x = q->x;
    p->y = q->y;
    p->ey = q->ey;
  } else {
    xl = q->pl->x;
    xr = q->x;
    yl = q->pl->y;
    yr = q->y;
    eyl = q->pl->ey;
    eyr = q->ey;
    if(fabs(yr - yl) < YEPS){
      /* linear approximation was used in integration in function cumulate */
      if(fabs(eyr - eyl) > EYEPS*fabs(eyr + eyl)){
	p->x = xl + ((xr - xl)/(eyr - eyl))
	       * (-eyl + sqrt((1. - prop)*eyl*eyl + prop*eyr*eyr));
      } else {
	p->x = xl + (xr - xl)*prop;
      }
      p->ey = ((p->x - xl)/(xr - xl)) * (eyr - eyl) + eyl;
      p->y = logshift(p->ey, env->ymax);
    } else {
      /* piece was integrated exactly in function cumulate */
      p->x = xl + ((xr - xl)/(yr - yl))
	      * (-yl + logshift(((1.-prop)*eyl + prop*eyr), env->ymax));
      p->y = ((p->x - xl)/(xr - xl)) * (yr - yl) + yl;
      p->ey = expshift(p->y, env->ymax);
    }
  }

  /* guard against imprecision yielding point outside interval */
  if ((p->x < xl) || (p->x > xr)) error("exit 1");

  return;
}

/* *********************************************************************** */

int test(ENVELOPE *env, POINT *p, SEXP myldens, METROPOLIS *metrop, SEXP rho)

/* to perform rejection, squeezing, and metropolis tests */

/* *env          : envelope attributes */
/* *p            : point to be tested */
/* myldens       : R function to evaluate log-density */
/* *metrop       : data required for metropolis step */
/* rho           : R environment in which the logdensity is evaluated */

{
  double u,y,ysqueez,ynew,yold,znew,zold,w;
  POINT *ql,*qr;
  
  /* for rejection test */
  u = u_random() * p->ey;
  y = logshift(u,env->ymax);

  if(!(metrop->on) && (p->pl->pl != NULL) && (p->pr->pr != NULL)){
    /* perform squeezing test */
    if(p->pl->f){
      ql = p->pl;
    } else {
      ql = p->pl->pl;
    }
    if(p->pr->f){
      qr = p->pr;
    } else {
      qr = p->pr->pr;
    }
    ysqueez = (qr->y * (p->x - ql->x) + ql->y * (qr->x - p->x))
               /(qr->x - ql->x);
    if(y <= ysqueez){
      /* accept point at squeezing step */
      return 1;
    }
  }

  /* evaluate log density at point to be tested */
  ynew = perfunc(myldens,env,p->x,rho);
  
  /* perform rejection test */
  if(!(metrop->on) || ((metrop->on) && (y >= ynew))){
    /* update envelope */
    p->y = ynew;
    p->ey = expshift(p->y,env->ymax);
    p->f = 1;
    if(update(env,p,myldens,metrop,rho)){
      /* envelope violation without metropolis */
      return -1;
    }
    /* perform rejection test */
    if(y >= ynew){
      /* reject point at rejection step */
      return 0;
    } else {
      /* accept point at rejection step */
      return 1;
    }
  }

  /* continue with metropolis step */
  yold = metrop->yprev;
  /* find envelope piece containing metrop->xprev */
  ql = env->p;
  while(ql->pl != NULL)ql = ql->pl;
  while(ql->pr->x < metrop->xprev)ql = ql->pr;
  qr = ql->pr;
  /* calculate height of envelope at metrop->xprev */
  w = (metrop->xprev - ql->x)/(qr->x - ql->x);
  zold = ql->y + w*(qr->y - ql->y);
  znew = p->y;
  if(yold < zold)zold = yold;
  if(ynew < znew)znew = ynew;
  w = ynew-znew-yold+zold;
  if(w > 0.0)w = 0.0;

  if(w > -YCEIL){
    w = exp(w);
  } else {
    w = 0.0;
  }
  u = u_random();
  if(u > w){
    /* metropolis says dont move, so replace current point with previous */
    /* markov chain iterate */
    p->x = metrop->xprev;
    p->y = metrop->yprev;
    p->ey = expshift(p->y,env->ymax);
    p->f = 1;
    p->pl = ql;
    p->pr = qr;
  } else {
    /* trial point accepted by metropolis, so update previous markov */
    /* chain iterate */
    metrop->xprev = p->x;
    metrop->yprev = ynew;
  }
  return 1;
}

/* *********************************************************************** */

int update(ENVELOPE *env, POINT *p, SEXP myldens, METROPOLIS *metrop, SEXP rho)

/* to update envelope to incorporate new point on log density*/
/* *env          : envelope attributes */
/* *p            : point to be incorporated */
/* myldens       : R function to evaluate log density */
/* *metrop       : for metropolis step */
/* rho           : R environment in which the logdensity is evaluated */

{
  POINT *m,*ql,*qr,*q;

  if(!(p->f) || (env->cpoint > env->npoint - 2)){
    /* y-value has not been evaluated or no room for further points */
    /* ignore this point */
    return 0;
  }

  /* copy working POINT p to a new POINT q */
  q = env->p + env->cpoint++;
  q->x = p->x;
  q->y = p->y;
  q->f = 1;

  /* allocate an unused POINT for a new intersection */
  m = env->p + env->cpoint++;
  m->f = 0;
  if((p->pl->f) && !(p->pr->f)){
    /* left end of piece is on log density; right end is not */
    /* set up new intersection in interval between p->pl and p */
    m->pl = p->pl;
    m->pr = q;
    q->pl = m;
    q->pr = p->pr;
    m->pl->pr = m;
    q->pr->pl = q;
  } else if (!(p->pl->f) && (p->pr->f)){
    /* left end of interval is not on log density; right end is */
    /* set up new intersection in interval between p and p->pr */
    m->pr = p->pr;
    m->pl = q;
    q->pr = m;
    q->pl = p->pl;
    m->pr->pl = m;
    q->pl->pr = q;
  } else {
    /* this should be impossible */
      error("exit 10");
  }

  /* now adjust position of q within interval if too close to an endpoint */
  if(q->pl->pl != NULL){
    ql = q->pl->pl;
  } else {
    ql = q->pl;
  }
  if(q->pr->pr != NULL){
    qr = q->pr->pr;
  } else {
    qr = q->pr;
  }
  if (q->x < (1. - XEPS) * ql->x + XEPS * qr->x){
    /* q too close to left end of interval */
    q->x = (1. - XEPS) * ql->x + XEPS * qr->x;
    q->y = perfunc(myldens,env,q->x,rho);
  } else if (q->x > XEPS * ql->x + (1. - XEPS) * qr->x){
    /* q too close to right end of interval */
    q->x = XEPS * ql->x + (1. - XEPS) * qr->x;
    q->y = perfunc(myldens,env,q->x,rho);
  }

  /* revise intersection points */
  if(meet(q->pl,env,metrop)){
    /* envelope violation without metropolis */
    return 1;
  }
  if(meet(q->pr,env,metrop)){
    /* envelope violation without metropolis */
    return 1;
  }
  if(q->pl->pl != NULL){
    if(meet(q->pl->pl->pl,env,metrop)){
      /* envelope violation without metropolis */
      return 1;
    }
  }
  if(q->pr->pr != NULL){
    if(meet(q->pr->pr->pr,env,metrop)){
      /* envelope violation without metropolis */
      return 1;
    }
  }

  /* exponentiate and integrate new envelope */
  cumulate(env);

  return 0;
}

/* *********************************************************************** */

void cumulate(ENVELOPE *env)

/* to exponentiate and integrate envelope */
/* *env     : envelope attributes */

{
  POINT *q,*qlmost;

  qlmost = env->p;
  /* find left end of envelope */
  while(qlmost->pl != NULL)qlmost = qlmost->pl;

  /* find maximum y-value: search envelope */
  env->ymax = qlmost->y;
  for(q = qlmost->pr; q != NULL; q = q->pr){
    if(q->y > env->ymax)env->ymax = q->y;
  }

  /* exponentiate envelope */
  for(q = qlmost; q != NULL; q = q->pr){
    q->ey = expshift(q->y,env->ymax);
  }

  /* integrate exponentiated envelope */
  qlmost->cum = 0.;
  for(q = qlmost->pr; q != NULL; q = q->pr){
    q->cum = q->pl->cum + area(q);
  }

  return;
}

/* *********************************************************************** */

int meet (POINT *q, ENVELOPE *env, METROPOLIS *metrop)
/* To find where two chords intersect */
/* q         : to store point of intersection */
/* *env      : envelope attributes */
/* *metrop   : for metropolis step */

{
  double gl,gr,grl,dl,dr;
  int il,ir,irl;

  if(q->f){
    /* this is not an intersection point */
      error("exit 30");
  }

  /* calculate coordinates of point of intersection */
  if ((q->pl != NULL) && (q->pl->pl->pl != NULL)){
    /* chord gradient can be calculated at left end of interval */
    gl = (q->pl->y - q->pl->pl->pl->y)/(q->pl->x - q->pl->pl->pl->x);
    il = 1;
  } else {
    /* no chord gradient on left */
    il = 0;
  }
  if ((q->pr != NULL) && (q->pr->pr->pr != NULL)){
    /* chord gradient can be calculated at right end of interval */
    gr = (q->pr->y - q->pr->pr->pr->y)/(q->pr->x - q->pr->pr->pr->x);
    ir = 1;
  } else {
    /* no chord gradient on right */
    ir = 0;
  }
  if ((q->pl != NULL) && (q->pr != NULL)){
    /* chord gradient can be calculated across interval */
    grl = (q->pr->y - q->pl->y)/(q->pr->x - q->pl->x);
    irl = 1;
  } else {
    irl = 0;
  }

  if(irl && il && (gl<grl)){
    /* convexity on left exceeds current threshold */
    if(!(metrop->on)){
      /* envelope violation without metropolis */
      return 1;
    }
    /* adjust left gradient */
    gl = gl + (1.0 + *(env->convex)) * (grl - gl);
  }

  if(irl && ir && (gr>grl)){
    /* convexity on right exceeds current threshold */
    if(!(metrop->on)){
      /* envelope violation without metropolis */
      return 1;
    }
    /* adjust right gradient */
    gr = gr + (1.0 + *(env->convex)) * (grl - gr);
  }

  if(il && irl){
    dr = (gl - grl) * (q->pr->x - q->pl->x);
    if(dr < YEPS){
      /* adjust dr to avoid numerical problems */
      dr = YEPS;
    }
  }

  if(ir && irl){
    dl = (grl - gr) * (q->pr->x - q->pl->x);
    if(dl < YEPS){
      /* adjust dl to avoid numerical problems */
      dl = YEPS;
    }
  }

  if(il && ir && irl){
    /* gradients on both sides */
    q->x = (dl * q->pr->x + dr * q->pl->x)/(dl + dr);
    q->y = (dl * q->pr->y + dr * q->pl->y + dl * dr)/(dl + dr);
  } else if (il && irl){
    /* gradient only on left side, but not right hand bound */
    q->x = q->pr->x;
    q->y = q->pr->y + dr;
  } else if (ir && irl){
    /* gradient only on right side, but not left hand bound */
    q->x = q->pl->x;
    q->y = q->pl->y + dl;
  } else if (il){
    /* right hand bound */
    q->y = q->pl->y + gl * (q->x - q->pl->x);
  } else if (ir){
    /* left hand bound */
    q->y = q->pr->y - gr * (q->pr->x - q->x);
  } else {
    /* gradient on neither side - should be impossible */
      error("exit 31");
  }
  if(((q->pl != NULL) && (q->x < q->pl->x)) ||
     ((q->pr != NULL) && (q->x > q->pr->x))){
    /* intersection point outside interval (through imprecision) */
      error("exit 32");
  }
  /* successful exit : intersection has been calculated */
  return 0;
}

/* *********************************************************************** */

double area(POINT *q)

/* To integrate piece of exponentiated envelope to left of POINT q */

{
  double a;

  if(q->pl == NULL){
    /* this is leftmost point in envelope */
      error("exit 1");
  } else if(q->pl->x == q->x){
    /* interval is zero length */
    a = 0.;
  } else if (fabs(q->y - q->pl->y) < YEPS){
    /* integrate straight line piece */
    a = 0.5*(q->ey + q->pl->ey)*(q->x - q->pl->x);
  } else {
    /* integrate exponential piece */
    a = ((q->ey - q->pl->ey)/(q->y - q->pl->y))*(q->x - q->pl->x);
  }
  return a;
}

/* *********************************************************************** */

double expshift(double y, double y0)

/* to exponentiate shifted y without underflow */
{
  if(y - y0 > -2.0 * YCEIL){
    return exp(y - y0 + YCEIL);
  } else {
    return 0.0;
  }
}

/* *********************************************************************** */

double logshift(double y, double y0)

/* inverse of function expshift */
{
  return (log(y) + y0 - YCEIL);
}

/* *********************************************************************** */

double perfunc(SEXP myldens, ENVELOPE *env, double x, SEXP rho)

/* to evaluate log density and increment count of evaluations */

/* myldens : R function to evaluate log density */
/* *env    : envelope attributes */
/* x       : point at which to evaluate log density */
/* rho     : R environment in which the logdensity is evaluated */
{
  double y;
  SEXP R_fcall, arg;

  /* evaluate logdensity function */
  PROTECT(R_fcall = lang2(myldens, R_NilValue));
  PROTECT(arg = NEW_NUMERIC(1));
  NUMERIC_POINTER(arg)[0] = x;
  SETCADR(R_fcall, arg);
  y = REAL(eval(R_fcall, rho))[0];
  UNPROTECT(2);

  /* increment count of function evaluations */
  (*(env->neval))++;

  return y;
}

/* *********************************************************************** */

void display(FILE *f, ENVELOPE *env)

/* to display envelope - for debugging only */
{
  POINT *q;

  /* print envelope attributes */
  fprintf(f,"========================================================\n");
  fprintf(f,"envelope attributes:\n");
  fprintf(f,"points in use = %d, points available = %d\n",
          env->cpoint,env->npoint);
  fprintf(f,"function evaluations = %d\n",*(env->neval));
  fprintf(f,"ymax = %f, p = %x\n",env->ymax,env->p);
  fprintf(f,"convexity adjustment = %f\n",*(env->convex));
  fprintf(f,"--------------------------------------------------------\n");

  /* find leftmost POINT */
  q = env->p;
  while(q->pl != NULL)q = q->pl;

  /* now print each POINT from left to right */
  for(q = env->p; q != NULL; q = q->pr){
    fprintf(f,"point at %x, left at %x, right at %x\n",q,q->pl,q->pr);
    fprintf(f,"x = %f, y = %f, ey = %f, cum = %f, f = %d\n",
            q->x,q->y,q->ey,q->cum,q->f);
  }
  fprintf(f,"========================================================\n");

  return;
}

/* *********************************************************************** */

double u_random()

/* to return a standard uniform random number */
{
    double x;
    GetRNGstate();
    x = unif_rand();
    PutRNGstate();
    return x;

}

/* *********************************************************************** */

