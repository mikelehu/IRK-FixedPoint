/*-----------------------------------------------------------------------------*/
/*									       */
/*                                GaussCommon.h				       */
/*									       */
/* ----------------------------------------------------------------------------*/

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include <def.h>
#include <sys/stat.h>
#include <Problems.h> 
#include <GaussCoefficients.h>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <mathlink.h>


void print_u 
( int neq,  val_type *u
);

void InitStat
( ode_sys *system, gauss_method *gsmethod, solver_stat *thestatptr
);

val_type NormalizedDistance
( int neq, int ns, toptions *options,
  val_type *z, val_type *zold
);

int StatYinit
( ode_sys *system, gauss_method *method,
 solver_stat *thestatptr
);

void RemoveDigits
(val_type *x, int m
);

void Default_Stage_init
( solution *u, val_type *z, ode_sys *system,
  gauss_method *method,solver_stat *thestatptr,  toptions *options
);

void Interpolated_Stage_init
( solution *u, val_type *z, ode_sys *system,
  gauss_method *method,solver_stat *thestatptr,  toptions *options
);

void StopCriterion
( ode_sys *system,  gauss_method *method,
 int *D0,bool *cont,val_type *DMin,
  val_type *Y,  val_type *Yold
);


void Fixed_point_Step
( ode_sys *system,   solution *u, 
  val_type tn,  val_type h, 
  toptions *options,  gauss_method *method,
 solver_stat *thestatptr
);


int General_FP_It
( ode_sys *system,  solution *u,  val_type tn,
  val_type h,  gauss_method *method,solver_stat *thestatptr
);

int Partitioned_FP_It
( ode_sys *system,  solution *u,  val_type tn,
  val_type h,  gauss_method *method,solver_stat *thestatptr
);

void MyOutput
( ode_sys *system,  gauss_method *method,
  val_type t, val_type h, solution *u,
 solver_stat *thestatptr,
  parameters *params, toptions *options,FILE *loga
);


void CompensatedSummation 
( gauss_method *gsmethod,
 solution *u,
  ode_sys *system,  toptions *options,
  solver_stat *thestatptr
);


void IRKFP 
(val_type t0, val_type t1, val_type h,
  gauss_method *gsmethod, solution *u,
  ode_sys *system, toptions *options,
 solver_stat *thestatptr
);




