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
(const int neq, const val_type *u
);

void InitStat
(const ode_sys *system,const gauss_method *gsmethod, solver_stat *thestatptr
);

val_type NormalizedDistance
(const int neq,const int ns,const toptions *options,
 const val_type *z,const val_type *zold
);

int StatYinit
(const ode_sys *system,const gauss_method *method,
 solver_stat *thestatptr
);

void RemoveDigits
(val_type *x,const int m
);

void Default_Stage_init
(const solution *u, val_type *z,const ode_sys *system,
 const gauss_method *method,solver_stat *thestatptr,const  toptions *options
);

void Interpolated_Stage_init
(const solution *u, val_type *z,const ode_sys *system,
 const gauss_method *method,solver_stat *thestatptr,const  toptions *options
);

void StopCriterion
(const ode_sys *system, const gauss_method *method,
 int *D0,bool *cont,val_type *DMin,
 const val_type *Y, const val_type *Yold
);


void Fixed_point_Step
(const ode_sys *system, const  solution *u, 
 const val_type tn, const val_type h, 
 const toptions *options, const gauss_method *method,
 solver_stat *thestatptr
);


int General_FP_It
(const ode_sys *system, const solution *u, const val_type tn,
 const val_type h, const gauss_method *method,solver_stat *thestatptr
);

int Partitioned_FP_It
(const ode_sys *system, const solution *u, const val_type tn,
 const val_type h, const gauss_method *method,solver_stat *thestatptr
);

void MyOutput
(const ode_sys *system, const gauss_method *method,
 const val_type t, val_type h,const solution *u,
 solver_stat *thestatptr,
 const parameters *params,const toptions *options,FILE *loga
);


void CompensatedSummation 
(const gauss_method *gsmethod,
 solution *u,
 const ode_sys *system, const toptions *options,
 const solver_stat *thestatptr
);


void IRKFP 
(val_type t0, val_type t1, val_type h,
 const gauss_method *gsmethod, solution *u,
 const ode_sys *system, toptions *options,
 solver_stat *thestatptr
);




