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
#include <GaussUserProblem.h> 
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

int statYinit
(const ode_sys *system,const gauss_method *method,
 solver_stat *thestatptr
);

void RemoveDigitsFcn
(val_type *x,const int m
);

int Yi_init
(const solution *u, val_type *z,const ode_sys *system,
 const gauss_method *method,solver_stat *thestatptr,const  toptions *options
);


void StopCriterion
(const ode_sys *system, const gauss_method *method,
 int *D0,bool *cont,val_type *DMin,
 const val_type *Y, const val_type *Yold
);


void Fixed_point_it
(const ode_sys *system, const  solution *u, 
 const val_type tn, const val_type h, 
 const toptions *options, const gauss_method *method,
 solver_stat *thestatptr
);

int It_Jacobi
(const ode_sys *system, const solution *u, const val_type tn,
 const val_type h, const gauss_method *method,solver_stat *thestatptr
);


void TheOutput
(const ode_sys *system, const gauss_method *method,
 const val_type t, const solution *u,
 solver_stat *thestatptr,
 const parameters *params,const toptions *options,FILE *loga
);

void TheOutput2
(const ode_sys *system, const gauss_method *method,
 const val_type t, const solution *u,
 const solution *u2,solver_stat *thestatptr,
 const parameters *params,const toptions *options,FILE *loga
);


void CompensatedSummation 
(const gauss_method *gsmethod,
 val_type *u0,solution *u,
 const ode_sys *system, const toptions *options,
 const solver_stat *thestatptr
);


void RKG 
(const gauss_method *gsmethod,
 solution *u,
 const ode_sys *system, toptions *options,
 void RKG_Step (), solver_stat *thestatptr
);


void RKG2 
(const gauss_method *gsmethod, const gauss_method *gsmethod2,
 solution *u, solution *u2,
 const ode_sys *system, toptions *options,
 toptions *options2,
 void RKG_Step (), solver_stat *thestatptr, solver_stat *thestatptr2
);


void select_gauss
(gauss_method *gsmethodptr, gauss_method *gsmethod2ptr,
 solution *uptr,solution *u2ptr,
 ode_sys *systemptr,
 toptions *optionsptr,toptions *options2ptr,
 solver_stat *thestatptr,solver_stat *thestat2ptr
);

void select_odefun
(const int codfun, ode_sys *system
);



