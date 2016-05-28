/*------------------------------------------------------------------------------*/
/*										*/
/*                                GaussCommon.h					*/
/*										*/
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
(int neq,val_type *u
);


void InitStat
(ode_sys *system,gauss_method *gsmethod, solver_stat *thestatptr
 );

val_type NormalizedDistance
(int neq,int ns,toptions *options,val_type *z,val_type *zold
 );

int Yi_init
(solution *u, val_type *z,ode_sys *system,
 gauss_method *method,solver_stat *thestatptr, toptions *options
);

void RemoveDigitsFcn
(ode_sys *system,gauss_method *gsmethod,val_type *z, int m
);

void UpdateDMin
(ode_sys *system,gauss_method *method,
 bool *D0,bool *cont,val_type *DMin,val_type *Y, val_type *Yold
);


void Fixed_point_it
(ode_sys *system, solution *u, val_type tn,val_type h, 
 toptions *options,gauss_method *method,solver_stat *thestatptr
);

int It_Jacobi
(ode_sys *system, solution *u, val_type tn,val_type h, 
 gauss_method *method,solver_stat *thestatptr
);

int It_Jacobi_Classic
(ode_sys *system, solution *u, val_type tn,val_type h, 
 gauss_method *method,solver_stat *thestatptr
);

int It_Seidel
(ode_sys *system, solution *u, val_type tn,val_type h, 
 gauss_method *method,solver_stat *thestatptr
);

int It_Seidel_Classic
(ode_sys *system, solution *u, val_type tn,val_type h, 
 gauss_method *method,solver_stat *thestatptr
);


void TheOutput
(ode_sys *system,val_type t,solution *u,solver_stat *thestatptr,
 parameters *params,toptions *options,FILE *loga
);

void TheOutput2
(ode_sys *system,val_type t,solution *u,solution *u2,solver_stat *thestatptr,
 parameters *params,toptions *options,FILE *loga
);

void RKG 
(gauss_method *gsmethod,
 solution *u,
 ode_sys *system, toptions *options,
 void RKG_Step (), solver_stat *thestatptr
);

void RKG2 
(gauss_method *gsmethod, gauss_method *gsmethod2,
 solution *u, solution *u2,
 ode_sys *system,toptions *options,toptions *options2,
 void RKG_Step (), solver_stat *thestatptr, solver_stat *thestatptr2
);

void select_gauss
(gauss_method *gsmethodptr, gauss_method *gsmethod2ptr,
 solution *uptr,solution *u2ptr,ode_sys *systemptr,
 toptions *optionsptr,toptions *options2ptr,
 solver_stat *thestatptr,solver_stat *thestat2ptr
);




