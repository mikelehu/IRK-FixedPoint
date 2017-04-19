void mathIRKFP P((int neq, int ns, double t0, double tend, double *u, long ulen, double *ee, long elen,  double h,double *rpar, long prelem, int *ipar, long pielem,int approximation,int threads,int iteration,int rdigits,const char *filename,int sampling,int codfun)); 

:Begin:
:Function:       mathIRKFP
:Pattern:  mathIRKFP [neq_Integer,ns_Integer,t0_Real,tend_Real,u_List,ee_List,h_Real,rpar_List,ipar_List,approximation_Integer, 
 threads_Integer,iteration_Integer,rdigits_Integer,filename_String,sampling_Integer, codfun_Integer]
:Arguments:      {neq,ns,t0,tend,u,ee,h,rpar,ipar,approximation,threads,iteration,rdigits,filename, sampling,codfun}
:ArgumentTypes:  {Integer,Integer, Real, Real, RealList,RealList,Real,RealList, IntegerList, Integer, Integer,Integer, Integer, String, Integer,Integer}
:ReturnType:     Manual
:End:

:Evaluate: mathGauss::usage = "..."


