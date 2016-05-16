void mathGauss P((int neq, int n, int ns, double t0, double tend, double *u, long ulen, double *ee, long elen,  double h,double *rpar, long prelem, int *ipar, long pielem,int approximation,int threads,int algorithm,int rdigits1,int rdigits2,const char *filename,int sampling,int codfun)); 

:Begin:
:Function:       mathGauss
:Pattern:  mathGauss [neq_Integer,n_Integer,ns_Integer,t0_Real,tend_Real,u_List,ee_List,h_Real,rpar_List,ipar_List,approximation_Integer, 
 threads_Integer,algorithm_Integer,rdigits1_Integer,rdigits2_Integer,filename_String,sampling_Integer, codfun_Integer]
:Arguments:      {neq,n,ns,t0,tend,u,ee,h,rpar,ipar,approximation,threads,algorithm,rdigits1,rdigits2,filename, sampling,codfun}
:ArgumentTypes:  {Integer,Integer,Integer, Real, Real, RealList,RealList,Real,RealList, IntegerList, Integer, Integer,Integer, Integer, Integer, String, Integer,Integer}
:ReturnType:     Manual
:End:

:Evaluate: mathGauss::usage = "azalpenak"


