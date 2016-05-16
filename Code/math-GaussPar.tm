void mathGaussPar P((int neq, int n, int ns, double t0, double tend, double *u, long ulen, double *ee, long elen,  double h,double *rpar, long prelem, int *ipar, long pielem,int approximation,int threads,int algoritmoa1,int algoritmoa2,int rdigits1,int rdigits2,const char *filename1,const char *filename2,int sampling,int codfun)); 

:Begin:
:Function:       mathGaussPar
:Pattern:  mathGaussPar [neq_Integer,n_Integer,ns_Integer,t0_Real,tend_Real,u_List,ee_List,h_Real,rpar_List,ipar_List,approximation_Integer, threads_Integer,algoritmoa1_Integer,algoritmoa2_Integer,rdigits1_Integer,rdigits2_Integer,filename1_String,filename2_String,sampling_Integer, codfun_Integer]
:Arguments:{neq,n,ns,t0,tend,u,ee,h,rpar,ipar,approximation,threads,algoritmoa1,algoritmoa2,rdigits1,rdigits2,filename1,filename2,sampling,codfun}
:ArgumentTypes:  {Integer,Integer, Integer, Real, Real, RealList,RealList,Real,RealList, IntegerList, Integer, Integer,Integer, Integer,Integer,Integer,String,String,Integer,Integer}
:ReturnType:     Manual
:End:

:Evaluate: mathGaussPar::usage = "azalpenak..."




