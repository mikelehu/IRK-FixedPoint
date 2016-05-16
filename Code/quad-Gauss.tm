
void quadGauss P((int neq, int n, int ns, double t0, double tend, const unsigned char *quadak, int qluz, const unsigned char *elsak, int elluz, double h, const unsigned char *realpars, int numrealpars, int *ipar, long ilen, int approximation, int threads, int algorithm, int rdigits1,int rdigits2,const char *myfilename, int sampling, int codfun));

:Begin:
:Function:       quadGauss
:Pattern:        quadGauss[neq_Integer, n_Integer, ns_Integer, t0_Real, t1_Real, dat_String, els_String, h_Real, rpars_String, ipars_List, approx_Integer, threads_Integer, alg_Integer, rdigits1_Integer, rdigits2_Integer, filename_String, sampling_Integer, codfun_Integer ]
:Arguments:      {neq,n,ns,t0,t1,dat,els, h,rpars,ipars,approx,threads,alg,rdigits1,rdigits2,filename, sampling,codfun}
:ArgumentTypes:  {Integer,Integer,Integer,Real,Real,ByteString,ByteString, Real, ByteString,IntegerList,Integer,Integer, Integer, Integer, Integer, String, Integer, Integer}
:ReturnType:     Manual
:End:

:Evaluate: quadGauss::usage = "irteera = quadGauss[neq_Integer, n_Integer, ns_Real, t0_Real, t1_Real, dat_String,els_String, h_Real, rpars_String, ipars_List,approx_Integer, threads_Integer, alg_Integer, rdigits_Integer, filename_String, sampling_Integer, codfun_Integer] "
