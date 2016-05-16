/*------------------------------------------------------------------------------*/
/*										*/
/*         GaussCoefficients.c							*/
/*		Implicit Runge-Kutta (based Gauss quadrature)  coeficients.     */
/*              			                                        */
/*              	Method coefficients: m,b,c.                             */
/*                      Interpolation coefficientes: nu.                        */     
/*              			                                        */
/*                                                                              */      
/* -----------------------------------------------------------------------------*/

#include <GaussCoefficients.h>


void GaussCoefficients (gauss_method *method,toptions *options)
{
     int i,j;
     int ns=method->ns;
     val_type sum;

     char mydir[20];
     FILE *fileM,*fileA,*fileB,*fileC,*fileNU;

     char filenameM[STRMAX];
     char filenameA[STRMAX];
     char filenameB[STRMAX];
     char filenameC[STRMAX];
     char filenameNU[STRMAX];

     method->m =(val_type *) malloc(ns*ns*sizeof(val_type));
     method->b = (val_type *) malloc(ns*sizeof(val_type));     
     method->hb = (val_type *) malloc(ns*sizeof(val_type));      
     method->c = (val_type *) malloc((ns)*sizeof(val_type));
     method->hc = (val_type *) malloc((ns)*sizeof(val_type));
     method->nu = (val_type *) malloc(ns*ns*sizeof(val_type));
     method->orderedindices=(int *) malloc(ns*sizeof(int));

#    if PREC ==2  //QUADRUPLEPRECISION
     int n;
     int width = 46;
     char buf[128];
     val_type aux;
#    endif


     switch (ns)
     { 
     case 6: 

           // orderindices (ascending order).
           method->orderedindices[0]=0;   
	   method->orderedindices[1]=5;
       	   method->orderedindices[2]=1;
	   method->orderedindices[3]=4;
	   method->orderedindices[4]=2;
	   method->orderedindices[5]=3;

//	   strcpy(mydir,"./CoefficientsData/S6/");
           strcpy(mydir,"../../CoefficientsData/S6/");

     break;

     case 8: 

           // orderindices (ascending order).
           method->orderedindices[0]=0;   
	   method->orderedindices[1]=7;
           method->orderedindices[2]=1;
	   method->orderedindices[3]=6;
	   method->orderedindices[4]=2;
	   method->orderedindices[5]=5;
	   method->orderedindices[6]=3;
	   method->orderedindices[7]=4;

//           strcpy(mydir,"./CoefficientsData/S8/");
           strcpy(mydir,"../../CoefficientsData/S8/");

     break;

     case 16: 

           // orderindices (ascending order).
           method->orderedindices[0]=0;   
	   method->orderedindices[1]=15;
       	   method->orderedindices[2]=1;
	   method->orderedindices[3]=14;
	   method->orderedindices[4]=2;
	   method->orderedindices[5]=13;
	   method->orderedindices[6]=3;
	   method->orderedindices[7]=12;
	   method->orderedindices[8]=4;
	   method->orderedindices[9]=11;
	   method->orderedindices[10]=5;
	   method->orderedindices[11]=10;
	   method->orderedindices[12]=6;
           method->orderedindices[13]=9;
           method->orderedindices[14]=7;
           method->orderedindices[15]=8;

//           strcpy(mydir,"./CoefficientsData/S16/");
           strcpy(mydir,"../../CoefficientsData/S16/");

     break;


     default:
          printf("Coefficients not defined\n");
     break;

     }

     strcpy(filenameM,mydir);
     strcpy(filenameA,mydir);
     strcpy(filenameB,mydir);
     strcpy(filenameC,mydir);
     strcpy(filenameNU,mydir);

#if PREC ==2  //QUADRUPLEPRECISION

     strcat(filenameM,"QMCoef.bin");
     strcat(filenameA,"QACoef.bin");
     strcat(filenameB,"QBCoef.bin");
     strcat(filenameC,"QCCoef.bin");
     strcat(filenameNU,"QNUCoef.bin");

#elif PREC ==3 //Float

     strcat(filenameM,"FMCoef.bin");
     strcat(filenameA,"FACoef.bin");
     strcat(filenameB,"FBCoef.bin");
     strcat(filenameC,"FCCoef.bin");
     strcat(filenameNU,"FNUCoef.bin");

#else      // DOUBLEPRECISION

     strcat(filenameM,"DMCoef.bin");
     strcat(filenameA,"DACoef.bin");
     strcat(filenameB,"DBCoef.bin");
     strcat(filenameC,"DCCoef.bin");
     strcat(filenameNU,"DNUCoef.bin");
 
#    endif

 
    fileM = fopen(filenameM,"rb");
    if (fileM == NULL) printf("File doesnt exists\n");

    fileA = fopen(filenameA,"rb");
    fileB = fopen(filenameB,"rb");
    fileC = fopen(filenameC,"rb");
    fileNU = fopen(filenameNU,"rb");

    fread(method->m, sizeof(val_type),ns*ns,fileM);
//    fread(method->a, sizeof(val_type),ns*ns,fileA);
    fread(method->b, sizeof(val_type),ns,fileB);
    fread(method->c, sizeof(val_type),ns,fileC);
    fread(method->nu, sizeof(val_type),ns*ns,fileNU);

/*---- Verify symplectic condition --------------------------------*/


#    if PREC ==2  //QUADRUPLEPRECISION
     for (i=0; i<ns; i++)
     {
        for (j=0; j<ns; j++)
        {
              aux=method->m[i*ns+j]+method->m[j*ns+i]-1.q;
              n = quadmath_snprintf(buf, sizeof buf, "%+-#*.30Qe", width, aux);
              if ((size_t) n < sizeof buf) printf("%s\n",buf);
        }

     }
#    else //DOUBLEPRECISION
     for (i=0; i<ns; i++)
     {
        printf("\n");
        for (j=0; j<ns; j++) printf ("%lg,", method->m[i*ns+j]+method->m[j*ns+i]-1.);
     }

#    endif

     printf("\n");

/*---- Calculate hb coefficients ----------------------------------*/

     sum=0.;
     for (i=1; i<ns-1; i++)
     {
        method->hb[i]=(options->h)*method->b[i];
        sum+=method->hb[i];
     }
     
     method->hb[0]=((options->h)-sum)/2;
     method->hb[ns-1]=((options->h)-sum)/2;


     fclose(fileM);
     fclose(fileA);
     fclose(fileB);
     fclose(fileC);
     fclose(fileNU);

     return;
}




