/*------------------------------------------------------------------------------*/
/*										*/
/*         GaussCoefficients.c							*/
/*		Implicit Runge-Kutta (based Gauss quadrature)  coeficients.     */
/*              			                                        */
/*              	Method coefficients: m,a,b,c.                           */
/*                      Interpolation coefficientes: nu.                        */     
/*              			                                        */
/*                                                                              */      
/* -----------------------------------------------------------------------------*/

#include <GaussCoefficients.h>


void GaussCoefficients
(const char *path, gauss_method *method,  val_type h)

{
     int i,j;
     int info;
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
     method->a =(val_type *) malloc(ns*ns*sizeof(val_type));
     method->b = (val_type *) malloc(ns*sizeof(val_type));     
     method->hb = (val_type *) malloc(ns*sizeof(val_type));      
     method->c = (val_type *) malloc((ns)*sizeof(val_type));
     method->hc = (val_type *) malloc((ns)*sizeof(val_type));
     method->nu = (val_type *) malloc(ns*ns*sizeof(val_type));
     method->orderedindices=(int *) malloc(ns*sizeof(int));

     strcpy(mydir,path);

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

           strcat(mydir,"S6/");

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

           strcat(mydir,"S8/");

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

           strcat(mydir,"S16/");

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


     strcat(filenameM,"DMCoef.bin");
     strcat(filenameA,"DACoef.bin");
     strcat(filenameB,"DBCoef.bin");
     strcat(filenameC,"DCCoef.bin");
     strcat(filenameNU,"DNUCoef.bin");

 
    fileM = fopen(filenameM,"rb");
    if (fileM == NULL) printf("File doesnt exists\n");
    fileA = fopen(filenameA,"rb");
    if (fileA == NULL) printf("File doesnt exists\n");
    fileB = fopen(filenameB,"rb");
    if (fileB == NULL) printf("File doesnt exists\n");
    fileC = fopen(filenameC,"rb");
    if (fileC == NULL) printf("File doesnt exists\n");
    fileNU = fopen(filenameNU,"rb");
    if (fileNU == NULL) printf("File doesnt exists\n");


    info=fread(method->m, sizeof(val_type),ns*ns,fileM);
    if (info == -1) printf("Error fread command\n");
    info=fread(method->a, sizeof(val_type),ns*ns,fileA);
    if (info == -1) printf("Error fread command\n");
    info=fread(method->b, sizeof(val_type),ns,fileB);
    if (info == -1) printf("Error fread command\n");
    info=fread(method->c, sizeof(val_type),ns,fileC);
    if (info == -1) printf("Error fread command\n");
    info=fread(method->nu, sizeof(val_type),ns*ns,fileNU);
    if (info == -1) printf("Error fread command\n");


/*---- Verify symplectic condition ---------------------------------------------*/

#    ifdef IOUT
     for (i=0; i<ns; i++)
     {
        printf("\n");
        for (j=0; j<ns; j++) printf ("%lg,", method->m[i*ns+j]+method->m[j*ns+i]-1.);
     }

     printf("\n");
#    endif

/*---- Calculate hb coefficients -----------------------------------------------*/

     sum=0.;
     for (i=1; i<ns-1; i++)
     {
        method->hb[i]=h*method->b[i];
        sum+=method->hb[i];
     }
     
     method->hb[0]=(h-sum)/2.;
     method->hb[ns-1]=(h-sum)/2.;

/*---- Calculate hc coefficients -----------------------------------------------*/

     for (i=0; i<ns; i++)
     {
        method->hc[i]=h*method->c[i];
     }

     fclose(fileM);
     fclose(fileA);
     fclose(fileB);
     fclose(fileC);
     fclose(fileNU);

     return;
}


