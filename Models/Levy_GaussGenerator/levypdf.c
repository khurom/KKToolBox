//******************************************************************************
//******************************************************************************
//	Filename:		levypdf.c
//	Function:		Program to compute and output a pdf of a Levy Process with
//					specified parameters gamma etc.
//	Author:			B. Hnat
//	Last edited:	10/03/2006 (K. Kiyani)
//******************************************************************************
//******************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define ABSV(a)  ((a)>0 ? (a) : (-a))
#define LOG2(a)  (((a)>1.0E-6) ? log((a))/log(2.0) : log(1.0E-6)/log(2.0))
#define LOG10(a) (((a)>1.0E-6) ? log10((a)) : log10(1.0E-6))
#define NEXP(a)  (((a)>1.0E-6) ? log10((a)) : log10(1.0E-6))

/* 1/f_s where f_s is a sampling frequncy */
#define NPOINT 100
#define INVPI  0.318309886
#define MAXQ   100
#define EPS    3.0e-11
#define DQ     1e-6
#define DX     0.1
#define X1     0.0
#define X2     10.0

FILE  *output;  
void gauleg(float x1, float x2, float x[], float w[], int n);
float func(float q, float gamma, float al, float x)
{
      double v1,v2;
      v1 = pow(q,al); v2 = q*x;
      return exp(-gamma*v1)*cos(v2);
}


/* This file will compute PDF of the difference based on the
given column in the data file. */
int main(int argc, const char * argv[])
{
    char      outfile[128], buf[1025], s[128];
    int       i, j, k, l;
    float     *xvec, *wvec;
    double    x0, x1, al, dp, v1,v2;
    double    x, y, q, dx, dq, gamma;
    
    /* Get output file name */
	printf("output file: %s", argv[1]);
    /*printf("\nOutput file: ");   scanf("%s", outfile);*/ sprintf(outfile,"%s",argv[1]);
    /*printf("Minimum x: ");       scanf("%s", buf);*/ x0 = atof(argv[2]);
    /*printf("Maximum x: ");       scanf("%s", buf);*/ x1 = atof(argv[3]);
    /*printf("Scaling index: ");   scanf("%s", buf);*/ al = atof(argv[4]);
    /*printf("Gamma parameter: "); scanf("%s", buf);*/ gamma = atof(argv[5]);

    /* Open input file and check number of columns */
    output = fopen(outfile, "w");

    xvec = (float*)malloc((NPOINT+1)*sizeof(float));
    wvec = (float*)malloc((NPOINT+1)*sizeof(float));

    printf("****************************\n");
    printf("** %f %f **\n", gamma, al);
    printf("** %f %f **\n", x0, x1);
    printf("****************************\n");

    fprintf(output,"%%# %e %e\n", gamma, al);
   /* Build Levy PDF(x) using Gauss-Legendre integration */
    x=x0-DX;
    while(x<x1)
    {  
       x+=DX; q=0; dp=0;
       gauleg(X1,X2,xvec,wvec,NPOINT);
       for(i=1;i<=NPOINT;i++) dp += (wvec[i]*func(xvec[i],gamma,al,x));
       y = INVPI*dp;
       fprintf(output,"%e %e\n", x, y);
    }

    printf("\n");
    fclose(output);  
    return 0;
}


//******************************************************************************
//	FUNCTIONS
//******************************************************************************

//------------------------------------------------------------------------------
//	gauleg	Given the lower and upper limits of integration 'x1' and 'x2', and 
//			given 'n', this routine returns arrays x[1..n] and w[1..n] 
//			containing the abscissas and weights of the Gauss-Legendre n-point
//			quadrature formula.
//------------------------------------------------------------------------------

void gauleg(float x1, float x2, float x[], float w[], int n)
{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++) {
		z=cos(3.141592654*(i-0.25)/(n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i]=xm-xl*z;
		x[n+1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i]=w[i];
	}
}
