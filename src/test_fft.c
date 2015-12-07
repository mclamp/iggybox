 /*  Based on :
  
 ============================================================================
 
 fourierd.c  -  Don Cross <dcross@intersrv.com>
 
 http://www.intersrv.com/~dcross/fft.html
 
 Contains definitions for doing Fourier transforms
 and inverse Fourier transforms.
 
 This module performs operations on arrays of 'double'.
 
 Revision history:
 
 1998 September 19 [Don Cross]
 Updated coding standards.
 Improved efficiency of trig calculations.
 
 ============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "datamodel.h"

#define  DDC_PI  (3.14159265358979323846)
#define CHECKPOINTER(p)  CheckPointer(p,#p)
#define FALSE 0
#define TRUE 1
int DFT(int dir,int m,double *x1,double *y1);
short FFT(short int dir,long m,double *x,double *y);

int IsPowerOfTwo ( unsigned x );
unsigned NumberOfBitsNeeded ( unsigned PowerOfTwo );
unsigned ReverseBits ( unsigned index, unsigned NumBits );
double Index_to_frequency ( unsigned NumSamples, unsigned Index );

void fft_float (
                unsigned  NumSamples,
                int       InverseTransform,
                float    *RealIn,
                float    *ImagIn,
                float    *RealOut,
                float    *ImagOut );
  



static void CheckPointer ( void *p, char *name )
{
  if ( p == NULL )
  {
    fprintf ( stderr, "Error in fft_double():  %s == NULL\n", name );
    exit(1);
  }
}

extern char **read_tokens(FILE *file,char tok, int *ntok);

int main(int argc, char**argv) {

  int p   = 9;
  int len = pow(2,p);
  
	double *inreal  = (double *)malloc(len*sizeof(double));
	double *inimag  = (double *)malloc(len*sizeof(double));
  
  double *outreal = (double *)malloc(len*sizeof(double));
	//double *outimag = (double *)malloc(len*sizeof(double));
  
  if (argc < 1) {
  //  printf("Usage: ./test_fft <datafile>\n");
   // exit(0);
  }
//  char *infile = "/Users/mclamp/cvs/pogview/src/woof2s";
  
  char *infile = argv[1];
  FILE *file = fopen(infile,"r");
  
  if (file == NULL) {
    printf("ERROR: Can't open file %s\n",infile);
    exit(0);
  }
  
  int i = 0;
  char **tokens = NULL;
  int ntok;
   
  while ((tokens = read_tokens(file,'\t',&ntok)) != NULL) {
    printf("Token %d\t%s\t%d\n",i,tokens[0],ntok);
    inreal[i] = atof(tokens[0]);
    
    int j = 0;
    
    while (j < ntok) {
      free(tokens[j]);
      j++;
    }
    ntok = 0;
    i++;
  }
  
  printf("i = %d\n",i);
  while (i > len) {
    printf("Mess %d %d %d\n",i,len,p);
    p = p+1;
    len = pow(2,p);
  }
  
  int chunk  = i;
  int offset = i-2;
  
  while (i < len) {
    inreal[i] = inreal[i-offset];
    
    if (i-offset > chunk) {
      offset += chunk;
    }
    i++;
  }
  
  
  //int i = 0;
	
	//while (i < len) {
		// I want 200 to be 2*pi
  //  float val = cos( i * 2 * DDC_PI / 10);
  //  inreal[i] = val;
  //  inimag[i] = 0.0;

//		i++;
//	}
	
  //fft_float(len,1,inreal,inimag,outreal,outimag);
  FFT(1,p,inreal,inimag);
  
  i = 0;
  
  while (i < len) {
    printf("inreal\t%d\t%f\n",i,inreal[i]);
    i++;
  }
  printf("\n");
         
  i = 0;
         
  while (i < len) {
    printf("outreal\t%d\t%f\n",i,outreal[i]);
    i++;
  }
  i = 0;
  while (i < len) {
      //printf("outimag\t%d\t%7.2f\n",i,outimag[i]);
    i++;
  }
  
  printf("\n");
                
  return 0;
}

/*
 Direct fourier transform
 */
int DFT(int dir,int m,double *x1,double *y1)
{
  long i,k;
  double arg;
  double cosarg,sinarg;
  double *x2=NULL,*y2=NULL;
  
  x2 = malloc(m*sizeof(double));
  y2 = malloc(m*sizeof(double));
  if (x2 == NULL || y2 == NULL)
    return(FALSE);
  
  for (i=0;i<m;i++) {
    x2[i] = 0;
    y2[i] = 0;
    arg = - dir * 2.0 * 3.141592654 * (double)i / (double)m;
    for (k=0;k<m;k++) {
      cosarg = cos(k * arg);
      sinarg = sin(k * arg);
      x2[i] += (x1[k] * cosarg - y1[k] * sinarg);
      y2[i] += (x1[k] * sinarg + y1[k] * cosarg);
    }
  }
  
  /* Copy the data back */
  if (dir == 1) {
    for (i=0;i<m;i++) {
      x1[i] = x2[i] / (double)m;
      y1[i] = y2[i] / (double)m;
    }
  } else {
    for (i=0;i<m;i++) {
      x1[i] = x2[i];
      y1[i] = y2[i];
    }
  }
  
  free(x2);
  free(y2);
  return(TRUE);
}

/*
 This computes an in-place complex-to-complex FFT 
 x and y are the real and imaginary arrays of 2^m points.
 dir =  1 gives forward transform
 dir = -1 gives reverse transform 
 */
short FFT(short int dir,long m,double *x,double *y)
{
  long n,i,i1,j,k,i2,l,l1,l2;
  double c1,c2,tx,ty,t1,t2,u1,u2,z;
  
  /* Calculate the number of points */
  n = 1;
  for (i=0;i<m;i++) 
    n *= 2;
  
  /* Do the bit reversal */
  i2 = n >> 1;
  j = 0;
  for (i=0;i<n-1;i++) {
    if (i < j) {
      tx = x[i];
      ty = y[i];
      x[i] = x[j];
      y[i] = y[j];
      x[j] = tx;
      y[j] = ty;
    }
    k = i2;
    while (k <= j) {
      j -= k;
      k >>= 1;
    }
    j += k;
  }
  
  /* Compute the FFT */
  c1 = -1.0; 
  c2 = 0.0;
  l2 = 1;
  for (l=0;l<m;l++) {
    l1 = l2;
    l2 <<= 1;
    u1 = 1.0; 
    u2 = 0.0;
    for (j=0;j<l1;j++) {
      for (i=j;i<n;i+=l2) {
        i1 = i + l1;
        t1 = u1 * x[i1] - u2 * y[i1];
        t2 = u1 * y[i1] + u2 * x[i1];
        x[i1] = x[i] - t1; 
        y[i1] = y[i] - t2;
        x[i] += t1;
        y[i] += t2;
      }
      z =  u1 * c1 - u2 * c2;
      u2 = u1 * c2 + u2 * c1;
      u1 = z;
    }
    c2 = sqrt((1.0 - c1) / 2.0);
    if (dir == 1) 
      c2 = -c2;
    c1 = sqrt((1.0 + c1) / 2.0);
  }
  
  /* Scaling for forward transform */
  if (dir == 1) {
    for (i=0;i<n;i++) {
      x[i] /= n;
      y[i] /= n;
    }
  }
  
  return(TRUE);
}
void fft_float (
                unsigned  NumSamples,
                int       InverseTransform,
                float    *RealIn,
                float    *ImagIn,
                float    *RealOut,
                float    *ImagOut ) {

  unsigned NumBits;    /* Number of bits needed to store indices */
  unsigned i, j, k, n;
  unsigned BlockSize, BlockEnd;
  
  double angle_numerator = 2.0 * DDC_PI;
  double tr, ti;     /* temp real, temp imaginary */
  
  if ( !IsPowerOfTwo(NumSamples) ) {
    fprintf (
             stderr,
             "Error in fft():  NumSamples=%u is not power of two\n",
             NumSamples );
    
    exit(1);
  }
  
  if ( InverseTransform ) {
    angle_numerator = -angle_numerator;
  }
  
  CHECKPOINTER ( RealIn );
  CHECKPOINTER ( RealOut );
  CHECKPOINTER ( ImagOut );
  
  NumBits = NumberOfBitsNeeded ( NumSamples );
  
  /*
   **   Do simultaneous data copy and bit-reversal ordering into outputs...
   */
  
  for ( i=0; i < NumSamples; i++ ) {
    j = ReverseBits ( i, NumBits );
    RealOut[j] = RealIn[i];
    ImagOut[j] = (ImagIn == NULL) ? 0.0 : ImagIn[i];
  }
  
  /*
   **   Do the FFT itself...
   */
  
  BlockEnd = 1;
  
  for ( BlockSize = 2; BlockSize <= NumSamples; BlockSize <<= 1 ) {
    double delta_angle = angle_numerator / (double)BlockSize;
    double sm2 = sin ( -2 * delta_angle );
    double sm1 = sin ( -delta_angle );
    double cm2 = cos ( -2 * delta_angle );
    double cm1 = cos ( -delta_angle );
    double w = 2 * cm1;
    double ar[3], ai[3];
    
    for ( i=0; i < NumSamples; i += BlockSize ) {

      ar[2] = cm2;
      ar[1] = cm1;
      
      ai[2] = sm2;
      ai[1] = sm1;
      
      for ( j=i, n=0; n < BlockEnd; j++, n++ ) {

        ar[0] = w*ar[1] - ar[2];
        ar[2] = ar[1];
        ar[1] = ar[0];
        
        ai[0] = w*ai[1] - ai[2];
        ai[2] = ai[1];
        ai[1] = ai[0];
        
        k = j + BlockEnd;
        tr = ar[0]*RealOut[k] - ai[0]*ImagOut[k];
        ti = ar[0]*ImagOut[k] + ai[0]*RealOut[k];
        
        RealOut[k] = RealOut[j] - tr;
        ImagOut[k] = ImagOut[j] - ti;
        
        RealOut[j] += tr;
        ImagOut[j] += ti;
      }
    }
    
    BlockEnd = BlockSize;
  }
  
  /*
   **   Need to normalize if inverse transform...
   */
  
  if ( InverseTransform ) {
    double denom = (double)NumSamples;
    
    for ( i=0; i < NumSamples; i++ ) {
      RealOut[i] /= denom;
      ImagOut[i] /= denom;
    }
  }
}

int IsPowerOfTwo ( unsigned x ) {
  if ( x < 2 )
    return 0;
  
  if ( x & (x-1) )        // Thanks to 'byang' for this cute trick!
    return 0;
  
  return 1;
}

unsigned NumberOfBitsNeeded ( unsigned PowerOfTwo ) {
  unsigned i;
  
  if ( PowerOfTwo < 2 ){
    fprintf (
             stderr,
             ">>> Error in fftmisc.c: argument %d to NumberOfBitsNeeded is too small.\n",
             PowerOfTwo );
    
    exit(1);
  }
  
  for ( i=0; ; i++ ) {
    if ( PowerOfTwo & (1 << i) )
      return i;
  }
}

unsigned ReverseBits ( unsigned index, unsigned NumBits ) {
  unsigned i, rev;
  
  for ( i=rev=0; i < NumBits; i++ ) {
    rev = (rev << 1) | (index & 1);
    index >>= 1;
  }
  
  return rev;
}

double Index_to_frequency ( unsigned NumSamples, unsigned Index ) {

  if ( Index >= NumSamples ) {
    return 0.0;
  } else if ( Index <= NumSamples/2 ) {
    return (double)Index / (double)NumSamples;
  }
  
  return -(double)(NumSamples-Index) / (double)NumSamples;
}









