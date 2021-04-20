#include "driver.h"

#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

void LoopFive( int, int, int, double *, int, double *, int, double *, int );
void LoopFour( int, int, int, double *, int, double *, int,  double *, int );
void LoopThree( int, int, int, double *, int, double *, double *, int );
void LoopTwo( int, int, int, double *, double *, double *, int );
void LoopOne( int, int, int, double *, double *, double *, int );
void Gemm_MRxNRKernel_Packed( int, double *, double *, double *, int );
void PackBlockA( int, int, double *, int, double * );
void PackPanelB( int, int, double *, int, double * );

/* Blocking parameters */
#define MC -1
#define MR -1
#define NC -1
#define NR -1
#define KC -1
  
/* Wrapper for GEMM function */
void MyGemm( int m, int n, int k, double *A, int ldA,
	     double *B, int ldB, double *C, int ldC )
{

  if (MC == -1 || MR == -1 || NC == -1 || NR == -1 || KC == -1 )
  {
    printf("Some of the blocking parameters are not set\n");
    exit(0);
  }
  if ( m % MR != 0 || MC % MR != 0 ){
    printf( "m and MC must be multiples of MR\n" );
    exit( 0 );
  }
  if ( n % NR != 0 || NC % NR != 0 ){
    printf( "n and NC must be multiples of NR\n" );
    exit( 0 );
  }

  LoopFive( m, n, k, A, ldA, B, ldB, C, ldC );
}

void LoopFive(  int m, int n, int k, 
                double *A, int ldA,
		            double *B, int ldB, 
                double *C, int ldC )
{
  for ( int j=0; j<n; j+=NC ) 
  {
    int jb = dmin( NC, n-j );    /* Last loop may not involve a full block */
    LoopFour( m, jb, k, A, ldA, &beta( 0,j ), ldB, &gamma( 0,j ), ldC );
  } 
}

void LoopFour(  int m, int n, int k, 
                double *A, int ldA, 
                double *B, int ldB,
	              double *C, int ldC )
{
  double *Btilde = ( double * ) _mm_malloc( KC * NC * sizeof( double ), 64 );
  
  for ( int p=0; p<k; p+=KC ) 
  {
    int pb = dmin( KC, k-p );    /* Last loop may not involve a full block */
    PackPanelB( pb, n, &beta( p, 0 ), ldB, Btilde );
    LoopThree( m, n, pb, &alpha( 0, p ), ldA, Btilde, C, ldC );
  }

  _mm_free( Btilde); 
}

void LoopThree( int m, int n, int k, 
                double *A, int ldA, 
                double *Btilde, 
                double *C, int ldC )
{
  double *Atilde = ( double * ) _mm_malloc( MC * KC * sizeof( double ), 64 );
       
  for ( int i=0; i<m; i+=MC ) {
    int ib = dmin( MC, m-i );    /* Last loop may not involve a full block */
    PackBlockA( ib, k, &alpha( i, 0 ), ldA, Atilde );
    LoopTwo( ib, n, k, Atilde, Btilde, &gamma( i,0 ), ldC );
  }

  _mm_free( Atilde);
}

void LoopTwo( int m, int n, int k, 
              double *Atilde, 
              double *Btilde, 
              double *C, int ldC )
{
  for ( int j=0; j<n; j+=NR ) {
    int jb = dmin( NR, n-j );
    LoopOne( m, jb, k, Atilde, &Btilde[ j*k ], &gamma( 0,j ), ldC );
  }
}

void LoopOne( int m, int n, int k, 
              double *Atilde, 
              double *MicroPanelB, 
              double *C, int ldC )
{
  for ( int i=0; i<m; i+=MR ) {
    int ib = dmin( MR, m-i );
    Gemm_MRxNRKernel_Packed( k, &Atilde[ i*k ], MicroPanelB, &gamma( i,0 ), ldC );
  }
}


/* DGEMM mircokernel 
  Computes C += AB where C is MR x NR, A is MR x KC, and B is KC x NR */
void Gemm_MRxNRKernel_Packed( int k, double *A, double *B, double *C, int ldC)
{

}

/* Pack a MC x KC block of A into Atilde */
void PackBlockA( int m, int k, double *A, int ldA, double *Atilde )
{

}

/* Pack a KC x NC block of B into Btilde */
void PackPanelB( int k, int n, double *B, int ldB, double *Btilde )
{

}
