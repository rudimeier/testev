/*******************************************************************************
*  Copyright (C) 2009-2014 Intel Corporation. All Rights Reserved.
*  The information and material ("Material") provided below is owned by Intel
*  Corporation or its suppliers or licensors, and title to such Material remains
*  with Intel Corporation or its suppliers or licensors. The Material contains
*  proprietary information of Intel or its suppliers and licensors. The Material
*  is protected by worldwide copyright laws and treaty provisions. No part of
*  the Material may be copied, reproduced, published, uploaded, posted,
*  transmitted, or distributed in any way without Intel's prior express written
*  permission. No license under any patent, copyright or other intellectual
*  property rights in the Material is granted to or conferred upon you, either
*  expressly, by implication, inducement, estoppel or otherwise. Any license
*  under such intellectual property rights must be express and approved by Intel
*  in writing.
*
********************************************************************************
*/
/*
   DSYEV Example.
   ==============

   Program computes all eigenvalues and eigenvectors of a real symmetric
   matrix A:

     1.96  -6.49  -0.47  -7.20  -0.65
    -6.49   3.80  -6.39   1.50  -6.34
    -0.47  -6.39   4.17  -1.51   2.67
    -7.20   1.50  -1.51   5.70   1.80
    -0.65  -6.34   2.67   1.80  -7.10

   Description.
   ============

   The routine computes all eigenvalues and, optionally, eigenvectors of an
   n-by-n real symmetric matrix A. The eigenvector v(j) of A satisfies

   A*v(j) = lambda(j)*v(j)

   where lambda(j) is its eigenvalue. The computed eigenvectors are
   orthonormal.

   Example Program Results.
   ========================

 DSYEV Example Program Results

 Eigenvalues
 -11.07  -6.23   0.86   8.87  16.09

 Eigenvectors (stored columnwise)
  -0.30  -0.61   0.40  -0.37   0.49
  -0.51  -0.29  -0.41  -0.36  -0.61
  -0.08  -0.38  -0.66   0.50   0.40
   0.00  -0.45   0.46   0.62  -0.46
  -0.80   0.45   0.17   0.31   0.16
*/
#include <stdlib.h>
#include <stdio.h>

/* DSYEV prototype */
extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
	double* w, double* work, int* lwork, int* info );

void print_matrix( int m, int n, double* a, int lda ) {
	int i, j;
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) {
			printf( " %6.2f", a[i+j*lda] );
		}
		printf( "\n" );
	}
}

/* Parameters */
#define N 5
#define LDA N

/* Main program */
int main() {
	/* Locals */
	int n = N, lda = LDA, info, lwork;
	double wkopt;
	double* work;
	/* Local arrays */
	double w[N];
	double a[LDA*N] = {
	1.96,  0.00,  0.00,  0.00,  0.00,
	-6.49,  3.80,  0.00,  0.00,  0.00,
	-0.47, -6.39,  4.17,  0.00,  0.00,
	-7.20,  1.50, -1.51,  5.70,  0.00,
	-0.65, -6.34,  2.67,  1.80, -7.10
	};
	
	printf( "inout matrix:\n" );
	print_matrix( n, n, a, lda );
	
	/* Query and allocate the optimal workspace */
	lwork = -1;
	dsyev_( "V", "U", &n, a, &lda, w, &wkopt, &lwork, &info );
	lwork = (int)wkopt;
	work = (double*)malloc( lwork*sizeof(double) );
	/* Solve eigenproblem */
	dsyev_( "V", "U", &n, a, &lda, w, work, &lwork, &info );
	/* Check for convergence */
	if( info > 0 ) {
		fprintf( stderr, "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}
	
	printf( "\nresults, eigenvalues:\n" );
	print_matrix( 1, n, w, 1 );
	printf( "\nresults, eigenvectors:\n" );
	print_matrix( n, n, a, lda );

	/* Free workspace */
	free( (void*)work );
	exit( 0 );
}
