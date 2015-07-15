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
/* DSYEVX prototype */
extern void dsyevx_( char* jobz, char* range, char* uplo, int* n, double* a,
	int* lda, double* vl, double* vu, int* il, int* iu, double* abstol,
	int* m, double* w, double* z, int* ldz, double* work, int* lwork,
	int* iwork, int* ifail, int* info );


void print_matrix( int m, int n, double* a, int lda ) {
	int i, j;
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) {
			printf( " %6.2f", a[(int64_t)i+j*lda] );
		}
		printf( "\n" );
	}
}

/* Parameters */
#define EXPERT 1

#define N 20
#define LDA N

#ifdef EXPERT
# define NSELECT 10
# define LDZ N
#endif

/* Main program */
int main() {
	/* Locals */
	int64_t i;
	char* jobz = "V";
	int n = N, lda = LDA, info, lwork;
	double wkopt;
	double* work;
#ifdef EXPERT
        int il, iu, m, ldz = LDZ;
        double abstol, vl, vu;
	int* iwork;
	int* ifail;
	double* z;
	/* Negative abstol means using the default value */
	abstol = -1.0;
	/* Set il, iu to compute NSELECT smallest eigenvalues */
	il = 1;
	iu = NSELECT;
	/* iwork dimension should be at least 5*n */
	iwork = malloc(sizeof(int)*5*N);
	ifail = malloc(sizeof(int)*N);
	z = malloc(sizeof(double)*LDZ*NSELECT);
#endif

#if 1
	double *w;
	double *a;

	a = malloc(sizeof(double)*LDA*N);
	w = malloc(sizeof(double)*N);

	for (i=0; i < (int64_t)LDA*N; i++) {
		a[i] = (double) random()/RAND_MAX;
	}
#endif
#if 0
	/* Local arrays */
	double w[N];
	double a[LDA*N] = {
	1.96,  0.00,  0.00,  0.00,  0.00,
	-6.49,  3.80,  0.00,  0.00,  0.00,
	-0.47, -6.39,  4.17,  0.00,  0.00,
	-7.20,  1.50, -1.51,  5.70,  0.00,
	-0.65, -6.34,  2.67,  1.80, -7.10
	};
#endif

	printf( "inout matrix:\n" );
	print_matrix( n, n, a, lda );

	/* Query and allocate the optimal workspace */
	lwork = -1;
#ifndef EXPERT
	dsyev_( jobz, "Upper", &n, a, &lda, w, &wkopt, &lwork, &info );
#else
	dsyevx_( jobz, "Indices", "Upper", &n, a, &lda, &vl, &vu, &il, &iu,
		&abstol, &m, w, z, &ldz, &wkopt, &lwork, iwork, ifail, &info );
#endif

	lwork = (int)wkopt;
	work = (double*)malloc(sizeof(double)*lwork);
	/* Solve eigenproblem */
#ifndef EXPERT
	dsyev_( jobz, "Upper", &n, a, &lda, w, work, &lwork, &info );
#else
	dsyevx_( jobz, "Indices", "Upper", &n, a, &lda, &vl, &vu, &il, &iu,
		&abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info );
#endif
	/* Check for convergence */
	if( info > 0 ) {
		fprintf( stderr, "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}

	{
#ifndef EXPERT
		int p_m = n;
		double* p_z = a;
		int p_ldz = lda;
#else
		int p_m = m;
		double* p_z = z;
		int p_ldz = ldz;
#endif
		printf( "\nThe total number of eigenvalues found:%2i\n", p_m );
		printf( "\nresults, eigenvalues:\n" );
		print_matrix( 1, p_m, w, 1 );
		printf( "\nresults, eigenvectors:\n" );
		print_matrix( n, p_m, p_z, p_ldz );
	}

	/* Free workspace */
	free( (void*)a );
	free( (void*)w );
	free( (void*)work );
#ifdef EXPERT
	free( (void*)iwork );
	free( (void*)ifail );
	free( (void*)z );
#endif
	exit( 0 );
}
