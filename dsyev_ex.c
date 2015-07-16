#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <limits.h>

/* DSYEV prototype */
extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
	double* w, double* work, int* lwork, int* info );
/* DSYEVX prototype */
extern void dsyevx_( char* jobz, char* range, char* uplo, int* n, double* a,
	int* lda, double* vl, double* vu, int* il, int* iu, double* abstol,
	int* m, double* w, double* z, int* ldz, double* work, int* lwork,
	int* iwork, int* ifail, int* info );

static inline void *xmalloc(const size_t size)
{
	void *ret = malloc(size);

	if (!ret && size) {
		fprintf(stderr, "error: cannot allocate %zu bytes\n", size);
		exit(1);
	}
	return ret;
}

/* create hardcoded matrix from the original example

   Results should be:
----- snip -----
input matrix size: 5

input matrix:
   1.96  -6.49  -0.47  -7.20  -0.65
   0.00   3.80  -6.39   1.50  -6.34
   0.00   0.00   4.17  -1.51   2.67
   0.00   0.00   0.00   5.70   1.80
   0.00   0.00   0.00   0.00  -7.10

memory malloced (kbytes): 1

The total number of eigenvalues found:  5

results, eigenvalues:
 -11.07  -6.23   0.86   8.87  16.09

results, eigenvectors:
  -0.30  -0.61   0.40  -0.37   0.49
  -0.51  -0.29  -0.41  -0.36  -0.61
  -0.08  -0.38  -0.66   0.50   0.40
  -0.00  -0.45   0.46   0.62  -0.46
  -0.80   0.45   0.17   0.31   0.16
----- snap -----
*/
static void hardcoded_matrix( int n, double* a )
{
#define TMP_N 5

	double A[TMP_N*TMP_N] = {
	1.96,  0.00,  0.00,  0.00,  0.00,
	-6.49,  3.80,  0.00,  0.00,  0.00,
	-0.47, -6.39,  4.17,  0.00,  0.00,
	-7.20,  1.50, -1.51,  5.70,  0.00,
	-0.65, -6.34,  2.67,  1.80, -7.10
	};

	if (n != TMP_N ) {
		fprintf(stderr, "error: n = %d, this example needs n = %d\n",
			n, TMP_N);
		exit(1);
	}

	memcpy(a, A, sizeof(double)*TMP_N*TMP_N);

#undef TMP_N
}

/* create a random matrix, fill only the upper triangle (should be sparse!) */
static void random_matrix_upper( int n, double* a )
{
	int i,j;
	for (j=0; j < n; j++) {
		for (i=0; i <= j; i++) {
			a[(size_t)j * n + i] = (double) random()/RAND_MAX;
		}
	}
}

static void print_matrix( int m, int n, double* a )
{
	int i, j;
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) {
			printf( " %6.2f", a[(size_t)j * m + i] );
		}
		printf( "\n" );
	}
}

/* Parameters */
#define EXPERT 1

/* Main program */
int main(int argc, char **argv) {
	/* Locals */
	char* jobz = "V";
	int n = 20;
	int lda;
	int info;
	int lwork;
	double wkopt;
	double* work;
	double *w;
	double *a;
#ifdef EXPERT
	/* Set il, iu to compute NSELECT smallest eigenvalues */
	int il = 1;
	int iu = 10;
	int nselect;
	int m;
	int ldz;
	/* Negative abstol means using the default value */
	double abstol = -1.0;
	double vl;
	double vu;
	int* iwork;
	int* ifail;
	double* z;
#endif

	/* parse command line */
	if (argc > 1) {
		errno = 0;
		long int tmp = strtol(argv[1], (char **) NULL, 10);
		if (errno != 0 || tmp <= 0 || tmp > INT_MAX) {
			fprintf(stderr, "error: bad parameter n = '%s', "
				"should be 0 < n <= INT_MAX\n", argv[1]);
			exit(1);
		}
		n = tmp;
	}

	printf( "input matrix size: %d\n", n );

	lda = n;
	a = xmalloc(sizeof(double)*lda*n);
	w = xmalloc(sizeof(double)*n);
#ifdef EXPERT
	if (n < iu) {
		iu = n;
	}
	nselect = iu - il + 1;
	ldz = n;
	iwork = xmalloc(sizeof(int)*5*n);
	ifail = xmalloc(sizeof(int)*n);
	z = xmalloc(sizeof(double)*ldz*nselect);
#endif

	random_matrix_upper( n, a );
	//hardcoded_matrix( n, a );

	printf( "\ninput matrix:\n" );
	print_matrix( n, n, a );

	/* Query and allocate the optimal workspace */
	lwork = -1;
#ifndef EXPERT
	dsyev_( jobz, "Upper", &n, a, &lda, w, &wkopt, &lwork, &info );
#else
	dsyevx_( jobz, "Indices", "Upper", &n, a, &lda, &vl, &vu, &il, &iu,
		&abstol, &m, w, z, &ldz, &wkopt, &lwork, iwork, ifail, &info );
#endif

	lwork = (int)wkopt;
	work = (double*)xmalloc(sizeof(double)*lwork);
	{
		size_t malloced = 0;
		malloced += sizeof(double)*lda*n + sizeof(double)*n + sizeof(double)*lwork;
#ifdef EXPERT
		malloced += sizeof(int)*5*n + sizeof(int)*n + sizeof(double)*ldz*nselect;
#endif
		printf( "\nmemory malloced (kbytes): %zu\n", malloced / 1024 + 1);
	}
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
#else
		int p_m = m;
		double* p_z = z;
#endif
		printf( "\nThe total number of eigenvalues found: %2d\n", p_m );
		printf( "\nresults, eigenvalues:\n" );
		print_matrix( 1, p_m, w );
		printf( "\nresults, eigenvectors:\n" );
		print_matrix( n, p_m, p_z );
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
