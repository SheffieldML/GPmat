#ifndef LAPACK_H
#define LAPACK_H

#ifdef _MSC_VER
/* For MSVC I'm using clapack/cblas */
#include "blaswrap.h"
#endif


//#define dgetrf_ dgetrf
using namespace std;

// **** LAPACK Operations ****

// compute all eigenvalues and, optionally, eigenvectors of a real symmetric matrix A
#ifndef _NOSYEV
extern "C" void dsyev_(const char* jobz,
		       const char* uplo,
		       const int& n,  
		       double *a,   
		       const int& lda,
		       double *w,        // eigenvalues
		       double *work,
		       const int& lwork,
		       int &info);
#endif // ndef _NOSYEV	    
// Solve A*X=B for X ... i.e. X=A^-1*B.
extern "C" void dsysv_(const char* uplo,
		      const int& n,
		      const int& nrhs,
		      const double *a,
		      const int& lda,
		      const int *ipiv,
		      double *b,
		      const int& ldb,
		      double *work,
		      const int& lwork,
		      const int& info);
// Compute an LU factorization of a general M by N matrix A.
extern "C" void dgetrf_(
			const int &m,	 // (input)
			const int &n,	 // (input)
			double *a,	 // a[n][lda] (input/output)
			const int &lda,	 // (input)
			int *ipiv,	 // ipiv[min(m,n)] (output)
			int &info	 // (output)
			);
// Compute the inverse using the LU factorization computed by dgetrf.
extern "C" void dgetri_(
			const int &n,	 // (input)
			double *a,	 // a[n][lda] (input/output)
			const int &lda,	 // (input)
			const int *ipiv, // ipiv[n] (input)
			double *work,	 // work[lwork] (workspace/output)
			const int &lwork, // (input)
			int &info	 // (output)
			);
// Compute the Cholesky factorization of a real symmetric positive definite matrix A.
extern "C" void dpotrf_(
			const char* t,  // whether upper or lower triangluar 'U' or 'L'
			const int &n,	// (input)
			double *a,	// a[n][lda] (input/output)
			const int &lda,	// (input)
			int &info	// (output)
			);
// Compute the inverse of a real symmetric positive definite matrix A using the Cholesky factorization computed by dpotrf.
extern "C" void dpotri_(
			const char* t,  // whether upper or lower triangular 'U' or 'L'
			const int& n,   // input
			double *a,      // a[n][lda]
			const int &lda, // (input)
			int &info       // (output)
			);

// ***** BLAS Level 1 operations *****

// Perform y <-> x
extern "C" void dswap_(const int& n, 
		      double *x, 
		      const int& incx, 
		      double *y, 
		      const int& incy);
// Perform y:= x
extern "C" void dcopy_(const int& n, 
		      const double *x, 
		      const int& incx, 
		      double *y, 
		      const int& incy);
// Perform y:= ay
extern "C" void dscal_(const int& n, 
		      const double& alpha, 
		      double* y, 
		      const int& incy);

// Perform y := ax + y
extern "C" void daxpy_(const int& n, 
		       const double& alpha, 
		       const double *x, 
		       const int& incx, 
		       double *y, 
		       const int& incy);
// Return xTy
extern "C" double ddot_(const int& n,
			const double *x,
			const int& incx,
			const double *y,
			const int& incy);
// Return xTx
extern "C" double dnrm2_(const int& n,
			 const double *x,
			 const int& incx);
// ***** BLAS Level 2 operations *****

// Perform one of the matrix-vector operations y:= alpha*op(A)*x + beta*y
extern "C" void dgemv_(
                       const char* trans, // N, T or C transformation to A.
		       const int& m, // rows of A
		       const int& n, // columns of A
		       const double& alpha, // prefactor on multiplication
		       const double *A, // elements of A
		       const int& lda, // first dimension of A.
		       const double *x, // elements of x
		       const int& incx, // increment of x 
		       const double& beta, //prefactor on y.
		       double *y, // elements of y (output stored here).
		       const int& incy // increment of y.
		       );

// Perform one of the matrix-vector operations y:= alpha*A*x + beta*y for A symmetric.
extern "C" void dsymv_(
		       const char* t, // whether upper or lower triangular stored.
		       const int& n, // order of A.
		       const double& alpha, // prefactor on multiplcation.
		       const double* A, // elements of A.
		       const int& lda, // first dimension of A.
		       const double *x, // elements of x.
		       const int& incx, // increment of x.
		       const double& beta, // prefactor on y.
		       double *y, // elements of y (output stored here).
		       const int& incy // increment of y.
		       );
// perform a rank one update of the matrix A, A:= alpha*xy' + A
extern "C" void dger_(const int& m, //
		      const int& n,
		      const double& alpha,
		      const double *x,
		      const int& incx,
		      const double *y,
		      const int& incy,
		      const double *A,
		      const int& lda);

// perform a rank one update of the symmetrix matrix A, A:= alpha*xx' + A
extern "C" void dsyr_(const char* type, //
		      const int& n,
		      const double& alpha,
		      const double *x,
		      const int& incx,
		      const double *A,
		      const int& lda);

// ***** BLAS Level 3 operations *****

// Perform one of the matrix-matrix operations C:= alpha*op(A)*op(B) + beta*C
extern "C" void dgemm_(
		       const char* transa, // N, T or C transformation to A.
		       const char* transb, //
		       const int &m, // rows of A
		       const int &n, // columns of B
		       const int &k, // columns of A rows of B.
		       const double &alpha, // prefactor on multiplication
		       const double *A, // elements of A
		       const int &lda, // first dimension of A.
		       const double *B, // elements of B
		       const int &ldb, // first dimension of B.
		       const double &beta, //prefactor on C.
		       double *C, // elements of C (output stored here).
		       const int &ldc // first dimension of C.
		       );


// Perform one of the symmetric rank K operations C:=alpha*A*A' + beta*C.
extern "C" void dsyrk_(
		       const char* type, 
		       const char* trans, 
		       const int& n, 
		       const int& k,
		       const double& alpha, 
		       const double* A, 
		       const int& lda,
		       const double& beta,
		       const double* C,
		       const int& ldc);

// Perform triangular matrix matrix operation.
extern "C" void dtrmm_(const char* side,
		       const char* type,
		       const char* trans,
		       const char* diag,
		       const int& m,
		       const int& n,
		       const double& alpha,
		       const double* A,
		       const int& lda,
		       double* B,
		       const int& ldb);
// Perform inverse triangular matrix matrix operation.
extern "C" void dtrsm_(const char* side,
		       const char* type,
		       const char* trans,
		       const char* diag,
		       const int& m,
		       const int& n,
		       const double& alpha,
		       const double* A,
		       const int& lda,
		       double* B,
		       const int& ldb);
		       
// Perform symmetric matrix matrix operation.
extern "C" void dsymm_(const char* side,
		       const char* uplo,
		       const int& m,
		       const int& n,
		       const double& alpha,
		       const double* A,
		       const int& lda,
		       const double* B,
		       const int& ldb,
		       const double& beta,
		       double* C,
		       const int& ldc);

#endif
