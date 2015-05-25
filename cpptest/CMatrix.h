#ifndef CMATRIX_H
#define CMATRIX_H
// Base matrix class for Sheffield ML code.
#include <iostream>
#include <vector>
#include "ndlutil.h"
#include "CMatlab.h"
#include "lapack.h"
#include "xpose.h"
using namespace std;


class CMatrix : public CMatinterface {
 public:
  // The default constructor.
  CMatrix()
    {
      nrows = 1;
      ncols = 1;
      symmetric = false;
      triangular = false;
      memAllocate();
      
    }
  // A constructor for creating a 1x1 CMatrix from a double.
  CMatrix(double val)
    {
      nrows = 1;
      ncols = 1;
      symmetric = false;
      triangular = false;
      memAllocate();
      vals[0] = val;
    }
  // The standard memory allocating constructor for creating a matrix o f size numRows*numCols.
  CMatrix(const int numRows, const int numCols) : nrows(numRows), ncols(numCols)
    {
      symmetric = false;
      triangular = false;
      memAllocate();
    }
  // Constructor which allocates memory and then fills the CMatrix with constant values.
  CMatrix(const int numRows, const int numCols, const double val) : nrows(numRows), ncols(numCols)
    {
      symmetric = false;
      triangular = false;
      memAllocate();
      setVals(val);
    }
  // Constructor for initialising a CMatrix from a double* array.
  CMatrix(const int numRows, const int numCols, double* inVals) : nrows(numRows), ncols(numCols)
    {
      symmetric = false;
      triangular = false;
      memAllocate();
      dcopy_(nrows*ncols, inVals, 1, vals, 1);
    }
  // Constructor for initialising a CMatrix from a double** array.
  CMatrix(const int numRows, const int numCols, double** inVals) : nrows(numRows), ncols(numCols)
    {
      symmetric = false;
      triangular = false;
      memAllocate();
      for(int i=0; i<nrows; i++)
	for(int j=0; j<ncols; j++)
	  vals[i+nrows*j] = inVals[i][j];
    }
	
  // Constructor for special initialisations such as identity or random matrices.
  CMatrix(const int numRows, const int numCols, char type) : nrows(numRows), ncols(numCols)
    {
      setSymmetric(false);
      setSymmetric(false);
      memAllocate();
      switch(type) {
      case 'I':
	// the identity
	assert(numRows==numCols);
	for(int i=0; i<nrows; i++)
	  for(int j=0; j<ncols; j++)
	    if(i==j)
	      vals[i+nrows*j] = 1.0;
	    else
	      vals[i+nrows*j] = 0.0;
	setSymmetric(true);
	break;
      default:
	assert(0);
      }
    }
  // Constructor for special initialisations where a value is also passed.
  CMatrix(const int numRows, const int numCols, char type, double val) : nrows(numRows), ncols(numCols)
    {
      setSymmetric(false);
      setSymmetric(false);
      memAllocate();
      switch(type) {
      case 'S':
	// a spherical covariance matrix
	for(int i=0; i<nrows; i++)
	  for(int j=0; j<ncols; j++)
	    if(i==j)
	      vals[i+nrows*j] = val;
	    else
	      vals[i+nrows*j] = 0.0;
	setSymmetric(true);
	break;
      default:
	assert(0);
      }
    }
  // The copy constructor, it performs a deep copy.
  CMatrix(const CMatrix& A) : nrows(A.nrows), ncols(A.ncols), symmetric(A.symmetric), triangular(A.triangular)
    {
      memAllocate();
      copy(A);
      
    }
  // The class destructor, it deallocates the memory.
  ~CMatrix()
    {
      memDeAllocate();
    }
  // Perform a deep copy of the matrix A, resizing if necessary.
  void deepCopy(const CMatrix& A)
    {
      resize(A.nrows, A.ncols);
      copy(A);
    }
  // Get the number of rows in the matrix.
  inline const int getRows() const
    {
      return nrows;
    }
  // Get the number of columns in the matrix.
  inline const int getCols() const
    {
      return ncols;
    }
  // Get the number of elements in the matrix (rows*cols).
  inline const int getNumElements() const
    {
      return nrows*ncols;
    }
  // Get the matrix element in the ith row and jth column (indexing from 0).
  inline const double getVal(const int i, const int j) const
    {
      assert(i>=0 && i<nrows);
      assert(j>=0 && j<ncols);
      return vals[i + nrows*j];
    } 
  // Get the ith element from the matrix.
  inline const double getVal(const int i) const
    {
      assert(i>=0 && i<nrows*ncols);
      return vals[i];
    }
  // Set all elements of the matrix to val.
  inline void setVals(const double val)
    {
      for(int i=0; i<nrows*ncols; i++)
	vals[i] = val;
    }
  // Set the ith element of the matrix to val.
  inline void setVal(const double val, const int i)
    {
      assert(i>=0 && i<nrows*ncols);
      vals[i] = val;
    }
  // Set the matrix element from the ith row and jth column to val.
  inline void setVal(const double val, const int i, const int j)
    {
      assert(i>=0 && i<nrows);
      assert(j>=0 && j<ncols);
      vals[i + nrows*j] = val;
    }
  // Returns true if the matrix has the same number of rows as columns.
  inline const bool isSquare() const
    {
      return nrows==ncols;
    }
  // Returns true if the matrix is symmetric.
  inline const bool isTriangular() const
    {
      return triangular;
    }
  inline void setTriangular(const bool val) 
    {
      assert((val && isSquare()) || !val);	
      triangular=val;
    }
  inline const bool isSymmetric() const
    {
      return symmetric;
    }
  inline void setSymmetric(const bool val) 
    {
      assert((val && isSquare()) || !val);	
      symmetric=val;
    }
  // Returns true if the matrix A has the same dimensions as the matrix.
  inline const bool dimensionsMatch(const CMatrix& A) const
    {
      return (nrows==A.nrows && ncols==A.ncols);
    }
  // Returns true if A has the same number of rows as the matrix.
  inline const bool rowsMatch(const CMatrix& A) const
    {
      return (nrows==A.nrows);
    }
  // Returns true if A has the same number of columns as the matrix.
  inline const bool colsMatch(const CMatrix& A) const
    {
      return (ncols==A.ncols);
    }
  // copy the upper part to the lower or vice versa.
  void copySymmetric(const char* type);
  void copyRowRow(const int i, const CMatrix& X, const int k);
  // Scale the matrix by a constant alpha.
  void scale(const double alpha)
    {
      dscal_(nrows*ncols, alpha, vals, 1); 
    }
  // Scale the jth column of the matrix.
  void scaleCol(const int j, const double alpha)
    {
      assert(j<ncols);
      dscal_(nrows, alpha, vals+j*nrows, 1);
    }
  // Scale the ith row of the matrix.
  void scaleRow(const int i, const double alpha)
    {
      assert(i<nrows);
      dscal_(ncols, alpha, vals+i, nrows);
    }
  // Level 1 BLAS axpy  y c:= alpha x + y
  void axpy(const CMatrix& x, const double alpha)
  {
    assert(x.ncols==ncols);
    assert(x.nrows==nrows);
    daxpy_(ncols*nrows, alpha, x.vals, 1, vals, 1);
  }
  // Level 1 BLAS axpy y(i, :) := alpha*x(k, :) + y(i, :);
  void axpyRowRow(const int i, const CMatrix& A, const int k, const double alpha)
    {
      assert(i<nrows);
      assert(k<A.nrows);
      assert(A.ncols==ncols);
      daxpy_(ncols, alpha, A.vals+k, A.nrows, vals+i, nrows);
    }
  // Level 1 BLAS axpy (i, :) := alpha*x(:, j)' + y(i, :);
  void axpyRowCol(const int i, const CMatrix& A, const int j, const double alpha)
    {
      assert(i<nrows);
      assert(j<A.ncols);
      assert(A.nrows==ncols);
      daxpy_(ncols, alpha, A.vals+j*A.nrows, 1, vals+i, nrows);
    }
  // Level 1 BLAS axpy (:, j) := alpha*x(:, k) + y(:, j);
  void axpyColCol(const int j, const CMatrix& A, const int k, const double alpha)
    {
      assert(j<ncols);
      assert(k<A.ncols);
      assert(A.nrows==nrows);
      daxpy_(nrows, alpha, A.vals+k*A.nrows, 1, vals+j*nrows, 1);
    }
  // Level 1 BLAS axpy (:, j) = alpha*x(i, :)' + y(:, j);
  void axpyColRow(const int j, const CMatrix& A, const int i, const double alpha)
    {
      assert(j<ncols);
      assert(i<A.nrows);
      assert(A.ncols==nrows);
      daxpy_(nrows, alpha, A.vals+i, A.nrows, vals+j*nrows, 1);
    }
  // Level 1 BLAS axpy diag(Y) = diag(Y) + alpha*x(i, :)'
  void axpyDiagRow(const CMatrix& A, const int i, const double alpha)
    {
      assert(isSquare());
      assert(i<A.nrows);
      assert(A.ncols==nrows);
      daxpy_(nrows, alpha, A.vals+i, A.nrows, vals, nrows+1);
    }
  // Level 1 BLAS axpy diag(Y) = diag(Y) + alpha*x(i, :)'
  void axpyDiagCol(const CMatrix& A, const int j, const double alpha)
    {
      assert(isSquare());
      assert(j<A.ncols);
      assert(A.nrows==nrows);
      daxpy_(nrows, alpha, A.vals+j*A.nrows, 1, vals, nrows+1);
    }
  // Level 2 BLAS Rank 1 update: ger, A = alpha*x*y' + A; 
  void ger(const CMatrix& x, const CMatrix& y, const double alpha)
  {
    assert(x.ncols==1);
    assert(y.ncols==1);
    assert(x.nrows==nrows);
    assert(y.nrows==ncols);
    dger_(nrows, ncols, alpha, x.vals, 1, y.vals, 1, vals, nrows);
  }
  // Level 2 BLAS Rank 1 update: A := alpha*x(k, :)'*y(i, :) + A;
  void gerRowRow(const CMatrix& x, const int i, const CMatrix& y, const int k, const double alpha)
    {
      assert(i<x.nrows);
      assert(k<y.nrows);
      assert(x.ncols==nrows);
      assert(y.ncols==ncols);
      dger_(nrows, ncols, alpha, x.vals+i, x.nrows, y.vals+k, y.nrows, vals, nrows);
    }
  // Level 2 BLAS Rank 1 update: A := alpha*x(:, j)*y(i, :) + A;
  void gerRowCol(const CMatrix& x, const int i, const CMatrix& y, int j, const double alpha)
    {
      assert(i<x.nrows);
      assert(j<y.ncols);
      assert(x.ncols==nrows);
      assert(y.nrows==ncols);
      dger_(nrows, ncols, alpha, x.vals+i, x.nrows, y.vals+j*y.nrows, 1, vals, nrows);
    }
  // Level 2 BLAS Rank 1 update: A := alpha*x(:, k)*y(:, j)' + A;
  void gerColCol(const CMatrix& x, const int j, const CMatrix& y, const int k, const double alpha)
    {
      assert(j<x.ncols);
      assert(k<y.ncols);
      assert(x.nrows==nrows);
      assert(y.nrows==ncols);
      dger_(nrows, ncols, alpha, x.vals+j*x.nrows, 1, y.vals+k*y.nrows, 1, vals, nrows);
    }
  // Level 2 BLAS Rank 1 update: A := alpha*x(i, :)'x(:, j)' + A;
  void gerColRow(const CMatrix& x, int j, const CMatrix& y, const int i, const double alpha)
    {
      assert(j<x.ncols);
      assert(i<y.nrows);
      assert(x.nrows==nrows);
      assert(y.ncols==ncols);
      dger_(nrows, ncols, alpha, x.vals+j*x.nrows, 1, y.vals+i, y.nrows, vals, nrows);
    }

  // Level 2 BLAS Rank 1 update: syr, A = alpha*x*x' + A; 
  void syr(const CMatrix& x, const double alpha, const char* trans)
  {
    assert(isSymmetric());
    assert(x.ncols==1);
    assert(x.nrows==nrows);
    dsyr_(trans, nrows, alpha, x.vals, 1, vals, nrows);
    copySymmetric(trans);
  }
  // Level 2 BLAS Rank 1 update: A := alpha*x(i, :)'*x(i, :) + A;
  void syrRow(const CMatrix& x, const int i, const double alpha, const char* trans)
    {
      assert(isSymmetric());
      assert(i<x.nrows);
      assert(x.ncols==nrows);
      dsyr_(trans, nrows, alpha, x.vals+i, x.nrows, vals, nrows);
      copySymmetric(trans);
    }
  // Level 2 BLAS Rank 1 update: A := alpha*x(:, j)x(:, j)' + A;
  void syrCol(const CMatrix& x, const int j, const double alpha, const char* trans)
    {
      assert(isSymmetric());
      assert(j<x.ncols);
      assert(x.ncols==nrows);
      dsyr_(trans, nrows, alpha, x.vals+j*x.nrows, 1, vals, nrows);
      copySymmetric(trans);
    }

  // Return the euclidean distance between two row vectors.
  double dist2Row(const int i, const CMatrix& A, const int k) const
    {
      assert(ncols==A.ncols);
      assert(k<A.nrows);
      assert(i<nrows);
      return norm2Row(i) + A.norm2Row(k) - 2.0*dotRowRow(i, A, k);
    }
  // Return the euclidean distance between two column vectors.
  double dist2Col(const int j, const CMatrix& A, const int k) const
    {
      assert(nrows==A.nrows);
      assert(k<A.ncols);
      assert(j<ncols);
      return norm2Col(j) + A.norm2Col(k) - 2.0*dotColCol(j, A, k);
    }
  // Return the norm of the ith row of the matrix.
  double normRow(const int i) const
    {
      return dnrm2_(ncols, vals+i, nrows);
    }
  // Return the squared norm of the ith row of the matrix.
  double norm2Row(const int i) const
    {
      double val=dnrm2_(ncols, vals+i, nrows);
      return val*val;
    }
  // Return the norm of the jth column of the matrix.
  double normCol(const int j) const
    {
      return dnrm2_(nrows, vals+j*nrows, 1);
    }
  // Return the squared norm of the jth column of the matrix.
  double norm2Col(const int j) const
    {
      double val=dnrm2_(nrows, vals+j*nrows, 1);
      return val*val;
    }
  // Return the inner product between the ith row of the matrix and the kth row of A.
  double dotRowRow(const int i, const CMatrix& A, const int k) const
    {
      return ddot_(ncols, A.vals+k, A.nrows, vals+i, nrows);
  
    }
  // Return the inner product between the ith row of the matrix and the jth column of A.
  double dotRowCol(const int i, const CMatrix& A, const int j) const
    {
      return ddot_(ncols, A.vals+j*A.nrows, 1, vals+i, nrows);
  
    }
  // Return the inner product between the jth column of the matrix and the kth column of A.
  double dotColCol(const int j, const CMatrix& A, const int k) const
    {
      return ddot_(nrows, A.vals+k*A.nrows, 1, vals+j*nrows, 1);
    }
  // Return the inner product between the jth column of the matrix and the ith row of A.
  double dotColRow(const int j, const CMatrix& A, const int i) const
    {
      return ddot_(nrows, A.vals+i, A.nrows, vals+j*nrows, 1);
    }
  // Swap the jth and the kth columns of the matrix.
  void swapCols(const int j, const int k)
    {
      assert(j<ncols && k<ncols);
      if(j!=k)
	dswap_(nrows, vals+j*nrows, 1, vals+k*nrows, 1);
    }
  // Swap the ith and the kth rows of the matrix.
  void swapRows(const int i, const int k)
    {
      assert(i<nrows && k<nrows);
      if(i!=k)
	dswap_(ncols, vals+i, nrows, vals+k, nrows);
    }
  // Add columns from A to the end of the matrix.
  void appendCols(const CMatrix& A)
    {
      assert(rowsMatch(A));
      int origNcols = ncols;
      memReAllocate(0, A.ncols);
      setMatrix(0, origNcols, A);
    }
  // Add rows from A to the end of the matrix.
  void appendRows(const CMatrix& A)
    {
      assert(colsMatch(A));
      int origNrows = nrows;
      memReAllocate(A.nrows, 0);
      setMatrix(origNrows, 0, A);
    }
  // Get the rows firstRow:lastRow and columns firstCol:lastCol and place in a  matrix C.
  void getMatrix(CMatrix& C, const int firstRow, const int lastRow, const int firstCol, const int lastCol) const
    {
      assert(firstRow>=0 && firstRow<=lastRow && lastRow<nrows);
      assert(firstCol>=0 && firstCol<=lastCol && lastCol<ncols);
      assert(C.nrows==lastRow-firstRow+1 && C.ncols==lastCol-firstCol+1);
      for(int j=0; j<C.ncols; j++)
	for(int i=0; i<C.nrows; i++)
	  C.vals[i+C.nrows*j] = vals[i+firstRow + nrows*(j+firstCol)];
    }
  // Get the rows in rows and columns firstCol:lastCol and place in a matrix C.
  void getMatrix(CMatrix& C, vector<int> rows, int firstCol, int lastCol)
    {
      assert(firstCol>=0 && firstCol<=lastCol && lastCol<ncols);
      assert(C.nrows==rows.size() && C.ncols==lastCol-firstCol+1);
      for(int i=0; i<C.nrows; i++)
	{
	  assert(rows[i]<nrows);
	  for(int j=0; j<C.ncols; j++)
	    C.vals[i+C.nrows*j] = vals[rows[i] + nrows*(j+firstCol)];
	}
    }
  // Get the rows firstRow:lastRow and columns in cols and place in matrix C.
  void getMatrix(CMatrix& C, int firstRow, int lastRow, vector<int> cols)
    {
      assert(firstRow>=0 && firstRow<=lastRow && lastRow<nrows);
      assert(C.nrows==lastRow-firstRow+1 && C.ncols==cols.size());
      for(int j=0; j<C.ncols; j++)
	{
	  assert(cols[j]<ncols);
	  for(int i=0; i<C.nrows; i++)
	    C.vals[i+C.nrows*j] = vals[i+firstRow + nrows*(cols[j])];
	}
    }
  // Get the rows from rows and columns from cols and place in matrix C.
  void getMatrix(CMatrix& C, vector<int> rows, vector<int> cols)
    {
      assert(C.nrows==rows.size() && C.ncols==cols.size());
      for(int i=0; i<C.nrows; i++)
	{
	  assert(rows[i]<C.nrows);
	  for(int j=0; j<C.ncols; j++)
	    {
	      assert(cols[j]<ncols);
	      C.vals[i+C.nrows*j] = vals[rows[i] + nrows*(cols[j])];
	    }
	}
    }
  // Place A's first row and column at row, col and the rest of the matrix follows.
  void setMatrix(const int row, const int col, const CMatrix& A)
    {
      assert(row+A.nrows <= nrows);
      assert(col+A.ncols <= ncols);
      for(int i=0; i<A.nrows; i++)
	for(int j=0; j<A.ncols; j++)
	  vals[i+row+nrows*(j+col)] = A.vals[i+A.nrows*j];
    }
  // Place the rows of A at the locations given by rows starting at column col.
  void setMatrix(const vector<int> rows, const int col, const CMatrix& A)
    {
      assert(rows.size()==A.nrows);
      assert(col+A.ncols<=ncols);
      for(int i = 0; i<A.nrows; i++)
	for(int j = col; j<A.ncols+col; j++)
	  {
	    assert(rows[i]<nrows);
	    vals[rows[i]+nrows*j]=A.vals[i+A.nrows*j];
	  }
    }
  // Place the columns of A at the locations given by cols starting at row row.
  void setMatrix(const int row, const vector<int> cols, const CMatrix& A)
    {
      assert(cols.size()==A.ncols);
      assert(row+A.nrows<=nrows);
      for(int i=row; i<A.nrows+row; i++)
	for(int j=0; j<A.ncols; j++)
	  {
	    assert(cols[i]<ncols);	    
	    vals[i+nrows*cols[j]]=A.vals[i+A.nrows*j];
	  }
    }
  // Place the rows and columns of A at rows and cols.
  void setMatrix(const vector<int> rows, const vector<int> cols, const CMatrix& A)
    {
      assert(cols.size()==A.ncols);
      assert(rows.size()==A.nrows);
      for(int i=0; i<A.nrows; i++)
	for(int j=0; j<A.ncols; j++)
	  {
	    assert(rows[i]<nrows);
	    assert(cols[j]<ncols);
	    vals[rows[i]+nrows*cols[j]]=A.vals[i+A.nrows*j];
	  }
    }
  // In place transpose of the matrix using Algorithm 380.
  void dtrans(const int lwork)
    {
      int res;
      int* work = new int[lwork];
      dtrans_(vals, nrows, ncols, nrows*ncols, work, lwork, res);
      assert(res==0);
    }
  // In place transpose of the matrix using Algorithm 513.
  void dtransr(const int lwork)
    {
      int res;
      int* work = new int[lwork];
      dtransr_(vals, nrows, ncols, nrows*ncols, work, lwork, res);
      assert(res==0);
    }
  // In place transpose of the matrix using Algorithm 467 (currently not working).
  void dxpose(const int lwork)
    {
      // this is algorithm 467 (it is supposed to be more efficient - but doesn't work!)
      bool* work = new bool[lwork];
      dxpose_(vals, nrows, ncols, nrows*ncols, work, lwork);
    }
  // Perform matrix transpose.
  void trans()
    {
      // if rows or columns are 1 dimensional then you don't need to move elements.
      if (nrows!=1 && ncols!=1)
	dtransr((nrows+ncols)/2); // this is algorithm 513 (467 doesn't seem to work)
      // if the matrix is square then you don't need to swap rows and columns.
      if (~isSquare())
	{
	  int temp = nrows;
	  nrows = ncols;
	  ncols = temp;
	}
    }
  // Multiply the elements of the matrix by the elements of A.
  void multiply(const CMatrix& A)
    {
      // if A is a row or column vector it is `replicated' before the operation.
      if(A.nrows==1)
	{
	  assert(A.ncols==ncols);
	  for(int i=0; i<nrows; i++)
	    for(int j=0; j<ncols; j++)
	      vals[i+nrows*j] *= A.vals[j];
	}
      else if(A.ncols==1)
	{
	  assert(A.nrows==nrows);
	  for(int j=0; j<ncols; j++)
	    for(int i=0; i<nrows; i++)
	      vals[i+nrows*j] *= A.vals[i];
	}
      else  
	{
	  assert(A.nrows==nrows && A.ncols == ncols);
	  for(int i=0; i<nrows*ncols; i++)
	    vals[i] *= A.vals[i];
	}    
    }
  void add(const double c)
    {
      for(int i=0; i<nrows*ncols; i++)
	vals[i] += c;
    }
  void add(const CMatrix& A)
    {
      // if A is a row or column vector it is `replicated' before the operation.
      if(A.nrows==1)
	{
	  assert(A.ncols==ncols);
	  for(int i=0; i<nrows; i++)
	    for(int j=0; j<ncols; j++)
	      vals[i+nrows*j] += A.vals[j];
	}
      else if(A.ncols==1)
	{
	  assert(A.nrows==nrows);
	  for(int j=0; j<ncols; j++)
	    for(int i=0; i<nrows; i++)
	      vals[i+nrows*j] += A.vals[i];
	}
      else  
	{
	  assert(A.nrows==nrows && A.ncols == ncols);
	  for(int i=0; i<nrows*ncols; i++)
	    vals[i] += A.vals[i];
	}    
    }
  void subtract(const double c)
    {
      for(int i=0; i<nrows*ncols; i++)
	vals[i] -= c;
    }
  void subtract(const CMatrix& A)
    {
      // if A is a row or column vector it is `replicated' before the operation.
      if(A.nrows==1)
	{
	  assert(A.ncols==ncols);
	  for(int i=0; i<nrows; i++)
	    for(int j=0; j<ncols; j++)
	      vals[i+nrows*j] -= A.vals[j];
	}
      else if(A.ncols==1)
	{
	  assert(A.nrows==nrows);
	  for(int j=0; j<ncols; j++)
	    for(int i=0; i<nrows; i++)
	      vals[i+nrows*j] -= A.vals[i];
	}
      else  
	{
	  assert(A.nrows==nrows && A.ncols == ncols);
	  for(int i=0; i<nrows*ncols; i++)
	    vals[i] -= A.vals[i];
	}    
    }
  void operator+=(const double c)
    {
      add(c);
    }
  void operator+=(const CMatrix& A)
    {
      add(A);
    }
  void operator-=(const double c)
    {
      subtract(c);
    }
  void operator-=(const CMatrix& A)
    {
      subtract(A);
    }
  void operator*=(const double c)
    {
      multiply(c);
    }
  void operator*=(const CMatrix& A)
    {
      multiply(A);
    }

  void operator-()
    {
      negate();
    }
  // element by element operations
  void invElements()
    {
      for(int i = 0; i<nrows*ncols; i++)
	vals[i] = 1/vals[i];
    }
  void exp()
    {
      for(int i = 0; i<nrows*ncols; i++)
	vals[i] = std::exp(vals[i]);
    }
  void sign()
    {
      for(int i = 0; i<nrows*ncols; i++)
	{
	  if(vals[i]>0)
	    vals[i]=1;
	  else
	    vals[i]=-1;
	}
    }
  void log()
    {
      for(int i = 0; i<nrows*ncols; i++)
	vals[i] = std::log(vals[i]);
    }

  // Lapack operations
  void lu();
  void inv();
  void chol();
  // Log determinant of a positive definite matrix where U is the Cholesky decomposition of the matrix.
  double logDet(CMatrix U);

  void chol(const char* type);
  void pdinv(CMatrix U);
  void pdinv();

  // LAPACK operations
  // Cholesky factorisation.
  void potrf(const char* type);
  // inverse based on Cholesky.
  void potri(const char* type);

  // BLAS operations
  // Level 1 BLAS operations.

  // Level 2 BLAS operations.
  // y:= alpha*op(A)*x + beta*y
  void gemv(const CMatrix& A, const CMatrix& x, const double alpha, const double beta, const char* trans);

  // y(i, :)' := alpha op(A)*x(k, :)' + beta*y(i, :)';
  void gemvRowRow(const int i, const CMatrix& A, const CMatrix& x, const int k, const double alpha, const double beta, const char* trans);
  // y(i, :)' := alpha op(A)*x(:, j) + beta*y(i, :)';
  void gemvRowCol(const int i, const CMatrix& A, const CMatrix& x, const int j, const double alpha, const double beta, const char* trans);
  // y(:, j) := alpha op(A)*x(:, k) + beta*y(:, j);
  void gemvColCol(const int j, const CMatrix& A, const CMatrix& x, const int k, const double alpha, const double beta, const char* trans);
  // y(:, j) := alpha op(A)*x(i, :)' + beta*y(:, j);
  void gemvColRow(const int j, const CMatrix& A, const CMatrix& x, const int i, const double alpha, const double beta, const char* trans);


  // y:= alpha*A*x + beta*y
  void symv(const CMatrix& A, const CMatrix& x, const double alpha, const double beta, const char* upperOrLower);

  // y(i, :)' := alpha A*x(k, :)' + beta*y(i, :)';
  void symvRowRow(const int i, const CMatrix& A, const CMatrix& x, const int k, const double alpha, const double beta, const char* upperOrLower);
  // y(i, :)' := alpha A*x(:, j) + beta*y(i, :)';
  void symvRowCol(const int i, const CMatrix& A, const CMatrix& x, const int j, const double alpha, const double beta, const char* upperOrLower);
  // y(:, j) := alpha A*x(:, k) + beta*y(:, j);
  void symvColCol(const int j, const CMatrix& A, const CMatrix& x, const int k, const double alpha, const double beta, const char* upperOrLower);
  // y(:, j) := alpha A*x(i, :)' + beta*y(:, j);
  void symvColRow(const int j, const CMatrix& A, const CMatrix& x, const int i, const double alpha, const double beta, const char* upperOrLower);


  // Level 3 BLAS operations.
  // C:= alpha*op(A)*op(B) + beta*C
  void gemm(const CMatrix& A, const CMatrix& B, const double alpha, const double beta, const char* transa, const char* transb);
  //  C:=alpha*op(A)*op(A)' + beta*C.
  void syrk(const CMatrix& A, const double alpha, const double beta, const char* type, const char* trans);

  void trmm(const CMatrix& B, const double alpha, const char* side, const char* type, const char* trans, const char* diag);

  void trsm(const CMatrix& B, const double alpha, const char* side, const char* type, const char* trans, const char* diag);

  // MATLAB interaction commands
  mxArray* toMxArray() const;
  void fromMxArray(const mxArray* matlabArray);
  void fromSparseMxArray(const mxArray* matlabArray);
  void fromFullMxArray(const mxArray* matlabArray);
  // sample all elements from a Gaussian with mean and variance.
  void randn(const double var, const double mean);
  // sample all elements from a standard normal.
  void randn();
  // set all elements of the matrix to zero.
  void zeros();
  // set all elements to their negative value.
  void negate();
  // sum all elements of the matrix.
  double sum() const;
  // check if the two matrices are identical to within a tolerance.
  bool equals(const CMatrix& A, const double tol=ndlutil::MATCHTOL) const;
  // find the maximum absolute difference between matrices.
  double maxAbsDiff(const CMatrix& X) const;
  // find the maximum element of the matrix.
  double max() const; 
  // resize a matrix.
  void resize(const int rows, const int cols)
    {
      if(rows!=nrows || cols!=ncols)
	{
	  memDeAllocate();
	  nrows = rows;
	  ncols = cols;
	  memAllocate();
	}
    }
  // io operations ...
  //  ostream& operator<<(ostream& os, CMatrix& A);

  // friend functions
  friend inline void swap(CMatrix& x, CMatrix& y);

 private:
  void memAllocate();
  void memDeAllocate();
  void memReAllocate(int rowIncrease, int colIncrease);
  void copy(const CMatrix& x);
  double* vals;
  int nrows;
  int ncols;
  enum{NEW, MATLAB, MALLOC};
  static const int allocation = NEW;
  bool symmetric;
  bool triangular;
};


CMatrix lu(const CMatrix& inMatrix);
CMatrix chol(const CMatrix& inMatrix);
CMatrix inv(const CMatrix& inMatrix);
// Normal matrix multiply.
CMatrix multiply(const CMatrix& A, const CMatrix& B);
// Matrix multiply but allowing transpose operations on the matrices.
CMatrix multiply(const CMatrix& A, const char* transa, const CMatrix& B, const char* transb);
// give the trace of the matrix.
double trace(const CMatrix& A);
// sum all elements in the matrix.
double sum(const CMatrix& A);
double max(const CMatrix& A);
double dist(const CMatrix& A, const CMatrix& B);
CMatrix pdinv(const CMatrix& A);
// Overload output operator for console display.
ostream& operator<<(ostream& os, const CMatrix& A);
//ostream& operator<<(CMatrix A);

CMatrix sumRow(const CMatrix&);
CMatrix meanRow(const CMatrix&);
CMatrix sumCol(const CMatrix&);
CMatrix meanCol(const CMatrix&);

double randn(const double mean, const double var);
double randn();
#endif
