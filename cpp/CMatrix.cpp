#include <cassert>
#include "CMatrix.h"
#include "lapack.h"
using namespace std;

void CMatrix::copyRowRow(const int i, const CMatrix& X, const int k)
{
  assert(X.ncols==ncols);
  assert(i<nrows);
  assert(k<X.nrows);
  dcopy_(ncols, X.vals+k, X.nrows, vals+i, nrows);
}
void CMatrix::copy(const CMatrix& x)
{
  // Level 1 Blas operation y <- x
  assert(x.ncols == ncols);
  assert(x.nrows == nrows);
  dcopy_(ncols*nrows, x.vals, 1, vals, 1);
  symmetric = x.symmetric;
}
void CMatrix::gemv(const CMatrix& A, const CMatrix& x, const double alpha, const double beta, const char* trans)
{
  // Level 2 Blas operation y <- alpha op(A)*x + beta*y
  assert(ncols==1);
  assert(((trans[0]=='n' || trans[0]=='N') && A.ncols==x.nrows && A.nrows == nrows)
	 ||
	 ((trans[0]=='t' || trans[0]=='T') && A.nrows==x.nrows && A.ncols == nrows));
  dgemv_(trans, A.nrows, A.ncols, alpha, A.vals, A.nrows, 
	 x.vals, 1, beta, vals, 1);
}
void CMatrix::gemvRowRow(const int i, const CMatrix& A, const CMatrix& x, const int k, const double alpha, const double beta, const char* trans)
{
  // Level 2 Blas operation y(i, :)' <- alpha op(A)*x(k, :)' + beta*y(i, :)';
  assert(i<nrows);
  assert(k<x.nrows);
  assert(((trans[0]=='n' || trans[0]=='N') && A.ncols==x.ncols && A.nrows == ncols)
	 ||
	 ((trans[0]=='t' || trans[0]=='T') && A.nrows==x.ncols && A.ncols == ncols));
  dgemv_(trans, A.nrows, A.ncols, alpha, A.vals, A.nrows, 
	 x.vals+k, x.nrows, beta, vals+i, nrows);
}
void CMatrix::gemvRowCol(const int i, const CMatrix& A, const CMatrix& x, const int j, const double alpha, const double beta, const char* trans)
{
  // Level 2 Blas operation y(i, :)' <- alpha op(A)*x(:, j) + beta*y(i, :)';
  assert(i<nrows);
  assert(j<x.ncols);
  assert(((trans[0]=='n' || trans[0]=='N') && A.ncols==x.nrows && A.nrows == ncols)
	 ||
	 ((trans[0]=='t' || trans[0]=='T') && A.nrows==x.nrows && A.ncols == ncols));
  dgemv_(trans, A.nrows, A.ncols, alpha, A.vals, A.nrows, 
	 x.vals+j*x.nrows, 1, beta, vals+i, nrows);
}
void CMatrix::gemvColCol(const int j, const CMatrix& A, const CMatrix& x, const int k, const double alpha, const double beta, const char* trans)
{
  // Level 2 Blas operation y(:, j) <- alpha op(A)*x(:, k) + beta*y(:, j);
  assert(j<ncols);
  assert(k<x.ncols);
  assert(((trans[0]=='n' || trans[0]=='N') && A.ncols==x.nrows && A.nrows == nrows)
	 ||
	 ((trans[0]=='t' || trans[0]=='T') && A.nrows==x.nrows && A.ncols == nrows));
  dgemv_(trans, A.nrows, A.ncols, alpha, A.vals, A.nrows, 
	 x.vals+k*x.nrows, 1, beta, vals+j*nrows, 1);
}
void CMatrix::gemvColRow(const int j, const CMatrix& A, const CMatrix& x, const int i, const double alpha, const double beta, const char* trans)
{
  // Level 2 Blas operation y(:, j) <- alpha op(A)*x(i, :)' + beta*y(:, j);
  assert(j<ncols);
  assert(i<x.nrows);
  assert(((trans[0]=='n' || trans[0]=='N') && A.ncols==x.ncols && A.nrows == nrows)
	 ||
	 ((trans[0]=='t' || trans[0]=='T') && A.nrows==x.ncols && A.ncols == nrows));
  dgemv_(trans, A.nrows, A.ncols, alpha, A.vals, A.nrows, 
	 x.vals+i, x.nrows, beta, vals+j*nrows, 1);
}  
void CMatrix::symv(const CMatrix& A, const CMatrix& x, const double alpha, const double beta, const char* upperLower)
{
  // Level 2 Blas operation, symmetric A,  y <- alpha A*x + beta*y
  assert(A.isSymmetric());
  assert(ncols==1);
  assert(x.ncols==1);
  assert(upperLower[0]=='u' || upperLower[0]=='U' || upperLower[0]=='l' || upperLower[0]=='L');
  assert(nrows==A.nrows);
  assert(nrows==x.nrows);

  
  dsymv_(upperLower, A.ncols, alpha, A.vals, A.nrows, 
	 x.vals, 1, beta, vals, 1);
}
void CMatrix::symvRowRow(const int i, const CMatrix& A, const CMatrix& x, const int k, const double alpha, const double beta, const char* upperLower)
{
  // Level 2 Blas operation, symmetric A,  y(i, :)' <- alpha A*x(k, :)' + beta*y(i, :)';
  assert(A.isSymmetric());
  assert(i<nrows);
  assert(k<x.nrows);
  assert(upperLower[0]=='u' || upperLower[0]=='U' || upperLower[0]=='l' || upperLower[0]=='L');
  assert(ncols==A.nrows);
  assert(ncols==x.ncols);
  dsymv_(upperLower, A.ncols, alpha, A.vals, A.nrows, 
	 x.vals+k, x.nrows, beta, vals+i, nrows);
}
void CMatrix::symvRowCol(const int i, const CMatrix& A, const CMatrix& x, const int j, const double alpha, const double beta, const char* upperLower)
{
  // Level 2 Blas operation, symmetric A,  y(i, :)' <- alpha A*x(:, j) + beta*y(i, :)';
  assert(A.isSymmetric());
  assert(i<nrows);
  assert(j<x.ncols);
  assert(upperLower[0]=='u' || upperLower[0]=='U' || upperLower[0]=='l' || upperLower[0]=='L');
  assert(ncols==A.nrows);
  assert(ncols==x.nrows);
  dsymv_(upperLower, A.ncols, alpha, A.vals, A.nrows, 
	 x.vals+j*x.nrows, 1, beta, vals+i, nrows);
}
void CMatrix::symvColCol(const int j, const CMatrix& A, const CMatrix& x, const int k, const double alpha, const double beta, const char* upperLower)
{
  // Level 2 Blas operation, symmetric A,  y(:, j) <- alpha A*x(:, k) + beta*y(:, j);
  assert(A.isSymmetric());
  assert(j<ncols);
  assert(k<x.ncols);
  assert(upperLower[0]=='u' || upperLower[0]=='U' || upperLower[0]=='l' || upperLower[0]=='L');
  assert(nrows==A.nrows);
  assert(A.ncols==A.nrows);
  assert(nrows==x.nrows);
  dsymv_(upperLower, A.ncols, alpha, A.vals, A.nrows, 
	 x.vals+k*x.nrows, 1, beta, vals+j*nrows, 1);
}
void CMatrix::symvColRow(const int j, const CMatrix& A, const CMatrix& x, const int i, const double alpha, const double beta, const char* upperLower)
{
  // Level 2 Blas operation, symmetric A,  y(:, j) <- alpha A*x(i, :)' + beta*y(:, j);
  assert(A.isSymmetric());
  assert(j<ncols);
  assert(i<x.nrows);
  assert(upperLower[0]=='u' || upperLower[0]=='U' || upperLower[0]=='l' || upperLower[0]=='L');
  assert(nrows==A.nrows);
  assert(nrows==x.ncols);
  dsymv_(upperLower, A.ncols, alpha, A.vals, A.nrows, 
	 x.vals+i, x.nrows, beta, vals+j*nrows, 1);
}  
  
void CMatrix::gemm(const CMatrix& A, const CMatrix& B, const double alpha, const double beta, const char* transa, const char* transb)
{
  setSymmetric(false);
  // Level 3 Blas operation C <- alpha op(A)*op(B) + beta*C
  int m = 0;
  int n = 0;
  int k = 0;
  switch(transa[0])
    {
    case 'n':
    case 'N':
      m = A.nrows;
      k = A.ncols;
      break;
    case 't':
    case 'T':
      m=A.ncols;
      k=A.nrows;
      break;
    default:
      assert(0);
    }
  switch(transb[0])
    {
    case 'n':
    case 'N':
      n = B.ncols;
      assert(k==B.nrows);
      break;
    case 't':
    case 'T':
      n = B.nrows;
      assert(k==B.ncols);
      break;
    default:
      assert(0);
    }	
  assert(n==ncols);
  assert(m==nrows);
  dgemm_(transa, transb, m, n, k, alpha, A.vals, A.nrows, 
	 B.vals, B.nrows, beta, vals, nrows);
  
}
void CMatrix::trmm(const CMatrix& A, const double alpha, const char* side, const char* type, const char* trans, const char* diag)
{
  assert(side[0]=='L' || side[0]=='l' || side[0]=='R' || side[0]=='r');
  assert(type[0]=='L' || type[0]=='l' || type[0]=='U' || type[0]=='u');
  assert(trans[0]=='N' || trans[0]=='n' || trans[0]=='T' || trans[0]=='t');
  assert(diag[0]=='N' || diag[0]=='n' || diag[0]=='U' || diag[0]=='u');
  
  assert(A.isTriangular());
  switch(side[0])
    {
    case 'L':
    case 'l':
      assert(A.nrows==nrows);
      break;
    case 'R':
    case 'r':
      assert(A.nrows==ncols);
      break;
    otherwise:
      cerr << "No such value for side.";
    }
  
  dtrmm_(side, type, trans, diag, nrows, ncols, alpha, A.vals, A.nrows, vals, nrows);
}
void CMatrix::trsm(const CMatrix& A, const double alpha, const char* side, const char* type, const char* trans, const char* diag)
{
  assert(side[0]=='L' || side[0]=='l' || side[0]=='R' || side[0]=='r');
  assert(type[0]=='L' || type[0]=='l' || type[0]=='U' || type[0]=='u');
  assert(trans[0]=='N' || trans[0]=='n' || trans[0]=='T' || trans[0]=='t');
  assert(diag[0]=='N' || diag[0]=='n' || diag[0]=='U' || diag[0]=='u');
  
  assert(A.isTriangular());
  switch(side[0])
    {
    case 'L':
    case 'l':
      assert(A.nrows==nrows);
      break;
    case 'R':
    case 'r':
      assert(A.nrows==ncols);
      break;
    otherwise:
      cerr << "No such value for side.";
    }
  
  dtrsm_(side, type, trans, diag, nrows, ncols, alpha, A.vals, A.nrows, vals, nrows);
}

void CMatrix::syrk(const CMatrix& A, const double alpha, const double beta, const char* type, const char* trans)
{
  assert(isSymmetric() || beta==0.0);  
  int n = 0;
  int k = 0;
  switch(trans[0])
    {
    case 'n':
    case 'N':
      n = ncols;
      k = A.ncols;
      assert(n==A.nrows);
      break;
    case 't':
    case 'T':
      n = nrows;
      k = A.nrows;
      assert(n==A.ncols);
      break;
    default:
      assert(0);
    }
      
  dsyrk_(type, trans, n, k, alpha, A.vals, A.nrows, beta, vals, nrows);
  copySymmetric(type);
}
void CMatrix::copySymmetric(const char* type)
{
  switch(type[0]) 
    {
    case 'U':
    case 'u':
      for(int i=0; i<nrows; i++)
	for(int j=0; j<i; j++)
	  vals[i + nrows*j] = vals[j + ncols*i];
      break;
    case 'L':
    case 'l':
      for(int j=0; j<ncols; j++)
	for(int i=0; i<j; i++)
	  vals[i + nrows*j] = vals[j + ncols*i];
      break;
    default:
      assert(0);
    }
}

void CMatrix::potrf(const char* type)
{
  assert(isSymmetric());
  int info;
  dpotrf_(type, nrows, vals, ncols, info);
  setSymmetric(false);
  setTriangular(true);
  if(info!=0) throw "Matrix not positive definite.";
}
void CMatrix::chol(const char* type)
{
  // type is either U or L for upper or lower triangular.
  potrf(type);
  switch((int)type[0]) {
  case 'L':
      for(int j=0; j<ncols; j++)
	for(int i=0; i<j; i++)
	  vals[i + nrows*j] = 0.0;
    break;
  case 'U':
      for(int i=0; i<nrows; i++)
	for(int j=0; j<i; j++)
	  vals[i + nrows*j] = 0.0;
    break;
  }
      
}
void CMatrix::chol()
{
  chol("U");
}
double CMatrix::logDet(CMatrix U)
{
  assert(isSymmetric()); /// actually should be positive definite.
  double logDet = 0.0;
  for(int i=0; i<U.getRows(); i++)
    logDet+=std::log(U.getVal(i, i));
  logDet *= 2;
  return logDet;
}

void CMatrix::potri(const char* type)
{
  assert(isSquare());
  int info;
  dpotri_(type, nrows, vals, ncols, info);
  if(info!=0) throw "Matrix not positive definite.";
}
void CMatrix::pdinv(CMatrix U)
{
  deepCopy(U);
  potri("U");
  for(int i=0; i<nrows; i++)
    for(int j=0; j<i; j++)
      vals[i + nrows*j] = vals[j + ncols*i];
}
void CMatrix::pdinv()
{
  potrf("U");
  potri("U");
  for(int i=0; i<nrows; i++)
    for(int j=0; j<i; j++)
      vals[i + nrows*j] = vals[j + ncols*i];
}

void CMatrix::lu()
{
  // this isn't really properly implemented yet ... need to return ipiv somehow.
  assert(isSquare());
  int info;
  int* ipiv = new int[nrows];
  // TODO should really check for errors here.
  dgetrf_(nrows, ncols, vals, ncols, ipiv, info);
  if(info!=0) throw "Matrix not full rank.";
  delete[] ipiv;
}

void CMatrix::inv()
{
  assert(isSquare());
  // create output matrix by lu decomposition of input
    
  int length = nrows;
  int* ipiv = new int[length];
  int info = 0;
    
  dgetrf_(nrows, ncols, 
	  vals, ncols, ipiv, info);
  if(info!=0) throw "Matrix not full rank.";
  int order = nrows;
  int lwork = order*16;
  double* work = new double[lwork];
  info = 0;
  dgetri_(order, vals, ncols, 
	  ipiv, work, lwork, info);
  // check for successful inverse
  if(info!=0) throw "Matrix not full rank.";
  delete[] work;
  delete[] ipiv;
}
mxArray* CMatrix::toMxArray() const
{
  int dims[2];
  dims[0] = nrows;
  dims[1] = ncols;
  mxArray* matlabArray = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  double* matlabVals = mxGetPr(matlabArray);
  dcopy_(ncols*nrows, vals, 1, matlabVals, 1);
  return matlabArray;
}
void CMatrix::fromMxArray(const mxArray* matlabArray) 
{
  // TODO implement true sparse matrices to handle this.
  if(mxIsSparse(matlabArray))
    {
      fromSparseMxArray(matlabArray);
    }
  else
    {
      fromFullMxArray(matlabArray);
    }
}
void CMatrix::fromFullMxArray(const mxArray* matlabArray)
{
  assert(!mxIsSparse(matlabArray));
  if(mxGetClassID(matlabArray) != mxDOUBLE_CLASS)
    {
      cerr << "mxArray is not a double matrix." << endl;
    }
  if(mxGetNumberOfDimensions(matlabArray) != 2)
    {
      cerr << "mxArray does not have 2 dimensions." << endl;
    }
  const int* dims = mxGetDimensions(matlabArray);
  resize(dims[0], dims[1]);
  double* matlabVals = mxGetPr(matlabArray);
  dcopy_(ncols*nrows, matlabVals, 1, vals, 1);
}
void CMatrix::fromSparseMxArray(const mxArray* matlabArray)
{
  assert(mxIsSparse(matlabArray));
  if(mxGetClassID(matlabArray) != mxDOUBLE_CLASS)
    {
      cerr << "mxArray is not a double matrix." << endl;
    }
  if(mxGetNumberOfDimensions(matlabArray) != 2)
    {
      cerr << "mxArray does not have 2 dimensions." << endl;
    }
  const int* dims = mxGetDimensions(matlabArray);
  resize(dims[0], dims[1]);
  double* matlabVals = mxGetPr(matlabArray);
  setVals(0.0);
  int* matlabIr = mxGetIr(matlabArray);
  int* matlabJc = mxGetJc(matlabArray);
  int nnz = matlabJc[getCols()];
  for(int j=0; j<getCols(); j++)
    for(int i=matlabJc[j]; i<matlabJc[j+1]; i++)
      {
	setVal(matlabVals[i], matlabIr[i], j);
      } 
}    
void CMatrix::randn(const double var, const double mean)
{
  dgrand_(0);
  double sd = sqrt(var);
  double val;
  for(int i=0; i<nrows*ncols; i++)
    {
      val = dgrand_(0);
      vals[i]=val*sd+mean;
    }
}
double CMatrix::sum() const
{
  // matrix should be square
  double sum = 0.0;
  for(int i=0; i<nrows*ncols; i++)
    sum += vals[i];
  return sum;
  
}
bool CMatrix::equals(const CMatrix& A, const double tol) const
{
  if(nrows != A.nrows)
    return false;
  if(ncols != A.ncols)
    return false;
  if(maxAbsDiff(A)>tol)
    return false;
  return true;
}
double CMatrix::maxAbsDiff(const CMatrix& X) const
{
  assert(dimensionsMatch(X));
  double max=0.0;
  double diff=0.0;
  for(int i=0; i<nrows*ncols; i++)
    {
      diff = abs(vals[i]-X.vals[i]);
      if(diff>max)
	max = diff;
    }
  return max;
} 

double CMatrix::max() const
{
  // matrix should be square
  double max = vals[0];
  double val = 0.0;
  for(int i=1; i<nrows*ncols; i++)
    val = vals[i];
    if(val > max)
      max = val;
  return max;
}
void CMatrix::randn()
{
  randn(1.0, 0.0);
}
void CMatrix::negate()
{
  scale(-1.0);
}
void CMatrix::zeros()
{
  for(int i=0; i<nrows*ncols; i++)
    {
      vals[i]=0.0;
    }
}
void CMatrix::memAllocate()
{
  switch(allocation)
    {
    case NEW:
      vals = new double[nrows*ncols];
      break;
    case MALLOC:
      vals = (double *)malloc(sizeof(double)*nrows*ncols);
      break;
    default:
      cerr << "Memory allocation method not known.";
    }
}
void CMatrix::memDeAllocate()
{
  switch(allocation)
    {
    case NEW:
      delete[] vals;
      break;
    case MALLOC:
      free(vals);
      break;
    case MATLAB:
      // do nothing MATLAB should handle it.
      break;
    default:
      cerr << "Error allocation method not known";
    }
}
void CMatrix::memReAllocate(int rowIncrease, int colIncrease)
{
  int newRows = nrows+rowIncrease;
  int newCols = ncols+colIncrease;
  int minRows = newRows;
  int minCols = newCols;
  double* newVals;
  
  if(nrows<minRows) minRows = nrows;
  if(ncols<minCols) minCols = ncols;
  assert(newRows>0);
  assert(newCols>0);
  switch(allocation)
    {
    case NEW:
      newVals = new double[newRows*newCols];
      for(int i=0; i<minRows; i++)
	for(int j=0; j<minCols; j++)
	  newVals[i + newRows*j] = vals[i + nrows*j];
      // set remaining elements to zero.
      for(int i=minRows; i<newRows; i++)
	for(int j=0; j<newCols; j++)
	  newVals[i+newRows*j] = 0.0;
      for(int j=minCols; j<newCols; j++)
	for(int i=0; i<newRows; i++)
	  newVals[i+newRows*j] = 0.0;
      delete []vals;
      vals = newVals;
      nrows=newRows;
      ncols=newCols;
      break;
    case MALLOC:
      break;
    case MATLAB:
      break;
    default:
      cerr << "Memory allocation method not known";
    }
  
}
ostream& operator<<(ostream& out, const CMatrix& A)
{
  for(int i = 0; i < A.getRows(); i++){
    for(int j = 0; j < A.getCols(); j++){
      if (A.getVal(i, j) > ndlutil::DISPEPS || A.getVal(i, j) < -ndlutil::DISPEPS)
        out << A.getVal(i, j) << " ";
      else
        out << 0 << " ";
    }
    out << endl;
  }
  return out;
}

CMatrix lu(const CMatrix& inMatrix) 
{
  
  CMatrix outMatrix(inMatrix);
  outMatrix.lu();
  return outMatrix;
}

CMatrix chol(const CMatrix& inMatrix) 
{
  
  CMatrix outMatrix(inMatrix);
  outMatrix.chol("U");
  return outMatrix;
}

CMatrix inv(const CMatrix& inMatrix)
{
  CMatrix outMatrix(inMatrix);
  outMatrix.inv();
  return outMatrix;
}

CMatrix multiply(const CMatrix& A, const CMatrix& B)
{
  CMatrix C(A.getRows(), B.getCols());
  C.gemm(A, B, 1.0, 0.0, "n", "n");
  return C;
}

CMatrix multiply(const CMatrix& A, const char* transa, const CMatrix& B, const char* transb)
{
  int n=0;
  switch(transa[0])
    {
    case 'n':
    case 'N':
      n = A.getRows();
      break;
    case 't':
    case 'T':
      n = A.getCols();
      break;
    default:
      assert(0);
    }
  int m = 0;
  switch(transb[0])
    {
    case 'n':
    case 'N':
      m = B.getCols();
      break;
    case 't':
    case 'T':
      m = B.getRows();
      break;
    default:
      assert(0);
    }
  CMatrix C(n, m);
  C.gemm(A, B, 1.0, 0.0, transa, transb);
  return C;
}
double trace(const CMatrix& A)
{
  // matrix should be square
  assert(A.isSquare());
  double tr = 0;
  for(int i=0; i<A.getRows(); i++)
    tr += A.getVal(i, i);
  return tr;
}
double sum(const CMatrix& A)
{
  return A.sum();
}
double max(const CMatrix& A)
{
  // matrix should be square
  return A.max();
}
CMatrix pdinv(const CMatrix& A)
{
  CMatrix B(A);
  B.pdinv();
  return B;
}

CMatrix sumRow(const CMatrix& A)
{
  double* s = new double[A.getCols()];
  for(int j=0; j<A.getCols(); j++)
    {
      s[j] = 0.0;
      for(int i=0; i<A.getRows(); i++)
	{
	  s[j] += A.getVal(i, j);
	}
    }
  CMatrix S(1, A.getCols(), s);
  return S;
} 
CMatrix meanRow(const CMatrix& A)
{
  CMatrix M = sumRow(A);
  for(int j=0; j<A.getCols(); j++)
    M.setVal(M.getVal(0, j)/A.getRows(), 0, j);
  return M;
}
CMatrix sumCol(const CMatrix& A)
{
  double* s = new double[A.getRows()];
  for(int i=0; i<A.getRows(); i++)
    {
      s[i] = 0.0;
      for(int j=0; j<A.getCols(); j++)
	{
	  s[i] += A.getVal(i, j);
	}
    }
  CMatrix S(A.getRows(), 1, s);
  return S;
} 
CMatrix meanCol(const CMatrix& A)
{
  CMatrix M = sumCol(A);
  for(int i=0; i<A.getRows(); i++)
    M.setVal(M.getVal(i, 0)/A.getCols(), i, 0);
  return M;
}

double randn(const double mean, const double var)
{
  dgrand_(0);
  double sd = sqrt(var);
  double val = dgrand_(0);
  val*=sd;
  val+=mean;
  return val;
}
double randn()
{
  return randn(0.0, 1.0);
}
inline void swap(CMatrix& x, CMatrix& y)
{
  // Level 1 Blas operation y <-> x
  assert(y.ncols==x.ncols);
  assert(y.nrows==x.nrows);   
  dswap_(x.ncols*x.nrows, x.vals, 1, y.vals, 1);
}
