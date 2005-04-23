#include <iostream>
#include "CMatrix.h"
using namespace std;
const double MAXDIFF=1e-13;

int testDot();
int testNorm();
int testInv();
int testCholesky();
int testRandn();
int testGemm();
int testSyrk();
int testTrmm();
int testTrsm();
int testGemv();
int testGer();
int testSyr();
int testSymv();
int testAxpy();
int testScale();

int main()
{
  int fail = 0;
  fail += testInv();
  fail += testCholesky();
  fail += testRandn();

  // Level 3 Blas
  fail += testGemm();
  fail += testSyrk();
  fail += testTrmm();
  fail += testTrsm();
  // Level 2 Blas
  fail += testGemv();
  fail += testSymv();
  fail += testGer();
  fail += testSyr();

  // Level 1 Blas
  fail += testAxpy();
  fail += testScale();
  fail += testDot();
  fail += testNorm();


  cout << endl << "Total number of failures: " << fail << endl;
}

int testDot()
{
  cout << "Test dot not yet implmenented" << endl;
  return 0;
}
int testNorm()
{
  cout << "Test norm not yet implemented" << endl;
  return 0;
}
int testInv()
{
  // implicitly tests dgetrf.
  int fail = 0;
  // Matrix inverse test
  CMatrix A;
  A.readMatlabFile("testInv.mat", "A");
  CMatrix Ainv;
  Ainv.readMatlabFile("testInv.mat", "Ainv");
  A.inv();
  if(Ainv.equals(A))
      cout << "Matrix inverse matches." << endl;
  else
    {
      cout << "FAILURE: Matrix inverse." << endl;
      fail++;
    }
  return fail;
}
int testCholesky()
{
  // implicitly tests dpotrf
  int fail = 0;
  // Lower Cholesky test.
  CMatrix C;
  C.readMatlabFile("testCholesky.mat", "C");
  C.setSymmetric(true);
  CMatrix L;
  L.readMatlabFile("testCholesky.mat", "L");
  C.chol("L");
  if(L.equals(C))
    {
      cout << "Lower Cholesky factor matches."<< endl;
    }
  else
    {
      cout << "FAILURE: Lower Cholesky" << endl;
      fail++;
    }

  // Upper Cholesky test.
  C.trans();
  CMatrix U;
  U.readMatlabFile("testCholesky.mat", "U");
  if(U.equals(C))
    cout << "Upper Cholesky factor matches."<< endl;
  else
    cout << "FAILURE: Upper Cholesky." << endl;
  return fail;
}
int testRandn()
{
  int fail = 0;
  int nrows = 1000;
  int ncols = 10;
  // Random number test.
  CMatrix normRand(nrows, ncols);
  normRand.randn();
  CMatrix mean = meanCol(meanRow(normRand));
  if(abs(mean.getVal(0))<1e-2)
    cout << "randn mean matches." << endl;
  else
    {
      cout << "POSSIBLE FAILURE: randn mean " << mean.getVal(0) << "." << endl;
      fail++;
    }
  normRand *= normRand;
  CMatrix var = meanCol(meanRow(normRand));
  if(abs(var.getVal(0)-1)<1e-2)
    cout << "randn variance matches." << endl;
  else
    {
      cout << "POSSIBLE FAILURE: randn variance " << var.getVal(0) << "." << endl;
      fail++;
    }
  
  // TODO Check variance as well ... not sure of best way of checking randn!
  return fail;
} 
int testGemm()
{
  int fail = 0;
  // gemm test
  CMatrix D;
  D.readMatlabFile("testGemm.mat", "D");
  CMatrix E;
  E.readMatlabFile("testGemm.mat", "E");
  CMatrix F;
  F.readMatlabFile("testGemm.mat", "F");
  CMatrix G;
  G.readMatlabFile("testGemm.mat", "G");
  CMatrix H;
  H.readMatlabFile("testGemm.mat", "H");
  CMatrix alph;
  alph.readMatlabFile("testGemm.mat", "alpha");
  double alpha = alph.getVal(0);
  CMatrix bet;
  bet.readMatlabFile("testGemm.mat", "beta");
  double beta = bet.getVal(0);
  
  // "n" "n" test
  CMatrix GEMM1;
  GEMM1.readMatlabFile("testGemm.mat", "GEMM1");
  F.gemm(D, E, alpha, beta, "n", "N");
  if(F.equals(GEMM1))
    cout << "gemm nn matches." << endl;
  else
    {
      cout << "FAILURE: gemm nn." << endl;
      fail++;
    }

  // "t" "t" test
  CMatrix GEMM2;
  GEMM2.readMatlabFile("testGemm.mat", "GEMM2");
  G.gemm(D, E, alpha, beta, "t", "T");
  if(G.equals(GEMM2))
    cout << "gemm tt matches." << endl;
  else
    {
      cout << "FAILURE: gemm tt." << endl;
      fail++;
    }
  CMatrix GEMM3;
  GEMM3.readMatlabFile("testGemm.mat", "GEMM3");
  GEMM1.gemm(D, H, alpha, beta, "N", "t");
  if(GEMM1.equals(GEMM3))
    cout << "gemm nt matches." << endl;
  else
    {
      cout << "FAILURE: gemm tt." << endl;
      fail++;
    }
  CMatrix GEMM4;
  GEMM4.readMatlabFile("testGemm.mat", "GEMM4");
  GEMM2.gemm(D, H, alpha, beta, "T", "n");
  if(GEMM2.equals(GEMM4))
    cout << "gemm tn matches." << endl;
  else
    {
      cout << "FAILURE: gemm tt." << endl;
      fail++;
    }
  return fail;
}
int testSyrk()
{
  int fail = 0;
  // syrk test
  CMatrix A;
  A.readMatlabFile("testSyrk.mat", "A");
  CMatrix C;
  C.readMatlabFile("testSyrk.mat", "C");
  C.setSymmetric(true);
  CMatrix D;
  D.readMatlabFile("testSyrk.mat", "D");
  D.setSymmetric(true);
  CMatrix alph;
  alph.readMatlabFile("testSyrk.mat", "alpha");
  double alpha = alph.getVal(0);
  CMatrix bet;
  bet.readMatlabFile("testSyrk.mat", "beta");
  double beta = bet.getVal(0);
  CMatrix SYRK1;
  SYRK1.readMatlabFile("testSyrk.mat", "SYRK1");
  CMatrix SYRK2;
  SYRK2.readMatlabFile("testSyrk.mat", "SYRK2");
  CMatrix F;

  F.deepCopy(C);
  F.syrk(A, alpha, beta, "u", "N");
  if(F.equals(SYRK1))
    cout << "syrk un matches." << endl;
  else
    {
      cout << "FAILURE: syrk un." << endl;
      fail++;
    }
  F.deepCopy(C);
  F.syrk(A, alpha, beta, "L", "n");
  if(F.equals(SYRK1))
    cout << "syrk ln matches." << endl;
  else
    {
      cout << "FAILURE: syrk ln." << endl;
      fail++;
    }
  F.deepCopy(D);
  F.syrk(A, alpha, beta, "U", "t");
  if(F.equals(SYRK2))
    cout << "syrk ut matches." << endl;
  else
    {
      cout << "FAILURE: syrk ut." << endl;
      fail++;
    }
  F.deepCopy(D);
  F.syrk(A, alpha, beta, "l", "T");
  if(F.equals(SYRK2))
    cout << "syrk lt matches." << endl;
  else
    {
      cout << "FAILURE: syrk lt." << endl;
      fail++;
    }
  return fail;
}
int testTrmm()
{
  int fail = 0;
  CMatrix B;
  B.readMatlabFile("testTrmm.mat", "B");
  CMatrix alph;
  alph.readMatlabFile("testTrmm.mat", "alpha");
  double alpha = alph.getVal(0);
  CMatrix L;
  L.readMatlabFile("testTrmm.mat", "L");
  L.setTriangular(true);
  CMatrix L2;
  L2.readMatlabFile("testTrmm.mat", "L2");
  L2.setTriangular(true);
  CMatrix U;
  U.readMatlabFile("testTrmm.mat", "U");
  U.setTriangular(true);
  CMatrix U2;
  U2.readMatlabFile("testTrmm.mat", "U2");
  U2.setTriangular(true);
  CMatrix F;

  CMatrix TRMM1;
  TRMM1.readMatlabFile("testTrmm.mat", "TRMM1");
  F.deepCopy(B);
  F.trmm(L, alpha, "L", "L", "N", "N");
  if(F.equals(TRMM1))
    cout << "trmm llnn matches." << endl;
  else
    {
      cout << "FAILURE: trmm llnn." << endl;
      fail++;
    }

  CMatrix TRMM2;
  TRMM2.readMatlabFile("testTrmm.mat", "TRMM2");
  F.deepCopy(B);
  F.trmm(L, alpha, "L", "L", "T", "N");
  if(F.equals(TRMM2))
    cout << "trmm lltn matches." << endl;
  else
    {
      cout << "FAILURE: trmm lltn." << endl;
      fail++;
    }
  
  CMatrix TRMM3;
  TRMM3.readMatlabFile("testTrmm.mat", "TRMM3");
  F.deepCopy(B);
  F.trmm(L2, alpha, "R", "L", "N", "N");
  if(F.equals(TRMM3))
    cout << "trmm rlnn matches." << endl;
  else
    {
      cout << "FAILURE: trmm rlnn." << endl;
      fail++;
    }
  CMatrix TRMM4;
  TRMM4.readMatlabFile("testTrmm.mat", "TRMM4");
  F.deepCopy(B);
  F.trmm(L2, alpha, "R", "L", "T", "N");
  if(F.equals(TRMM4))
    cout << "trmm rltn matches." << endl;
  else
    {
      cout << "FAILURE: trmm rltn." << endl;
      fail++;
    }

  CMatrix TRMM5;
  TRMM5.readMatlabFile("testTrmm.mat", "TRMM5");
  F.deepCopy(B);
  F.trmm(L, alpha, "L", "L", "N", "U");
  if(F.equals(TRMM5))
    cout << "trmm llnu matches." << endl;
  else
    {
      cout << "FAILURE: trmm llnu." << endl;
      fail++;
    }

  CMatrix TRMM6;
  TRMM6.readMatlabFile("testTrmm.mat", "TRMM6");
  F.deepCopy(B);
  F.trmm(L, alpha, "L", "L", "T", "U");
  if(F.equals(TRMM6))
    cout << "trmm lltu matches." << endl;
  else
    {
      cout << "FAILURE: trmm lltu." << endl;
      fail++;
    }
  
  CMatrix TRMM7;
  TRMM7.readMatlabFile("testTrmm.mat", "TRMM7");
  F.deepCopy(B);
  F.trmm(L2, alpha, "R", "L", "N", "U");
  if(F.equals(TRMM7))
    cout << "trmm rlnu matches." << endl;
  else
    {
      cout << "FAILURE: trmm rlnu." << endl;
      fail++;
    }
  CMatrix TRMM8;
  TRMM8.readMatlabFile("testTrmm.mat", "TRMM8");
  F.deepCopy(B);
  F.trmm(L2, alpha, "R", "L", "T", "U");
  if(F.equals(TRMM8))
    cout << "trmm rltu matches." << endl;
  else
    {
      cout << "FAILURE: trmm rltu." << endl;
      fail++;
    }
  CMatrix TRMM9;
  TRMM9.readMatlabFile("testTrmm.mat", "TRMM9");
  F.deepCopy(B);
  F.trmm(U, alpha, "L", "U", "N", "N");
  if(F.equals(TRMM9))
    cout << "trmm lunn matches." << endl;
  else
    {
      cout << "FAILURE: trmm lunn." << endl;
      fail++;
    }

  CMatrix TRMM10;
  TRMM10.readMatlabFile("testTrmm.mat", "TRMM10");
  F.deepCopy(B);
  F.trmm(U, alpha, "L", "U", "T", "N");
  if(F.equals(TRMM10))
    cout << "trmm lutn matches." << endl;
  else
    {
      cout << "FAILURE: trmm lutn." << endl;
      fail++;
    }
  
  CMatrix TRMM11;
  TRMM11.readMatlabFile("testTrmm.mat", "TRMM11");
  F.deepCopy(B);
  F.trmm(U2, alpha, "R", "U", "N", "N");
  if(F.equals(TRMM11))
    cout << "trmm runn matches." << endl;
  else
    {
      cout << "FAILURE: trmm runn." << endl;
      fail++;
    }
  CMatrix TRMM12;
  TRMM12.readMatlabFile("testTrmm.mat", "TRMM12");
  F.deepCopy(B);
  F.trmm(U2, alpha, "R", "U", "T", "N");
  if(F.equals(TRMM12))
    cout << "trmm rutn matches." << endl;
  else
    {
      cout << "FAILURE: trmm rutn." << endl;
      fail++;
    }

  CMatrix TRMM13;
  TRMM13.readMatlabFile("testTrmm.mat", "TRMM13");
  F.deepCopy(B);
  F.trmm(U, alpha, "L", "U", "N", "U");
  if(F.equals(TRMM13))
    cout << "trmm lunu matches." << endl;
  else
    {
      cout << "FAILURE: trmm lunu." << endl;
      fail++;
    }

  CMatrix TRMM14;
  TRMM14.readMatlabFile("testTrmm.mat", "TRMM14");
  F.deepCopy(B);
  F.trmm(U, alpha, "L", "U", "T", "U");
  if(F.equals(TRMM14))
    cout << "trmm lutu matches." << endl;
  else
    {
      cout << "FAILURE: trmm lutu." << endl;
      fail++;
    }
  
  CMatrix TRMM15;
  TRMM15.readMatlabFile("testTrmm.mat", "TRMM15");
  F.deepCopy(B);
  F.trmm(U2, alpha, "R", "U", "N", "U");
  if(F.equals(TRMM15))
    cout << "trmm runu matches." << endl;
  else
    {
      cout << "FAILURE: trmm runu." << endl;
      fail++;
    }
  CMatrix TRMM16;
  TRMM16.readMatlabFile("testTrmm.mat", "TRMM16");
  F.deepCopy(B);
  F.trmm(U2, alpha, "R", "U", "T", "U");
  if(F.equals(TRMM16))
    cout << "trmm rutu matches." << endl;
  else
    {
      cout << "FAILURE: trmm rutu." << endl;
      fail++;
    }

  
  return fail;
}
int testTrsm()
{
  double tolInv = 1e-8;
  int fail = 0;
  CMatrix B;
  B.readMatlabFile("testTrsm.mat", "B");
  CMatrix alph;
  alph.readMatlabFile("testTrsm.mat", "alpha");
  double alpha = alph.getVal(0);
  CMatrix L;
  L.readMatlabFile("testTrsm.mat", "L");
  L.setTriangular(true);
  CMatrix L2;
  L2.readMatlabFile("testTrsm.mat", "L2");
  L2.setTriangular(true);
  CMatrix U;
  U.readMatlabFile("testTrsm.mat", "U");
  U.setTriangular(true);
  CMatrix U2;
  U2.readMatlabFile("testTrsm.mat", "U2");
  U2.setTriangular(true);
  CMatrix F;

  CMatrix TRSM1;
  TRSM1.readMatlabFile("testTrsm.mat", "TRSM1");
  F.deepCopy(B);
  F.trsm(L, alpha, "L", "L", "N", "N");
  double absDiff=F.maxAbsDiff(TRSM1);
  if(absDiff<tolInv)
    cout << "trsm llnn matches." << endl;
  else
    {
      cout << "FAILURE: trsm llnn, absolute difference " << absDiff << "." << endl;
      fail++;
    }

  CMatrix TRSM2;
  TRSM2.readMatlabFile("testTrsm.mat", "TRSM2");
  F.deepCopy(B);
  F.trsm(L, alpha, "L", "L", "T", "N");
  absDiff=F.maxAbsDiff(TRSM2);
  if(absDiff<tolInv)
    cout << "trsm lltn matches." << endl;
  else
    {
      cout << "FAILURE: trsm lltn, absolute difference " << absDiff << "." << endl;
      fail++;
    }
  
  CMatrix TRSM3;
  TRSM3.readMatlabFile("testTrsm.mat", "TRSM3");
  F.deepCopy(B);
  F.trsm(L2, alpha, "R", "L", "N", "N");
  absDiff=F.maxAbsDiff(TRSM3);
  if(absDiff<tolInv)
    cout << "trsm rlnn matches." << endl;
  else
    {
      cout << "FAILURE: trsm rlnn, absolute difference " << absDiff << "." << endl;
      fail++;
    }
  CMatrix TRSM4;
  TRSM4.readMatlabFile("testTrsm.mat", "TRSM4");
  F.deepCopy(B);
  F.trsm(L2, alpha, "R", "L", "T", "N");
  absDiff=F.maxAbsDiff(TRSM4);
  if(absDiff<tolInv)
    cout << "trsm rltn matches." << endl;
  else
    {
      cout << "FAILURE: trsm rltn, absolute difference " << absDiff << "." << endl;
      fail++;
    }

  CMatrix TRSM5;
  TRSM5.readMatlabFile("testTrsm.mat", "TRSM5");
  F.deepCopy(B);
  F.trsm(L, alpha, "L", "L", "N", "U");
  absDiff=F.maxAbsDiff(TRSM5);
  if(absDiff<tolInv)
    cout << "trsm llnu matches." << endl;
  else
    {
      cout << "FAILURE: trsm llnu, absolute difference " << absDiff << "." << endl;
      fail++;
    }

  CMatrix TRSM6;
  TRSM6.readMatlabFile("testTrsm.mat", "TRSM6");
  F.deepCopy(B);
  F.trsm(L, alpha, "L", "L", "T", "U");
  absDiff=F.maxAbsDiff(TRSM6);
  if(absDiff<tolInv)
    cout << "trsm lltu matches." << endl;
  else
    {
      cout << "FAILURE: trsm lltu, absolute difference " << absDiff << "." << endl;
      fail++;
    }
  
  CMatrix TRSM7;
  TRSM7.readMatlabFile("testTrsm.mat", "TRSM7");
  F.deepCopy(B);
  F.trsm(L2, alpha, "R", "L", "N", "U");
  absDiff=F.maxAbsDiff(TRSM7);
  if(absDiff<tolInv)
    cout << "trsm rlnu matches." << endl;
  else
    {
      cout << "FAILURE: trsm rlnu, absolute difference " << absDiff << "." << endl;
      fail++;
    }
  CMatrix TRSM8;
  TRSM8.readMatlabFile("testTrsm.mat", "TRSM8");
  F.deepCopy(B);
  F.trsm(L2, alpha, "R", "L", "T", "U");
  absDiff=F.maxAbsDiff(TRSM8);
  if(absDiff<tolInv)
    cout << "trsm rltu matches." << endl;
  else
    {
      cout << "FAILURE: trsm rltu, absolute difference " << absDiff << "." << endl;
      fail++;
    }
  CMatrix TRSM9;
  TRSM9.readMatlabFile("testTrsm.mat", "TRSM9");
  F.deepCopy(B);
  F.trsm(U, alpha, "L", "U", "N", "N");
  absDiff=F.maxAbsDiff(TRSM9);
  if(absDiff<tolInv)
    cout << "trsm lunn matches." << endl;
  else
    {
      cout << "FAILURE: trsm lunn, absolute difference " << absDiff << "." << endl;
      fail++;
    }

  CMatrix TRSM10;
  TRSM10.readMatlabFile("testTrsm.mat", "TRSM10");
  F.deepCopy(B);
  F.trsm(U, alpha, "L", "U", "T", "N");
  absDiff=F.maxAbsDiff(TRSM10);
  if(absDiff<tolInv)
    cout << "trsm lutn matches." << endl;
  else
    {
      cout << "FAILURE: trsm lutn, absolute difference " << absDiff << "." << endl;
      fail++;
    }
  
  CMatrix TRSM11;
  TRSM11.readMatlabFile("testTrsm.mat", "TRSM11");
  F.deepCopy(B);
  F.trsm(U2, alpha, "R", "U", "N", "N");
  absDiff=F.maxAbsDiff(TRSM11);
  if(absDiff<tolInv)
    cout << "trsm runn matches." << endl;
  else
    {
      cout << "FAILURE: trsm runn, absolute difference " << absDiff << "." << endl;
      fail++;
    }
  CMatrix TRSM12;
  TRSM12.readMatlabFile("testTrsm.mat", "TRSM12");
  F.deepCopy(B);
  F.trsm(U2, alpha, "R", "U", "T", "N");
  absDiff=F.maxAbsDiff(TRSM12);
  if(absDiff<tolInv)
    cout << "trsm rutn matches." << endl;
  else
    {
      cout << "FAILURE: trsm rutn, absolute difference " << absDiff << "." << endl;
      fail++;
    }

  CMatrix TRSM13;
  TRSM13.readMatlabFile("testTrsm.mat", "TRSM13");
  F.deepCopy(B);
  F.trsm(U, alpha, "L", "U", "N", "U");
  absDiff=F.maxAbsDiff(TRSM13);
  if(absDiff<tolInv)
    cout << "trsm lunu matches." << endl;
  else
    {
      cout << "FAILURE: trsm lunu, absolute difference " << absDiff << "." << endl;
      fail++;
    }

  CMatrix TRSM14;
  TRSM14.readMatlabFile("testTrsm.mat", "TRSM14");
  F.deepCopy(B);
  F.trsm(U, alpha, "L", "U", "T", "U");
  absDiff=F.maxAbsDiff(TRSM14);
  if(absDiff<tolInv)
    cout << "trsm lutu matches." << endl;
  else
    {
      cout << "FAILURE: trsm lutu, absolute difference " << absDiff << "." << endl;
      fail++;
    }
  
  CMatrix TRSM15;
  TRSM15.readMatlabFile("testTrsm.mat", "TRSM15");
  F.deepCopy(B);
  F.trsm(U2, alpha, "R", "U", "N", "U");
  absDiff = F.maxAbsDiff(TRSM15);
  if(absDiff<tolInv)
    cout << "trsm runu matches." << endl;
  else
    {
      cout << "FAILURE: trsm runu, absolute difference " << absDiff << "." << endl;
      fail++;
    }
  CMatrix TRSM16;
  TRSM16.readMatlabFile("testTrsm.mat", "TRSM16");
  F.deepCopy(B);
  F.trsm(U2, alpha, "R", "U", "T", "U");
  absDiff=F.maxAbsDiff(TRSM16);
  if(absDiff<tolInv)
    cout << "trsm rutu matches." << endl;
  else
    {
      cout << "FAILURE: trsm rutu, absolute difference " << absDiff << "." << endl;
      fail++;
    }

  
  return fail;
}


int testAxpy()
{
  int fail = 0;
  CMatrix iMat;
  iMat.readMatlabFile("testAxpy.mat", "i");
  int i = (int)iMat.getVal(0) - 1;
  CMatrix jMat;
  jMat.readMatlabFile("testAxpy.mat", "j");
  int j = (int)jMat.getVal(0) - 1;
  CMatrix kMat;
  kMat.readMatlabFile("testAxpy.mat", "k");
  int k = (int)kMat.getVal(0) - 1;
  CMatrix alphaMat;
  alphaMat.readMatlabFile("testAxpy.mat", "alpha");
  double alpha = alphaMat.getVal(0);
  CMatrix A;
  A.readMatlabFile("testAxpy.mat", "A");
  CMatrix B;
  B.readMatlabFile("testAxpy.mat", "B");
  CMatrix C;
  C.readMatlabFile("testAxpy.mat", "C");
  CMatrix D;
  D.readMatlabFile("testAxpy.mat", "D");
  CMatrix AXPY1;
  AXPY1.readMatlabFile("testAxpy.mat", "AXPY1");
  CMatrix F;
  F.deepCopy(B);
  F.axpyRowRow(i, A, k, alpha);
  if(F.equals(AXPY1))
    cout << "axpyRowRow matches." << endl;
  else
    {
      cout << "FAILURE: axpyRowRow." << endl;
      fail++;
    }
  CMatrix AXPY2;
  AXPY2.readMatlabFile("testAxpy.mat", "AXPY2");
  F.deepCopy(B);
  F.axpyRowCol(i, C, j, alpha);
  if(F.equals(AXPY2))
    cout << "axpyRowCol matches." << endl;
  else
    {
      cout << "FAILURE: axpyRowCol." << endl;
      fail++;
    }
  CMatrix AXPY3;
  AXPY3.readMatlabFile("testAxpy.mat", "AXPY3");
  F.deepCopy(B);
  F.axpyColCol(j, A, k, alpha);
  if(F.equals(AXPY3))
    cout << "axpyColCol matches." << endl;
  else
    {
      cout << "FAILURE: axpyColCol." << endl;
      fail++;
    }
  CMatrix AXPY4;
  AXPY4.readMatlabFile("testAxpy.mat", "AXPY4");
  F.deepCopy(B);
  F.axpyColRow(j, C, i, alpha);
  if(F.equals(AXPY4))
    cout << "axpyColRow matches." << endl;
  else
    {
      cout << "FAILURE: axpyColRow." << endl;
      fail++;
    }
  
  CMatrix AXPY5;
  AXPY5.readMatlabFile("testAxpy.mat", "AXPY5");
  F.deepCopy(D);
  F.axpyDiagRow(C, i, alpha);
  if(F.equals(AXPY5))
    cout << "axpyDiagRow matches." << endl;
  else
    {
      cout << "FAILURE: axpyDiagRow." << endl;
      fail++;
    }
  CMatrix AXPY6;
  AXPY6.readMatlabFile("testAxpy.mat", "AXPY6");
  F.deepCopy(D);
  F.axpyDiagCol(B, j, alpha);
  if(F.equals(AXPY6))
    cout << "axpyDiagCol matches." << endl;
  else
    {
      cout << "FAILURE: axpyDiagCol." << endl;
      fail++;
    }
  return fail;
}

int testGemv()
{
  int fail = 0;
  CMatrix iMat;
  iMat.readMatlabFile("testGemv.mat", "i");
  int i = (int)iMat.getVal(0) - 1;
  CMatrix jMat;
  jMat.readMatlabFile("testGemv.mat", "j");
  int j = (int)jMat.getVal(0) - 1;
  CMatrix kMat;
  kMat.readMatlabFile("testGemv.mat", "k");
  int k = (int)kMat.getVal(0) - 1;
  CMatrix alphaMat;
  alphaMat.readMatlabFile("testGemv.mat", "alpha");
  double alpha = alphaMat.getVal(0);
  CMatrix betaMat;
  betaMat.readMatlabFile("testGemv.mat", "beta");
  double beta = betaMat.getVal(0);
  CMatrix A;
  A.readMatlabFile("testGemv.mat", "A");
  CMatrix B;
  B.readMatlabFile("testGemv.mat", "B");
  CMatrix C;
  C.readMatlabFile("testGemv.mat", "C");
  CMatrix D;
  D.readMatlabFile("testGemv.mat", "D");
  CMatrix E;
  E.readMatlabFile("testGemv.mat", "E");
  CMatrix GEMV1;
  GEMV1.readMatlabFile("testGemv.mat", "GEMV1");
  CMatrix GEMV2;
  GEMV2.readMatlabFile("testGemv.mat", "GEMV2");
  CMatrix GEMV3;
  GEMV3.readMatlabFile("testGemv.mat", "GEMV3");
  CMatrix GEMV4;
  GEMV4.readMatlabFile("testGemv.mat", "GEMV4");
  CMatrix GEMV5;
  GEMV5.readMatlabFile("testGemv.mat", "GEMV5");
  CMatrix GEMV6;
  GEMV6.readMatlabFile("testGemv.mat", "GEMV6");
  CMatrix GEMV7;
  GEMV7.readMatlabFile("testGemv.mat", "GEMV7");
  CMatrix GEMV8;
  GEMV8.readMatlabFile("testGemv.mat", "GEMV8");
  
  CMatrix F(B);
  F.gemvRowRow(i, C, D, k, alpha, beta, "n");
  if(F.equals(GEMV1))
    cout << "gemvRowRow n matches." << endl;
  else
  {
      cout << "FAILURE: gemvRowRow n." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.gemvRowRow(i, E, D, k, alpha, beta, "t");
  if(F.equals(GEMV2))
    cout << "gemvRowRow t matches." << endl;
  else
  {
      cout << "FAILURE: gemvRowRow t." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.gemvRowCol(i, C, A, k, alpha, beta, "n");
  if(F.equals(GEMV3))
    cout << "gemvRowCol n matches." << endl;
  else
  {
      cout << "FAILURE: gemvRowCol n." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.gemvRowCol(i, E, A, k, alpha, beta, "t");
  if(F.equals(GEMV4))
    cout << "gemvRowCol t matches." << endl;
  else
  {
      cout << "FAILURE: gemvRowCol t." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.gemvColCol(j, C, A, k, alpha, beta, "n");
  if(F.equals(GEMV5))
    cout << "gemvColCol n matches." << endl;
  else
  {
      cout << "FAILURE: gemvColCol n." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.gemvColCol(j, E, A, k, alpha, beta, "t");
  if(F.equals(GEMV6))
    cout << "gemvColCol t matches." << endl;
  else
  {
      cout << "FAILURE: gemvColCol t." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.gemvColRow(j, C, D, k, alpha, beta, "n");
  if(F.equals(GEMV7))
    cout << "gemvColRow n matches." << endl;
  else
  {
      cout << "FAILURE: gemvColRow n." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.gemvColRow(j, E, D, k, alpha, beta, "t");
  if(F.equals(GEMV8))
    cout << "gemvColRow t matches." << endl;
  else
  {
      cout << "FAILURE: gemvColRow t." << endl;
      fail++;
  }
  return fail;
}
int testGer()
{
  int fail = 0;
  CMatrix iMat;
  iMat.readMatlabFile("testGer.mat", "i");
  int i = (int)iMat.getVal(0) - 1;
  CMatrix jMat;
  jMat.readMatlabFile("testGer.mat", "j");
  int j = (int)jMat.getVal(0) - 1;
  CMatrix kMat;
  kMat.readMatlabFile("testGer.mat", "k");
  int k = (int)kMat.getVal(0) - 1;
  CMatrix alphaMat;
  alphaMat.readMatlabFile("testGer.mat", "alpha");
  double alpha = alphaMat.getVal(0);
  CMatrix x;
  x.readMatlabFile("testGer.mat", "x");
  CMatrix y;
  y.readMatlabFile("testGer.mat", "y");
  CMatrix A;
  A.readMatlabFile("testGer.mat", "A");
  CMatrix B;
  B.readMatlabFile("testGer.mat", "B");
  CMatrix C;
  C.readMatlabFile("testGer.mat", "C");
  CMatrix D;
  D.readMatlabFile("testGer.mat", "D");
  CMatrix GER1;
  GER1.readMatlabFile("testGer.mat", "GER1");
  CMatrix GER2;
  GER2.readMatlabFile("testGer.mat", "GER2");
  CMatrix GER3;
  GER3.readMatlabFile("testGer.mat", "GER3");
  CMatrix GER4;
  GER4.readMatlabFile("testGer.mat", "GER4");
  CMatrix GER5;
  GER5.readMatlabFile("testGer.mat", "GER5");
  
  CMatrix F;
  F.deepCopy(A);
  F.ger(x, y, alpha);
  if(F.equals(GER1))
    cout << "ger matches." << endl;
  else
  {
      cout << "FAILURE: ger." << endl;
      fail++;
  }
  F.deepCopy(A);
  F.gerRowRow(C, i, B, k, alpha);
  if(F.equals(GER2))
    cout << "gerRowRow matches." << endl;
  else
  {
      cout << "FAILURE: gerRowRow." << endl;
      fail++;
  }
  F.deepCopy(A);
  F.gerRowCol(C, i, D, j, alpha);
  if(F.equals(GER3))
    cout << "gerRowCol matches." << endl;
  else
    {
      cout << "FAILURE: gerRowCol." << endl;
      fail++;
    }
  F.deepCopy(A);
  F.gerColCol(B, j, D, k, alpha);
  if(F.equals(GER4))
    cout << "gerColCol matches." << endl;
  else
    {
      cout << "FAILURE: gerColCol." << endl;
      fail++;
    }
  F.deepCopy(A);
  F.gerColRow(B, j, B, i, alpha);
  if(F.equals(GER5))
    cout << "gerColRow matches." << endl;
  else
    {
      cout << "FAILURE: gerColRow." << endl;
      fail++;
    }
  return fail;
  
}

int testSyr()
{
  int fail = 0;
  CMatrix iMat;
  iMat.readMatlabFile("testSyr.mat", "i");
  int i = (int)iMat.getVal(0) - 1;
  CMatrix jMat;
  jMat.readMatlabFile("testSyr.mat", "j");
  int j = (int)jMat.getVal(0) - 1;
  CMatrix alphaMat;
  alphaMat.readMatlabFile("testSyr.mat", "alpha");
  double alpha = alphaMat.getVal(0);
  CMatrix x;
  x.readMatlabFile("testSyr.mat", "x");
  CMatrix A;
  A.readMatlabFile("testSyr.mat", "A");
  A.setSymmetric(true);
  CMatrix B;
  B.readMatlabFile("testSyr.mat", "B");
  CMatrix SYR1;
  SYR1.readMatlabFile("testSyr.mat", "SYR1");
  CMatrix SYR2;
  SYR2.readMatlabFile("testSyr.mat", "SYR2");
  CMatrix SYR3;
  SYR3.readMatlabFile("testSyr.mat", "SYR3");
  
  CMatrix F;
  F.deepCopy(A);
  F.syr(x, alpha, "U");
  if(F.equals(SYR1))
    cout << "syr u matches." << endl;
  else
  {
      cout << "FAILURE: syr u." << endl;
      fail++;
  }
  F.deepCopy(A);
  F.syrRow(B, i, alpha, "U");
  if(F.equals(SYR2))
    cout << "syrRow u matches." << endl;
  else
  {
      cout << "FAILURE: syrRow u." << endl;
      fail++;
  }
  F.deepCopy(A);
  F.syrRow(B, i, alpha, "L");
  if(F.equals(SYR2))
    cout << "syrRow l matches." << endl;
  else
  {
      cout << "FAILURE: syrRow l." << endl;
      fail++;
  }
  F.deepCopy(A);
  F.syrCol(B, j, alpha, "u");
  if(F.equals(SYR3))
    cout << "syrCol u matches." << endl;
  else
    {
      cout << "FAILURE: syrCol u." << endl;
      fail++;
    }
  F.deepCopy(A);
  F.syrCol(B, j, alpha, "l");
  if(F.equals(SYR3))
    cout << "syrCol l matches." << endl;
  else
    {
      cout << "FAILURE: syrCol l." << endl;
      fail++;
    }
  return fail;
  
}

int testSymv()
{
  int fail = 0;
  CMatrix iMat;
  iMat.readMatlabFile("testSymv.mat", "i");
  int i = (int)iMat.getVal(0) - 1;
  CMatrix jMat;
  jMat.readMatlabFile("testSymv.mat", "j");
  int j = (int)jMat.getVal(0) - 1;
  CMatrix kMat;
  kMat.readMatlabFile("testSymv.mat", "k");
  int k = (int)kMat.getVal(0) - 1;
  CMatrix alphaMat;
  alphaMat.readMatlabFile("testSymv.mat", "alpha");
  double alpha = alphaMat.getVal(0);
  CMatrix betaMat;
  betaMat.readMatlabFile("testSymv.mat", "beta");
  double beta = betaMat.getVal(0);
  CMatrix A;
  A.readMatlabFile("testSymv.mat", "A");
  CMatrix B;
  B.readMatlabFile("testSymv.mat", "B");
  CMatrix C;
  C.readMatlabFile("testSymv.mat", "C");
  C.setSymmetric(true);
  CMatrix D;
  D.readMatlabFile("testSymv.mat", "D");
  CMatrix SYMV1;
  SYMV1.readMatlabFile("testSymv.mat", "SYMV1");
  CMatrix SYMV2;
  SYMV2.readMatlabFile("testSymv.mat", "SYMV2");
  CMatrix SYMV3;
  SYMV3.readMatlabFile("testSymv.mat", "SYMV3");
  CMatrix SYMV4;
  SYMV4.readMatlabFile("testSymv.mat", "SYMV4");
  CMatrix SYMV5;
  SYMV5.readMatlabFile("testSymv.mat", "SYMV5");
  CMatrix SYMV6;
  SYMV6.readMatlabFile("testSymv.mat", "SYMV6");
  CMatrix SYMV7;
  SYMV7.readMatlabFile("testSymv.mat", "SYMV7");
  CMatrix SYMV8;
  SYMV8.readMatlabFile("testSymv.mat", "SYMV8");
  
  CMatrix F(B);
  F.symvRowRow(i, C, D, k, alpha, beta, "u");
  if(F.equals(SYMV1))
    cout << "symvRowRow u matches." << endl;
  else
  {
      cout << "FAILURE: symvRowRow u." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.symvRowRow(i, C, D, k, alpha, beta, "l");
  if(F.equals(SYMV2))
    cout << "symvRowRow l matches." << endl;
  else
  {
      cout << "FAILURE: symvRowRow l." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.symvRowCol(i, C, A, k, alpha, beta, "u");
  if(F.equals(SYMV3))
    cout << "symvRowCol u matches." << endl;
  else
  {
      cout << "FAILURE: symvRowCol u." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.symvRowCol(i, C, A, k, alpha, beta, "l");
  if(F.equals(SYMV4))
    cout << "symvRowCol l matches." << endl;
  else
  {
      cout << "FAILURE: symvRowCol l." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.symvColCol(j, C, A, k, alpha, beta, "u");
  if(F.equals(SYMV5))
    cout << "symvColCol u matches." << endl;
  else
  {
      cout << "FAILURE: symvColCol u." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.symvColCol(j, C, A, k, alpha, beta, "l");
  if(F.equals(SYMV6))
    cout << "symvColCol l matches." << endl;
  else
  {
      cout << "FAILURE: symvColCol l." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.symvColRow(j, C, D, k, alpha, beta, "u");
  if(F.equals(SYMV7))
    cout << "symvColRow u matches." << endl;
  else
  {
      cout << "FAILURE: symvColRow u." << endl;
      fail++;
  }
  F.deepCopy(B);
  F.symvColRow(j, C, D, k, alpha, beta, "l");
  if(F.equals(SYMV8))
    cout << "symvColRow l matches." << endl;
  else
  {
      cout << "FAILURE: symvColRow l." << endl;
      fail++;
  }
  return fail;
}

int testScale()
{
  int fail = 0;
  CMatrix iMat;
  iMat.readMatlabFile("testScale.mat", "i");
  int i = (int)iMat.getVal(0) - 1;
  CMatrix jMat;
  jMat.readMatlabFile("testScale.mat", "j");
  int j = (int)jMat.getVal(0) - 1;
  CMatrix alphaMat;
  alphaMat.readMatlabFile("testScale.mat", "alpha");
  double alpha = alphaMat.getVal(0);
  CMatrix A;
  A.readMatlabFile("testScale.mat", "A");
  CMatrix B;
  B.readMatlabFile("testScale.mat", "B");
  CMatrix C;
  C.readMatlabFile("testScale.mat", "C");
  CMatrix SCALE1;
  SCALE1.readMatlabFile("testScale.mat", "SCALE1");
  CMatrix SCALE2;
  SCALE2.readMatlabFile("testScale.mat", "SCALE2");
  CMatrix SCALE3;
  SCALE3.readMatlabFile("testScale.mat", "SCALE3");
  
  A.scale(alpha);
  if(A.equals(SCALE1))
    cout << "scale matches." << endl;
  else
  {
      cout << "FAILURE: scale." << endl;
      fail++;
  }
  B.scaleRow(i, alpha);
  if(B.equals(SCALE2))
    cout << "scaleRow matches." << endl;
  else
  {
      cout << "FAILURE: scaleRow." << endl;
      fail++;
  }
  C.scaleCol(j, alpha);
  if(C.equals(SCALE3))
    cout << "scaleCol matches." << endl;
  else
  {
      cout << "FAILURE: scaleCol." << endl;
      fail++;
  }
  return fail;
  

}
/*
  assert(max(A-Ainv)
  cout << "A = " << endl << A;

  // create a positive definite matrix
  double b[9] = {18.0455,    4.9297,    3.6417,    
		 4.9297,    5.9624,   -4.7233,    
		 3.6417,   -4.7233,    7.8203};
  CMatrix B(nrows, ncols, b);
  cout << "Matrix B = " << endl << B;
  CMatrix C(B);
//CMatrix C(nrows, ncols);
  //  C.deepCopy(B);
  //C.chol();
  cout << "Cholesky decomposition of B = " << endl << C;
  CMatrix D(nrows, ncols);
  D.deepCopy(B);
  D.pdinv();
  cout << "Inverse of B using pdinv = " << endl << D;
  CMatrix E(nrows, ncols); 
  E.deepCopy(B);
  E.inv();
  cout << "Inverse of B using inv = " << endl << E;
  CMatrix F(nrows, ncols);
  F.gemm(B, E, 1.0, 0.0, "n", "n");
  cout << "B multiplied by its inverse = " << endl << F;
  CMatrix G(nrows, ncols);
  G.gemm(B, D, 1.0, 0.0, "n", "n");
  cout << "B multiplied by its pdinverse = " << endl << G;
  double trA = trace(A);
  cout << "Trace of A is " << trA << endl;
  cout << "B is " << B << endl;

  double x[60] = {-0.4326,    0.2944,   -1.6041,
		  -1.6656,   -1.3362,   0.2573,
		  0.1253,    0.7143,   -1.0565,
		  0.2877,    1.6236,    1.4151,
		  -1.1465,   -0.6918,   -0.8051,
		  1.1909,    0.8580,    0.5287,
		  1.1892,    1.2540,    0.2193,
		  -0.0376,   -1.5937,   -0.9219,
		  0.3273,   -1.4410,   -2.1707,
		  0.1746,    0.5711,   -0.0592,
		  -0.1867,   -0.3999,   -1.0106,
		  0.7258,    0.6900,    0.6145,
		  -0.5883,    0.8156,    0.5077,
		  2.1832,    0.7119,    1.6924,
		  -0.1364,    1.2902,    0.5913,
		  0.1139,    0.6686,   -0.6436,
		  1.0668,    1.1908,    0.3803,
		  0.0593,   -1.2025,   -1.0091,
		  -0.0956,   -0.0198,   -0.0195,
		  -0.8323,   -0.1567,   -0.0482};
  nrows=3;
  ncols=20;
  CMatrix X(nrows, ncols, x);
  X.trans();
  ncols = 3;
  nrows = 20;
  cout << "X: " << endl << X << endl;
  vector<int> rows;
  rows.push_back(0);
  rows.push_back(3);
  rows.push_back(19);
  cout << "Extracting  rows." << endl;
  CMatrix X2(rows.size(), ncols);
  X.getMatrix(X2, rows, 0, ncols-1);
  cout << "Extracted matrix: " << endl << X2;
  CMatrix X3(X.getRows(), X.getCols());
  X3.deepCopy(X);
  X3.appendCols(X);
  X3.appendCols(X);
  cout << X3.getCols() << endl;
  cout << X3.getRows() << endl;
  //  X3.appendRows(X);
  cout << X3;
  vector<int> cols;
  cols.push_back(0);
  cols.push_back(4);
  cols.push_back(8);
  CMatrix X4(nrows, cols.size());
  X3.getMatrix(X4, 0, nrows-1, cols);
  X4-=X;
  cout<<X4;
  CMatrix cov(ncols, ncols);
  cov.gemm(X, X, 1.0, 0.0, "t", "n");
  cout<< "Covariance " << endl << cov;
  
  // now transpose x and use gemm with the same operation.
  X.trans();
  CMatrix cov2(ncols, ncols);
  cov2.gemm(X, X, 1.0, 0.0, "n", "t");
  cout<< "Covariance 2" << endl << cov2;
  
  // this should be the same as the other 2 but it isn't"!
  cov.dsyrk(X, "U", 1.0, 0.0, "n");
  cout << "Covariance 3" << endl << cov;
  cout << "Covariance 4" << endl << multiply(X, "n", X, "t") << endl;
  CMatrix X5(X.getRows(), X.getCols());
  X5.deepCopy(X);
  X5.trans();
  cout << "Covariance 5" << endl << multiply(X, X5) << endl;
  cout << "Sum over all X: " << sum(X) << endl;
  
  cout << "Kernel 1" << endl << multiply(X, "t", X, "n") << endl;
 
  //  X.trans();
  cout << "Row sum of X: " << endl <<sumRow(X) << endl;
  cout << "Row mean of X: " << endl <<meanRow(X) << endl;
  cout << "Col sum of X: " << endl << sumCol(X) << endl;
  cout << "Col mean of X: " << endl << meanCol(X) << endl;

  CMatrix XBar(X.getRows(), X.getCols());
  XBar.deepCopy(X);
  XBar-=meanRow(X);
  cout << "Centred X: " << endl << XBar << endl;
  cout << "Mean of XBar: " << endl << meanRow(XBar) << endl;

  CMatrix randMat(100000, 10);
  randMat.randn();
  //  cout << "random matrix: " << endl << randMat << endl;
  cout << "mean of matrix: " << endl << meanRow(randMat) << endl;
  randMat *= randMat;
  cout << "var of matrix: " << endl << meanRow(randMat) << endl;

  // test transpose and swap column and row.
  double x2[21] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,16, 17, 18, 19, 20};
  nrows=3;
  ncols=7;
  CMatrix X6(nrows, ncols, x2);
  CMatrix X7(nrows, ncols, x2);
  X5.trans();
  cout << X6;
  cout << endl << " swap rows 1 and 3 " << endl;
  X6.swapRows(0, 2);
  cout << X6;
  X6.trans();
  cout <<endl << " Transpose of matrix " << endl;
  cout << X6;
  cout << endl << " swap columns 1 and 3 " << endl;
  X6.swapCols(0, 2);
  cout << X6;
  X6.trans();
  cout << endl << " Transpose of matrix " << endl;
  cout << X6;
  X6.subtract(X7);
  assert(abs(max(X6))<EPS);

  // Take x7, scale 3rd column by 2 transpose scale 3rd row by 0.5.
  X7.trans();
  CMatrix X8(X7);
  X8.scale(2.0);
  cout <<  X7.dotCol(2, X8, 0);
*/
