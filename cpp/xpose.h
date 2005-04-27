#ifndef XPOSE_H
#define XPOSE_H

// This file contains header information for various fortran functions
// that are used.

/********************** From the book *****************************/
// http://jin.ece.uiuc.edu/routines/routines.html
// 		*   COMPUTATION OF SPECIAL FUNCTIONS   *
// 		*                                      *
// 		*   Shanjie Zhang and Jianming Jin     *
// 		*                                      *
// 		*   Copyright 1996 by John Wiley &     *
// 		*              Sons, Inc.              *
//
//      All the programs and subroutines contained in this diskette
// are copyrighted.  However, we give permission to the reader who
// purchases this book to incorporate any of these programs into his
// or her programs provided that the copyright is acknowledged.

// computes the psi function putting the answer in y.
extern "C" void psi_(const double& x, // input value
		      double& y);      // output
// computes the gamma function.
extern "C" void lgama_(const int& type, // 1 for gamma 0 for log gamma.
		       const double& x, // input value.
		       double& y);      // output.


/********************** William Cody ****************************/
// This is William Cody's scaled complementary erf implementation 
extern "C" double derfcx_(const double& x);

/********************** ACM Algorithms **************************/
// this is algorithm 488 for GRAND, Gaussian sampling.
extern "C" double dgrand_(const int& N); // N is a dummy parameter

// this is algorithm 467 (which I can't get to work).
extern "C" void dxpose_(double* A, 
		       const int& nrows,
		       const int& ncols,
		       const int& numElements,
		       bool* work,
		       const int& lwork);

// this is algorithm 380 (which should be slower).
extern "C" void dtrans_(double* A,
		       const int& ncols,
		       const int& nrows,
		       const int& numElements,
		       int* work,
		       const int& lwork,
		       int& iok);

// this is algorithm 513 (a revised version of 380).
extern "C" void dtransr_(double* A,
		       const int& ncols,
		       const int& nrows,
		       const int& numElements,
		       int* work,
		       const int& lwork,
		       int& iok);


#endif
