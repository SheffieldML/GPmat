// tihis is algorithm 488 for GRAND, Gaussian sampling.
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
