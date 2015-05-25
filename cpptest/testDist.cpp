#include "CDist.h"

int main()
{
  int nrows = 1;
  int ncols = 9;
  CGaussianDist prior;
  double x[9] = {18.0455,    4.9297,    3.6417,    
		 4.9297,    5.9624,   -4.7233,    
		 3.6417,   -4.7233,    7.8203};
  CMatrix X(nrows, ncols, x);
  cout << "Log probability: ";
  cout << prior.logProb(x[0]) << endl;

}

