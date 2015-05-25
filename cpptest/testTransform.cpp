#include "CTransform.h"

int main()
{
  int ncols = 9;
  CNegLogLogitTransform trans;
  double x[9] = {1e-10, 0.001, 0.02, 0.5, 1, 2, 4, 100, 1e16};
  double a[9];
  double newx[9];
  double g[9];
  for(int i=0; i<9; i++)
    {
      a[i] = trans.xtoa(x[i]);
      newx[i] = trans.atox(a[i]);
      g[i] = trans.gradfact(a[i]);
    }
  cout << "x:    " << "a:   " << "newx:   " << "g:   " << endl;
  for(int i=0; i<9; i++)
    {
      cout << x[i] << "  " << a[i] << "  "<< newx[i] << "  " << g[i] << endl;
    }

}

