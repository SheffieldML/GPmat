/* Copyright (c) 2010 Antti Honkela <antti.honkela@tkk.fi>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. */
#include <math.h>
#include "mex.h"
#include "matrix.h"

mxArray *one;

double xrandn(void)
{
  mxArray *tmp;
  double val;

  mexCallMATLAB(1, &tmp, 1, &one, "randn");
  val = mxGetScalar(tmp);
  mxDestroyArray(tmp);
  
  return val;
}


double xrand(void)
{
  mxArray *tmp;
  double val;

  mexCallMATLAB(1, &tmp, 1, &one, "rand");
  val = mxGetScalar(tmp);
  mxDestroyArray(tmp);

  return val;
}


/* Based on Marsaglia and Tsang, "A Simple Method for generating gamma
 * variables", ACM Transactions on Mathematical Software 26(3):363-372
 * (2000).
 */
double xgamrnd(const double a, const double b)
{
  double d,c,x,v,u,multiplier=1.0;
  d = (a<1.0 ? a+1.0 : a) - 1.0/3.0;
  c = 1.0/sqrt(9.0*d);
  for(;;) {
    do {
      x=xrandn();
      v=1.+c*x;
    } while(v <= 0.);
    v = v*v*v;
    u = xrand();
    if (u < 1.0-0.0331*(x*x)*(x*x))
      break;
    if (log(u) <  0.5*x*x + d*(1.0-v+log(v)))
      break;
  }
  if (a < 1.0) {
    u = xrand();
    multiplier = pow(u, 1.0 / a);
  }

  return multiplier*b*d*v;
}



void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int dim1, dim2, dim1_1, dim2_1, dim1_2, dim2_2, i;
  int inc1, inc2;
  double *in1, *in2, *out, *signs;
  double v1, v2, temp;

  /******************** input variables ********************/

  in1     = mxGetPr(prhs[0]);
  in2     = mxGetPr(prhs[1]);
  dim1_1  = (int)mxGetM(prhs[0]);
  dim2_1  = (int)mxGetN(prhs[0]);
  dim1_2  = (int)mxGetM(prhs[1]);
  dim2_2  = (int)mxGetN(prhs[1]);

  dim1 = fmax(dim1_1, dim1_2);
  dim2 = fmax(dim2_1, dim2_2);

  /* Check that inputs are either a scalar and a matrix or matrices of
     the same size */
  if (!(((dim1_1 == 1 && dim2_1 == 1) || ((dim1_1 == dim1 && dim2_1 == dim2))) &&
	((dim1_2 == 1 && dim2_2 == 1) || ((dim1_2 == dim1 && dim2_2 == dim2)))))
    mexErrMsgTxt("xgamrnd: invalid input sizes");

  /* Increments to walk through the inputs */
  if (dim1_1 == 1 && dim2_1 == 1)
    inc1 = 0;
  else
    inc1 = 1;

  if (dim1_2 == 1 && dim2_2 == 1)
    inc2 = 0;
  else
    inc2 = 1;

  /* Create outputs, complex if single output, real otherwise */
  plhs[0] = mxCreateDoubleMatrix(dim1, dim2, mxREAL);
  out     = mxGetPr(plhs[0]);

  one = mxCreateDoubleScalar(1);

  /* Do the hard work */
  for (i=0; i<dim1*dim2; i++) {
    v1 = *in1;
    in1 += inc1;
    v2 = *in2;
    in2 += inc2;

    *out++ = xgamrnd(v1, v2);
  }

  mxDestroyArray(one);

  return;
}
