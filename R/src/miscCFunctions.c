#include <math.h>


extern void
_ClnDiffErfs(double *in1, double *in2, int *dim1, int *dim2, double *out, int *signs);

/*
 * This function is based on cleaning up f2c translation of Fortran
 * code found at Netlib/Specfun.  The original Fortran code was
 * written by W. J. Cody, 1990.  The translation was performed by
 * Antti Honkela, 2008.
 */

double calerf(double x, int jint)
{
  /* constants */
  double sqrpi = 5.6418958354775628695E-1; /* 1 / sqrt(M_PI); */
  double xinf = 1.79e308;
  double xneg = -26.628;
  double xsmall = 1.11e-16;
  double xbig = 26.543;
  double xhuge = 6.71e7;
  double xmax = 2.53e307;

  static double a[] = {3.16112374387056560e00,1.13864154151050156e02, 
		       3.77485237685302021e02,3.20937758913846947e03,
		       1.85777706184603153e-1};
  static double b[] = {2.36012909523441209e01,2.44024637934444173e02,
		       1.28261652607737228e03,2.84423683343917062e03};

  static double c[] = {5.64188496988670089E-1,8.88314979438837594E0,
		       6.61191906371416295E01,2.98635138197400131E02,
		       8.81952221241769090E02,1.71204761263407058E03,
		       2.05107837782607147E03,1.23033935479799725E03,
		       2.15311535474403846E-8};

  static double d[] = {1.57449261107098347E01,1.17693950891312499E02,
		       5.37181101862009858E02,1.62138957456669019E03,
		       3.29079923573345963E03,4.36261909014324716E03,
		       3.43936767414372164E03,1.23033935480374942E03};

  static double p[] = {3.05326634961232344E-1,3.60344899949804439E-1,
		       1.25781726111229246E-1,1.60837851487422766E-2,
		       6.58749161529837803E-4,1.63153871373020978E-2};
  
  static double q[] = {2.56852019228982242E00,1.87295284992346047E00,
		       5.27905102951428412E-1,6.05183413124413191E-2,
		       2.33520497626869185E-3};

  /* Local variables */
  double xden;
  double xnum;
  double y, del, ysq;
  double result;

  int i;

  /* ------------------------------------------------------------------ */

  /* This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x) */
  /*   for a real argument  x.  It contains three FUNCTION type */
  /*   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX), */
  /*   and one SUBROUTINE type subprogram, CALERF.  The calling */
  /*   statements for the primary entries are: */

  /*                   Y=ERF(X)     (or   Y=DERF(X)), */

  /*                   Y=ERFC(X)    (or   Y=DERFC(X)), */
  /*   and */
  /*                   Y=ERFCX(X)   (or   Y=DERFCX(X)). */

  /*   The routine  CALERF  is intended for internal packet use only, */
  /*   all computations within the packet being concentrated in this */
  /*   routine.  The function subprograms invoke  CALERF  with the */
  /*   statement */

  /*          CALL CALERF(ARG,RESULT,JINT) */

  /*   where the parameter usage is as follows */

  /*      Function                     Parameters for CALERF */
  /*       call              ARG                  Result          JINT */

  /*     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0 */
  /*     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1 */
  /*     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2 */

  /*   The main computation evaluates near-minimax approximations */
  /*   from "Rational Chebyshev approximations for the error function" */
  /*   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This */
  /*   transportable program uses rational functions that theoretically */
  /*   approximate  erf(x)  and  erfc(x)  to at least 18 significant */
  /*   decimal digits.  The accuracy achieved depends on the arithmetic */
  /*   system, the compiler, the intrinsic functions, and proper */
  /*   selection of the machine-dependent constants. */

  /* ******************************************************************* */
  /* ******************************************************************* */

  /* Explanation of machine-dependent constants */

  /*   XMIN   = the smallest positive floating-point number. */
  /*   XINF   = the largest positive finite floating-point number. */
  /*   XNEG   = the largest negative argument acceptable to ERFCX; */
  /*            the negative of the solution to the equation */
  /*            2*exp(x*x) = XINF. */
  /*   XSMALL = argument below which erf(x) may be represented by */
  /*            2*x/sqrt(pi)  and above which  x*x  will not underflow. */
  /*            A conservative value is the largest machine number X */
  /*            such that   1.0 + X = 1.0   to machine precision. */
  /*   XBIG   = largest argument acceptable to ERFC;  solution to */
  /*            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where */
  /*            W(x) = exp(-x*x)/[x*sqrt(pi)]. */
  /*   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to */
  /*            machine precision.  A conservative value is */
  /*            1/[2*sqrt(XSMALL)] */
  /*   XMAX   = largest acceptable argument to ERFCX; the minimum */
  /*            of XINF and 1/[sqrt(pi)*XMIN]. */

  /* ******************************************************************* */

  /* Error returns */

  /*  The program returns  ERFC = 0      for  ARG .GE. XBIG; */

  /*                       ERFCX = XINF  for  ARG .LT. XNEG; */
  /*      and */
  /*                       ERFCX = 0     for  ARG .GE. XMAX. */


  /*  Author: W. J. Cody */
  /*          Mathematics and Computer Science Division */
  /*          Argonne National Laboratory */
  /*          Argonne, IL 60439 */

  /*  Latest modification: March 19, 1990 */

  /* ------------------------------------------------------------------ */
  /* S    REAL */
  /* D    DOUBLE PRECISION */
  /* ------------------------------------------------------------------ */
  y = fabs(x);
  if (y <= 0.46875) {
    /* ------------------------------------------------------------------ */
    /*  Evaluate  erf  for  |X| <= 0.46875 */
    /* ------------------------------------------------------------------ */
    ysq = 0.0;
    if (y > xsmall) {
      ysq = y * y;
    }
    xnum = a[4] * ysq;
    xden = ysq;
    for (i = 0; i < 3; ++i) {
      xnum = (xnum + a[i]) * ysq;
      xden = (xden + b[i]) * ysq;
      /* L20: */
    }
    result = x * (xnum + a[3]) / (xden + b[3]);
    if (jint != 0) {
      result = 1.0 - result;
    }
    if (jint == 2) {
      result = exp(ysq) * result;
    }
    goto L800;
    /* ------------------------------------------------------------------ */
    /*  Evaluate  erfc  for 0.46875 <= |X| <= 4.0 */
    /* ------------------------------------------------------------------ */
  } else if (y <= 4.0) {
    xnum = c[8] * y;
    xden = y;
    for (i = 0; i < 7; ++i) {
      xnum = (xnum + c[i]) * y;
      xden = (xden + d[i]) * y;
      /* L120: */
    }
    result = (xnum + c[7]) / (xden + d[7]);
    if (jint != 2) {
      ysq = rint(y * 16.0) / 16.0;
      del = (y - ysq) * (y + ysq);
      result = exp(-ysq * ysq) * exp(-del) * result;
    }
    /* ------------------------------------------------------------------ */
    /*  Evaluate  erfc  for |X| > 4.0 */
    /* ------------------------------------------------------------------ */
  } else {
    result = 0.0;
    if (y >= xbig) {
      if (jint != 2 || y >= xmax) {
	goto L300;
      }
      if (y >= xhuge) {
	result = sqrpi / y;
	goto L300;
      }
    }
    ysq = 1.0 / (y * y);
    xnum = p[5] * ysq;
    xden = ysq;
    for (i = 0; i < 4; ++i) {
      xnum = (xnum + p[i]) * ysq;
      xden = (xden + q[i]) * ysq;
      /* L240: */
    }
    result = ysq * (xnum + p[4]) / (xden + q[4]);
    result = (sqrpi - result) / y;
    if (jint != 2) {
      ysq = rint(y * 16.0) / 16.0;
      del = (y - ysq) * (y + ysq);
      result = exp(-ysq * ysq) * exp(-del) * result;
    }
  }
  /* ------------------------------------------------------------------ */
  /*  Fix up for negative argument, erf, etc. */
  /* ------------------------------------------------------------------ */
 L300:
  if (jint == 0) {
    result = 0.5 - result + 0.5;
    if (x < 0.0) {
      result = -(result);
    }
  } else if (jint == 1) {
    if (x < 0.0) {
      result = 2.0 - result;
    }
  } else {
    if (x < 0.0) {
      if (x < xneg) {
	result = xinf;
      } else {
	ysq = rint(x * 16.0) / 16.0;
	del = (x - ysq) * (x + ysq);
	y = exp(ysq * ysq) * exp(del);
	result = y + y - result;
      }
    }
  }
 L800:
  return result;
  /* ---------- Last card of CALERF ---------- */
} /* calerf_ */


void
_ClnDiffErfs(double *in1, double *in2, int *dim1, int *dim2, double *out, int *signs) {
  const double xsmall = 2.2204e-16;

  int inc1, inc2, i;
  int d1 = *dim1;
  int d2 = *dim2;
  int dmax = (d1 > d2 ? d1 : d2);
  double v1, v2, temp;

  if (*dim1 == 1)
    inc1 = 0;
  else
    inc1 = 1;

  if (*dim2 == 1)
    inc2 = 0;
  else
    inc2 = 1;

  /* Do the hard work */
  for (i=0; i<dmax; i++) {
    v1 = *in1;
    in1 += inc1;
    v2 = *in2;
    in2 += inc2;

    if (v1 < v2) {
      temp = v1;
      v1 = v2;
      v2 = temp;
      *signs++ = -1;
    }
    else
      *signs++ = 1;
    
    if (v1 * v2 < 0)
      *out++ = log(erf(v1) - erf(v2));
    else if (fabs(v1-v2) < 3*xsmall)
      *out++ = -INFINITY;
    else if (v2 > 0)
      *out++ = log(calerf(v2, 2) - calerf(v1, 2) * exp(v2*v2 - v1*v1) + 1e-30) - v2*v2;
    else
      *out++ = log(calerf(-v1, 2) - calerf(-v2, 2) * exp(v1*v1 - v2*v2) + 1e-30) - v1*v1;
  }
}
