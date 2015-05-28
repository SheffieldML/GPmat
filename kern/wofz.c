/* wofz.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/
#include "mex.h"

#include "f2c.h"

int wofz(doublereal *xi, 
	  doublereal *yi, 
	  doublereal *u, 
	  doublereal *v, 
	  logical *flag__)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer i_dnnt(doublereal *);
    double exp(doublereal), cos(doublereal), sin(doublereal), pow_di(
	    doublereal *, integer *);

    /* Local variables */
    static logical a, b;
    static doublereal c__, h__;
    static integer i__, j, n;
    static doublereal x, y, h2, u1, v1, u2, v2, w1;
    static integer nu;
    static doublereal rx, ry, sx, sy, tx, ty;
    static integer np1, kapn;
    static doublereal xabs, yabs, daux, qrho, xaux, xsum, ysum, xabsq, xquad, 
	    yquad, qlambda;


/*  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES */
/*  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z), */
/*  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I */
/*  MEANS SQRT(-1). */
/*  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT */
/*  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT */
/*  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO */
/*  OF THE FUNCTION. */
/*  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION. */


/*  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS : */
/*     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF */
/*                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE */
/*                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION */
/*                FLOATING-POINT ARITHMETIC */
/*     RMAXEXP  = LN(RMAX) - LN(2) */
/*     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION */
/*                GONIOMETRIC FUNCTION (DCOS, DSIN, ...) */
/*  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL */
/*  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS */


/*  PARAMETER LIST */
/*     XI     = REAL      PART OF Z */
/*     YI     = IMAGINARY PART OF Z */
/*     U      = REAL      PART OF W(Z) */
/*     V      = IMAGINARY PART OF W(Z) */
/*     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL */
/*              OCCUR OR NOT; TYPE LOGICAL; */
/*              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING */
/*              MEANING : */
/*              FLAG=.FALSE. : NO ERROR CONDITION */
/*              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE */
/*                             BECOMES INACTIVE */
/*  XI, YI      ARE THE INPUT-PARAMETERS */
/*  U, V, FLAG  ARE THE OUTPUT-PARAMETERS */

/*  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI) */

/*  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE */
/*  PUT TO 0 UPON UNDERFLOW; */

/*  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF */
/*  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE. */







    *flag__ = FALSE_;

    xabs = abs(*xi);
    yabs = abs(*yi);
    x = xabs / 6.3f;
    y = yabs / 4.4f;


/*     THE FOLLOWING IF-STATEMENT PROTECTS */
/*     QRHO = (X**2 + Y**2) AGAINST OVERFLOW */

    if (xabs > 5e153 || yabs > 5e153) {
	goto L100;
    }

/* Computing 2nd power */
    d__1 = x;
/* Computing 2nd power */
    d__2 = y;
    qrho = d__1 * d__1 + d__2 * d__2;

/* Computing 2nd power */
    d__1 = xabs;
    xabsq = d__1 * d__1;
/* Computing 2nd power */
    d__1 = yabs;
    xquad = xabsq - d__1 * d__1;
    yquad = xabs * 2 * yabs;

    a = qrho < .085264;

    if (a) {

/*  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED */
/*  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297) */
/*  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED */
/*  ACCURACY */

	qrho = (1 - y * .85f) * sqrt(qrho);
	d__1 = qrho * 72 + 6;
	n = i_dnnt(&d__1);
	j = (n << 1) + 1;
	xsum = 1.f / j;
	ysum = 0.;
	for (i__ = n; i__ >= 1; --i__) {
	    j += -2;
	    xaux = (xsum * xquad - ysum * yquad) / i__;
	    ysum = (xsum * yquad + ysum * xquad) / i__;
	    xsum = xaux + 1.f / j;
/* L10: */
	}
	u1 = (xsum * yabs + ysum * xabs) * -1.12837916709551257388 + 1.f;
	v1 = (xsum * xabs - ysum * yabs) * 1.12837916709551257388;
	daux = exp(-xquad);
	u2 = daux * cos(yquad);
	v2 = -daux * sin(yquad);

	*u = u1 * u2 - v1 * v2;
	*v = u1 * v2 + v1 * u2;

    } else {

/*  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE */
/*  CONTINUED FRACTION */
/*  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED */
/*  ACCURACY */

/*  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED */
/*  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION */
/*  IS USED TO CALCULATE THE DERIVATIVES OF W(Z) */
/*  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED */
/*  TO OBTAIN THE REQUIRED ACCURACY */
/*  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED */
/*  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY */


	if (qrho > 1.f) {
	    h__ = 0.;
	    kapn = 0;
	    qrho = sqrt(qrho);
	    nu = (integer) (1442 / (qrho * 26 + 77) + 3);
	} else {
	    qrho = (1 - y) * sqrt(1 - qrho);
	    h__ = qrho * 1.88f;
	    h2 = h__ * 2;
	    d__1 = qrho * 34 + 7;
	    kapn = i_dnnt(&d__1);
	    d__1 = qrho * 26 + 16;
	    nu = i_dnnt(&d__1);
	}

	b = h__ > 0.f;

	if (b) {
	    qlambda = pow_di(&h2, &kapn);
	}

	rx = 0.f;
	ry = 0.f;
	sx = 0.f;
	sy = 0.f;

	for (n = nu; n >= 0; --n) {
	    np1 = n + 1;
	    tx = yabs + h__ + np1 * rx;
	    ty = xabs - np1 * ry;
/* Computing 2nd power */
	    d__1 = tx;
/* Computing 2nd power */
	    d__2 = ty;
	    c__ = .5f / (d__1 * d__1 + d__2 * d__2);
	    rx = c__ * tx;
	    ry = c__ * ty;
	    if (b && n <= kapn) {
		tx = qlambda + sx;
		sx = rx * tx - ry * sy;
		sy = ry * tx + rx * sy;
		qlambda /= h2;
	    }
/* L11: */
	}

	if (h__ == 0.f) {
	    *u = rx * 1.12837916709551257388;
	    *v = ry * 1.12837916709551257388;
	} else {
	    *u = sx * 1.12837916709551257388;
	    *v = sy * 1.12837916709551257388;
	}

	if (yabs == 0.f) {
/* Computing 2nd power */
	    d__1 = xabs;
	    *u = exp(-(d__1 * d__1));
	}

    }



/*  EVALUATION OF W(Z) IN THE OTHER QUADRANTS */


    if (*yi < 0.f) {

	if (a) {
	    u2 *= 2;
	    v2 *= 2;
	} else {
	    xquad = -xquad;


/*         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2) */
/*         AGAINST OVERFLOW */

	    if (yquad > 3537118876014220. || xquad > 708.503061461606) {
		goto L100;
	    }

	    w1 = exp(xquad) * 2;
	    u2 = w1 * cos(yquad);
	    v2 = -w1 * sin(yquad);
	}

	*u = u2 - *u;
	*v = v2 - *v;
	if (*xi > 0.f) {
	    *v = -(*v);
	}
    } else {
	if (*xi < 0.f) {
	    *v = -(*v);
	}
    }

    return 0;

L100:
    *flag__ = TRUE_;
    return 0;

} /* wofz_ */

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{  
  if(mxGetClassID(prhs[0]) != mxDOUBLE_CLASS)
    mexErrMsgTxt("Error input should be DOUBLE");  
  const int numDims = mxGetNumberOfDimensions(prhs[0]);
  const int numElements = mxGetNumberOfElements(prhs[0]);
  const int* dims = mxGetDimensions(prhs[0]);

  doublereal* xi = mxGetPr(prhs[0]);
  doublereal* yi;
  int i = 0;
  if(!mxIsComplex(prhs[0]))
  {
    yi = mxCalloc(numElements, sizeof(doublereal));
    for(i = 0; i<numElements; i++)
      yi[i] = 0;
  }
  else
    yi = mxGetPi(prhs[0]);
  logical* flag = mxCalloc(numElements, sizeof(logical));
  plhs[0] = mxCreateNumericArray(numDims, dims, mxDOUBLE_CLASS, mxCOMPLEX);
  double* u = mxGetPr(plhs[0]);
  double* v = mxGetPi(plhs[0]);
  for(i = 0; i<numElements; i++)
  {
    wofz(&xi[i], &yi[i], &u[i], &v[i], &flag[i]);
  }
  mxFree(flag);
}
