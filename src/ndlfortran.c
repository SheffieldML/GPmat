/* ndlfortran.f -- translated by f2c (version 20041007).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

#ifdef _FORTRAN_MAIN_FIX
    int MAIN__() {return 0;};
#endif

static integer c__0 = 0;

/* http://jin.ece.uiuc.edu/routines/routines.html */
/* 		**************************************** */
/* 		*           DISK TO ACCOMPANY          * */
/* 		*   COMPUTATION OF SPECIAL FUNCTIONS   * */
/* 		*                                      * */
/* 		*   Shanjie Zhang and Jianming Jin     * */
/* 		*                                      * */
/* 		*   Copyright 1996 by John Wiley &     * */
/* 		*              Sons, Inc.              * */
/* 		*                                      * */
/* 		**************************************** */

/* I. INTRODUCTION */

/*      As stated in the preface of our book "Computation of Special */
/* Functions,"  the purpose of this book is to share with the reader */
/* a set of computer programs (130 in total) which we have developed */
/* during the past several years for computing a variety of  special */
/* mathematical functions.  For your convenience, we attach to the */
/* book this diskette that contains all the computer programs */
/* listed or mentioned in the book. */
/*      In this diskette,  we place all the programs under directory */
/* SMF\PROGRAMS. In order to illustrate the use of these programs */
/* and facilitate your testing of the programs, we wrote a short */
/* simple main program for each program so that you can readily test */
/* them. */
/*      All the programs are written in FORTRAN-77 and tested on PCs */
/* and workstations. Therefore, they should run on any computer with */
/* implementation of the FORTRAN-77 standard. */
/*      Although we have made a great effort to test these programs, */
/* we would not be surprised  to find some errors in them.  We would */
/* appreciate it if you can bring to our attention any errors you find. */
/* You can do this by either writing us directly at the location */
/* (e-mail: j-jin1@uiuc.edu) or writing to the publisher, whose address */
/* appears on the back cover of the book.  However, we must note that */
/* all these programs are sold "as is," and we cannot guarantee to */
/* correct the errors reported by readers on any fixed schedule. */
/*      All the programs and subroutines  contained in this diskette */
/* are copyrighted.   However,  we give permission to the reader who */
/* purchases this book to incorporate any of these programs into his */
/* or her programs provided that the copyright is acknowledged. */

/*       ================================================== */
/*       Purpose: This program computes the psi function */
/*                using subroutine PSI */
/*       Input :  x  --- Argument of psi(x) */
/*       Output:  PS --- psi(x) */
/*       Examples: */
/*                   x          Psi(x) */
/*                ------------------------ */
/*                  .25      -4.227453533 */
/*                  .50      -1.963510026 */
/*                  .75      -1.085860880 */
/*                 1.00       -.577215665 */
/*                 1.25       -.227453533 */
/*                 1.50        .036489974 */
/*                 1.75        .247472454 */
/*                 2.00        .422784335 */
/*       ================================================== */

/* Subroutine */ int psi_(doublereal *x, doublereal *ps)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k, n;
    static doublereal s, a1, a2, a3, a4, a5, a6, a7, a8, x2, el, xa, pi;


/*       ====================================== */
/*       Purpose: Compute the psi function */
/*       Input :  x  --- Argument of psi(x) */
/*       Output:  PS --- psi(x) */
/*       ====================================== */

    xa = abs(*x);
    pi = 3.141592653589793;
    el = .5772156649015329;
    s = 0.;
    if (*x == (doublereal) ((integer) (*x)) && *x <= 0.f) {
	*ps = 1e300;
	return 0;
    } else if (xa == (doublereal) ((integer) xa)) {
	n = (integer) xa;
	i__1 = n - 1;
	for (k = 1; k <= i__1; ++k) {
/* L10: */
	    s += 1. / k;
	}
	*ps = -el + s;
    } else if (xa + .5f == (doublereal) ((integer) (xa + .5f))) {
	n = (integer) (xa - .5f);
	i__1 = n;
	for (k = 1; k <= i__1; ++k) {
/* L20: */
	    s += 1.f / (k * 2. - 1.);
	}
	*ps = -el + s * 2. - 1.386294361119891;
    } else {
	if (xa < 10.f) {
	    n = 10 - (integer) xa;
	    i__1 = n - 1;
	    for (k = 0; k <= i__1; ++k) {
/* L30: */
		s += 1. / (xa + k);
	    }
	    xa += n;
	}
	x2 = 1. / (xa * xa);
	a1 = -.08333333333333;
	a2 = .0083333333333333333;
	a3 = -.0039682539682539683;
	a4 = .0041666666666666667;
	a5 = -.0075757575757575758;
	a6 = .021092796092796093;
	a7 = -.083333333333333333;
	a8 = .4432598039215686;
	*ps = log(xa) - .5 / xa + x2 * (((((((a8 * x2 + a7) * x2 + a6) * x2 + 
		a5) * x2 + a4) * x2 + a3) * x2 + a2) * x2 + a1);
	*ps -= s;
    }
    if (*x < 0.f) {
	*ps = *ps - pi * cos(pi * *x) / sin(pi * *x) - 1. / *x;
    }
    return 0;
} /* psi_ */

/* http://jin.ece.uiuc.edu/routines/routines.html */
/* 		**************************************** */
/* 		*           DISK TO ACCOMPANY          * */
/* 		*   COMPUTATION OF SPECIAL FUNCTIONS   * */
/* 		*                                      * */
/* 		*   Shanjie Zhang and Jianming Jin     * */
/* 		*                                      * */
/* 		*   Copyright 1996 by John Wiley &     * */
/* 		*              Sons, Inc.              * */
/* 		*                                      * */
/* 		**************************************** */

/* I. INTRODUCTION */

/*      As stated in the preface of our book "Computation of Special */
/* Functions,"  the purpose of this book is to share with the reader */
/* a set of computer programs (130 in total) which we have developed */
/* during the past several years for computing a variety of  special */
/* mathematical functions.  For your convenience, we attach to the */
/* book this diskette that contains all the computer programs */
/* listed or mentioned in the book. */
/*      In this diskette,  we place all the programs under directory */
/* SMF\PROGRAMS. In order to illustrate the use of these programs */
/* and facilitate your testing of the programs, we wrote a short */
/* simple main program for each program so that you can readily test */
/* them. */
/*      All the programs are written in FORTRAN-77 and tested on PCs */
/* and workstations. Therefore, they should run on any computer with */
/* implementation of the FORTRAN-77 standard. */
/*      Although we have made a great effort to test these programs, */
/* we would not be surprised  to find some errors in them.  We would */
/* appreciate it if you can bring to our attention any errors you find. */
/* You can do this by either writing us directly at the location */
/* (e-mail: j-jin1@uiuc.edu) or writing to the publisher, whose address */
/* appears on the back cover of the book.  However, we must note that */
/* all these programs are sold "as is," and we cannot guarantee to */
/* correct the errors reported by readers on any fixed schedule. */
/*      All the programs and subroutines  contained in this diskette */
/* are copyrighted.   However,  we give permission to the reader who */
/* purchases this book to incorporate any of these programs into his */
/* or her programs provided that the copyright is acknowledged. */

/*       =================================================== */
/*       Purpose: This program computes the gamma function */
/*                ・x) for x > 0 using subroutine LGAMA */
/*       Examples: */
/*                  x           ・x) */
/*                ------------------------- */
/*                 0.5     .1772453851D+01 */
/*                 2.5     .1329340388D+01 */
/*                 5.0     .2400000000D+02 */
/*                 7.5     .1871254306D+04 */
/*                10.0     .3628800000D+06 */
/*       =================================================== */

/* Subroutine */ int lgama_(integer *kf, doublereal *x, doublereal *gl)
{
    /* Initialized data */

    static doublereal a[10] = { .08333333333333333,-.002777777777777778,
	    7.936507936507937e-4,-5.952380952380952e-4,8.417508417508418e-4,
	    -.001917526917526918,.00641025641025641,-.02955065359477124,
	    .1796443723688307,-1.3924322169059 };

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer k, n;
    static doublereal x0, x2, xp, gl0;


/*       ================================================== */
/*       Purpose: Compute gamma function ・x) or ln[・x)] */
/*       Input:   x  --- Argument of ・x) ( x > 0 ) */
/*                KF --- Function code */
/*                       KF=1 for ・x); KF=0 for ln[・x)] */
/*       Output:  GL --- ・x) or ln[・x)] */
/*       ================================================== */

    x0 = *x;
    if (*x == 1.f || *x == 2.f) {
	*gl = 0.;
	goto L20;
    } else if (*x <= 7.f) {
	n = (integer) (7 - *x);
	x0 = *x + n;
    }
    x2 = 1. / (x0 * x0);
    xp = 6.283185307179586477;
    gl0 = a[9];
    for (k = 9; k >= 1; --k) {
/* L10: */
	gl0 = gl0 * x2 + a[k - 1];
    }
    *gl = gl0 / x0 + log(xp) * .5 + (x0 - .5) * log(x0) - x0;
    if (*x <= 7.f) {
	i__1 = n;
	for (k = 1; k <= i__1; ++k) {
	    *gl -= log(x0 - 1.);
/* L15: */
	    x0 += -1.;
	}
    }
L20:
    if (*kf == 1) {
	*gl = exp(*gl);
    }
    return 0;
} /* lgama_ */

/*  This is William Cody's erf implementations */
/* Subroutine */ int calerf_(doublereal *arg, doublereal *result, integer *
	jint)
{
    /* Initialized data */

    static doublereal four = 4.;
    static doublereal one = 1.;
    static doublereal half = .5;
    static doublereal two = 2.;
    static doublereal zero = 0.;
    static doublereal sqrpi = .56418958354775628695;
    static doublereal thresh = .46875;
    static doublereal sixten = 16.;
    static doublereal xinf = 1.79e308;
    static doublereal xneg = -26.628;
    static doublereal xsmall = 1.11e-16;
    static doublereal xbig = 26.543;
    static doublereal xhuge = 6.71e7;
    static doublereal xmax = 2.53e307;
    static doublereal a[5] = { 3.1611237438705656,113.864154151050156,
	    377.485237685302021,3209.37758913846947,.185777706184603153 };
    static doublereal b[4] = { 23.6012909523441209,244.024637934444173,
	    1282.61652607737228,2844.23683343917062 };
    static doublereal c__[9] = { .564188496988670089,8.88314979438837594,
	    66.1191906371416295,298.635138197400131,881.95222124176909,
	    1712.04761263407058,2051.07837782607147,1230.33935479799725,
	    2.15311535474403846e-8 };
    static doublereal d__[8] = { 15.7449261107098347,117.693950891312499,
	    537.181101862009858,1621.38957456669019,3290.79923573345963,
	    4362.61909014324716,3439.36767414372164,1230.33935480374942 };
    static doublereal p[6] = { .305326634961232344,.360344899949804439,
	    .125781726111229246,.0160837851487422766,6.58749161529837803e-4,
	    .0163153871373020978 };
    static doublereal q[5] = { 2.56852019228982242,1.87295284992346047,
	    .527905102951428412,.0605183413124413191,.00233520497626869185 };

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), d_int(doublereal *);

    /* Local variables */
    static integer i__;
    static doublereal x, y, del, ysq, xden, xnum;

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

/*   Approximate values for some important machines are: */

/*                          XMIN       XINF        XNEG     XSMALL */

/*  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15 */
/*  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15 */
/*  IEEE (IBM/XT, */
/*    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8 */
/*  IEEE (IBM/XT, */
/*    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16 */
/*  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17 */
/*  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18 */
/*  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17 */
/*  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16 */


/*                          XBIG       XHUGE       XMAX */

/*  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293 */
/*  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465 */
/*  IEEE (IBM/XT, */
/*    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37 */
/*  IEEE (IBM/XT, */
/*    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307 */
/*  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75 */
/*  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307 */
/*  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38 */
/*  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307 */

/* ******************************************************************* */
/* ******************************************************************* */

/* Error returns */

/*  The program returns  ERFC = 0      for  ARG .GE. XBIG; */

/*                       ERFCX = XINF  for  ARG .LT. XNEG; */
/*      and */
/*                       ERFCX = 0     for  ARG .GE. XMAX. */


/* Intrinsic functions required are: */

/*     ABS, AINT, EXP */


/*  Author: W. J. Cody */
/*          Mathematics and Computer Science Division */
/*          Argonne National Laboratory */
/*          Argonne, IL 60439 */

/*  Latest modification: March 19, 1990 */

/* ------------------------------------------------------------------ */
/* S    REAL */
/* ------------------------------------------------------------------ */
/*  Mathematical constants */
/* ------------------------------------------------------------------ */
/* S    DATA FOUR,ONE,HALF,TWO,ZERO/4.0E0,1.0E0,0.5E0,2.0E0,0.0E0/, */
/* S   1     SQRPI/5.6418958354775628695E-1/,THRESH/0.46875E0/, */
/* S   2     SIXTEN/16.0E0/ */
/* ------------------------------------------------------------------ */
/*  Machine-dependent constants */
/* ------------------------------------------------------------------ */
/* S    DATA XINF,XNEG,XSMALL/3.40E+38,-9.382E0,5.96E-8/, */
/* S   1     XBIG,XHUGE,XMAX/9.194E0,2.90E3,4.79E37/ */
/* ------------------------------------------------------------------ */
/*  Coefficients for approximation to  erf  in first interval */
/* ------------------------------------------------------------------ */
/* S    DATA A/3.16112374387056560E00,1.13864154151050156E02, */
/* S   1       3.77485237685302021E02,3.20937758913846947E03, */
/* S   2       1.85777706184603153E-1/ */
/* S    DATA B/2.36012909523441209E01,2.44024637934444173E02, */
/* S   1       1.28261652607737228E03,2.84423683343917062E03/ */
/* ------------------------------------------------------------------ */
/*  Coefficients for approximation to  erfc  in second interval */
/* ------------------------------------------------------------------ */
/* S    DATA C/5.64188496988670089E-1,8.88314979438837594E0, */
/* S   1       6.61191906371416295E01,2.98635138197400131E02, */
/* S   2       8.81952221241769090E02,1.71204761263407058E03, */
/* S   3       2.05107837782607147E03,1.23033935479799725E03, */
/* S   4       2.15311535474403846E-8/ */
/* S    DATA D/1.57449261107098347E01,1.17693950891312499E02, */
/* S   1       5.37181101862009858E02,1.62138957456669019E03, */
/* S   2       3.29079923573345963E03,4.36261909014324716E03, */
/* S   3       3.43936767414372164E03,1.23033935480374942E03/ */
/* ------------------------------------------------------------------ */
/*  Coefficients for approximation to  erfc  in third interval */
/* ------------------------------------------------------------------ */
/* S    DATA P/3.05326634961232344E-1,3.60344899949804439E-1, */
/* S   1       1.25781726111229246E-1,1.60837851487422766E-2, */
/* S   2       6.58749161529837803E-4,1.63153871373020978E-2/ */
/* S    DATA Q/2.56852019228982242E00,1.87295284992346047E00, */
/* S   1       5.27905102951428412E-1,6.05183413124413191E-2, */
/* S   2       2.33520497626869185E-3/ */
/* ------------------------------------------------------------------ */
    x = *arg;
    y = abs(x);
    if (y <= thresh) {
/* ------------------------------------------------------------------ */
/*  Evaluate  erf  for  |X| <= 0.46875 */
/* ------------------------------------------------------------------ */
	ysq = zero;
	if (y > xsmall) {
	    ysq = y * y;
	}
	xnum = a[4] * ysq;
	xden = ysq;
	for (i__ = 1; i__ <= 3; ++i__) {
	    xnum = (xnum + a[i__ - 1]) * ysq;
	    xden = (xden + b[i__ - 1]) * ysq;
/* L20: */
	}
	*result = x * (xnum + a[3]) / (xden + b[3]);
	if (*jint != 0) {
	    *result = one - *result;
	}
	if (*jint == 2) {
	    *result = exp(ysq) * *result;
	}
	goto L800;
/* ------------------------------------------------------------------ */
/*  Evaluate  erfc  for 0.46875 <= |X| <= 4.0 */
/* ------------------------------------------------------------------ */
    } else if (y <= four) {
	xnum = c__[8] * y;
	xden = y;
	for (i__ = 1; i__ <= 7; ++i__) {
	    xnum = (xnum + c__[i__ - 1]) * y;
	    xden = (xden + d__[i__ - 1]) * y;
/* L120: */
	}
	*result = (xnum + c__[7]) / (xden + d__[7]);
	if (*jint != 2) {
	    d__1 = y * sixten;
	    ysq = d_int(&d__1) / sixten;
	    del = (y - ysq) * (y + ysq);
	    *result = exp(-ysq * ysq) * exp(-del) * *result;
	}
/* ------------------------------------------------------------------ */
/*  Evaluate  erfc  for |X| > 4.0 */
/* ------------------------------------------------------------------ */
    } else {
	*result = zero;
	if (y >= xbig) {
	    if (*jint != 2 || y >= xmax) {
		goto L300;
	    }
	    if (y >= xhuge) {
		*result = sqrpi / y;
		goto L300;
	    }
	}
	ysq = one / (y * y);
	xnum = p[5] * ysq;
	xden = ysq;
	for (i__ = 1; i__ <= 4; ++i__) {
	    xnum = (xnum + p[i__ - 1]) * ysq;
	    xden = (xden + q[i__ - 1]) * ysq;
/* L240: */
	}
	*result = ysq * (xnum + p[4]) / (xden + q[4]);
	*result = (sqrpi - *result) / y;
	if (*jint != 2) {
	    d__1 = y * sixten;
	    ysq = d_int(&d__1) / sixten;
	    del = (y - ysq) * (y + ysq);
	    *result = exp(-ysq * ysq) * exp(-del) * *result;
	}
    }
/* ------------------------------------------------------------------ */
/*  Fix up for negative argument, erf, etc. */
/* ------------------------------------------------------------------ */
L300:
    if (*jint == 0) {
	*result = half - *result + half;
	if (x < zero) {
	    *result = -(*result);
	}
    } else if (*jint == 1) {
	if (x < zero) {
	    *result = two - *result;
	}
    } else {
	if (x < zero) {
	    if (x < xneg) {
		*result = xinf;
	    } else {
		d__1 = x * sixten;
		ysq = d_int(&d__1) / sixten;
		del = (x - ysq) * (x + ysq);
		y = exp(ysq * ysq) * exp(del);
		*result = y + y - *result;
	    }
	}
    }
L800:
    return 0;
/* ---------- Last card of CALERF ---------- */
} /* calerf_ */

/* S    REAL FUNCTION ERF(X) */
doublereal derf_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer jint;
    extern /* Subroutine */ int calerf_(doublereal *, doublereal *, integer *)
	    ;
    static doublereal result;

/* -------------------------------------------------------------------- */

/* This subprogram computes approximate values for erf(x). */
/*   (see comments heading CALERF). */

/*   Author/date: W. J. Cody, January 8, 1985 */

/* -------------------------------------------------------------------- */
/* S    REAL             X, RESULT */
/* ------------------------------------------------------------------ */
    jint = 0;
    calerf_(x, &result, &jint);
/* S    ERF = RESULT */
    ret_val = result;
    return ret_val;
/* ---------- Last card of DERF ---------- */
} /* derf_ */

/* S    REAL FUNCTION ERFC(X) */
doublereal derfc_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer jint;
    extern /* Subroutine */ int calerf_(doublereal *, doublereal *, integer *)
	    ;
    static doublereal result;

/* -------------------------------------------------------------------- */

/* This subprogram computes approximate values for erfc(x). */
/*   (see comments heading CALERF). */

/*   Author/date: W. J. Cody, January 8, 1985 */

/* -------------------------------------------------------------------- */
/* S    REAL             X, RESULT */
/* ------------------------------------------------------------------ */
    jint = 1;
    calerf_(x, &result, &jint);
/* S    ERFC = RESULT */
    ret_val = result;
    return ret_val;
/* ---------- Last card of DERFC ---------- */
} /* derfc_ */

/* S    REAL FUNCTION ERFCX(X) */
doublereal derfcx_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer jint;
    extern /* Subroutine */ int calerf_(doublereal *, doublereal *, integer *)
	    ;
    static doublereal result;

/* ------------------------------------------------------------------ */

/* This subprogram computes approximate values for exp(x*x) * erfc(x). */
/*   (see comments heading CALERF). */

/*   Author/date: W. J. Cody, March 30, 1987 */

/* ------------------------------------------------------------------ */
/* S    REAL             X, RESULT */
/* ------------------------------------------------------------------ */
    jint = 2;
    calerf_(x, &result, &jint);
/* S    ERFCX = RESULT */
    ret_val = result;
    return ret_val;
/* ---------- Last card of DERFCX ---------- */
} /* derfcx_ */

/*     ALGORITHM 488 COLLECTED ALGORITHMS FROM ACM. */
/*     ALGORITHM APPEARED IN COMM. ACM, VOL. 17, NO. 12, */
/*     P. 704. */
doublereal dgrand_(integer *n)
{
    /* Initialized data */

    static doublereal d__[60] = { .67448975,.47585963,.383771164,.328611323,
	    .291142827,.263684322,.242508452,.225567444,.211634166,.199924267,
	    .189910758,.181225181,.1736014,.166841909,.160796729,.155349717,
	    .150409384,.145902577,.141770033,.137963174,.134441762,.13117215,
	    .128125965,.12527909,.122610883,.12010356,.117741707,.115511892,
	    .113402349,.11140272,.109503852,.107697617,.105976772,.104334841,
	    .102766012,.101265052,.099827234,.098448282,.097124309,.095851778,
	    .094627461,.093448407,.092311909,.091215482,.090156838,.089133867,
	    .088144619,.087187293,.086260215,.085361834,.084490706,.083645487,
	    .082824924,.082027847,.081253162,.080499844,.079766932,.079053527,
	    .078358781,.077681899 };
    static doublereal u = 0.f;

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal a;
    static integer i__;
    static doublereal v, w;
    extern doublereal rand_(integer *);

/* EXCEPT ON THE FIRST CALL GRAND RETURNS A */
/* PSEUDO-RANDOM NUMBER HAVING A GAUSSIAN (I.E. */
/* NORMAL) DISTRIBUTION WITH ZERO MEAN AND UNIT */
/* STANDARD DEVIATION.  THUS, THE DENSITY IS  F(X) = */
/* EXP(-0.5*X**2)/SQRT(2.0*PI). THE FIRST CALL */
/* INITIALIZES GRAND AND RETURNS ZERO. */
/* THE PARAMETER N IS DUMMY. */
/* GRAND CALLS A FUNCTION RAND, AND IT IS ASSUMED THAT */
/* SUCCESSIVE CALLS TO RAND(0) GIVE INDEPENDENT */
/* PSEUDO- RANDOM NUMBERS DISTRIBUTED UNIFORMLY ON (0, */
/* 1), POSSIBLY INCLUDING 0 (BUT NOT 1). */
/* THE METHOD USED WAS SUGGESTED BY VON NEUMANN, AND */
/* IMPROVED BY FORSYTHE, AHRENS, DIETER AND BRENT. */
/* ON THE AVERAGE THERE ARE 1.37746 CALLS OF RAND FOR */
/* EACH CALL OF GRAND. */
/* WARNING - DIMENSION AND DATA STATEMENTS BELOW ARE */
/*           MACHINE-DEPENDENT. */
/* DIMENSION OF D MUST BE AT LEAST THE NUMBER OF BITS */
/* IN THE FRACTION OF A FLOATING-POINT NUMBER. */
/* THUS, ON MOST MACHINES THE DATA STATEMENT BELOW */
/* CAN BE TRUNCATED. */
/* IF THE INTEGRAL OF SQRT(2.0/PI)*EXP(-0.5*X**2) FROM */
/* A(I) TO INFINITY IS 2**(-I), THEN D(I) = A(I) - */
/* A(I-1). */
/* END OF MACHINE-DEPENDENT STATEMENTS */
/* U MUST BE PRESERVED BETWEEN CALLS. */
/* INITIALIZE DISPLACEMENT A AND COUNTER I. */
    a = 0.f;
    i__ = 0;
/* INCREMENT COUNTER AND DISPLACEMENT IF LEADING BIT */
/* OF U IS ONE. */
L10:
    u += u;
    if (u < 1.f) {
	goto L20;
    }
    u += -1.f;
    ++i__;
    a -= d__[i__ - 1];
    goto L10;
/* FORM W UNIFORM ON 0 .LE. W .LT. D(I+1) FROM U. */
L20:
    w = d__[i__] * u;
/* FORM V = 0.5*((W-A)**2 - A**2). NOTE THAT 0 .LE. V */
/* .LT. LOG(2). */
    v = w * (w * .5f - a);
/* GENERATE NEW UNIFORM U. */
L30:
    u = rand_(&c__0);
/* ACCEPT W AS A RANDOM SAMPLE IF V .LE. U. */
    if (v <= u) {
	goto L40;
    }
/* GENERATE RANDOM V. */
    v = rand_(&c__0);
/* LOOP IF U .GT. V. */
    if (u > v) {
	goto L30;
    }
/* REJECT W AND FORM A NEW UNIFORM U FROM V AND U. */
    u = (v - u) / (1.f - u);
    goto L20;
/* FORM NEW U (TO BE USED ON NEXT CALL) FROM U AND V. */
L40:
    u = (u - v) / (1.f - v);
/* USE FIRST BIT OF U FOR SIGN, RETURN NORMAL VARIATE. */
    u += u;
    if (u < 1.f) {
	goto L50;
    }
    u += -1.f;
    ret_val = w - a;
    return ret_val;
L50:
    ret_val = a - w;
    return ret_val;
} /* dgrand_ */

/* Algorithm 467 from ACM (can't get this working at the moment). */
/* See http://pdp-10.trailing-edge.com/red405a2/11/uetp/lib/467.for */
/* Subroutine */ int dxpose_(doublereal *a, integer *n1, integer *n2, integer 
	*n12, logical *moved, integer *nwork)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, m, n, i1, i2, ip, ia1, ia2, idiv, iexp[8], nexp[8], 
	    mmia1, mmia2, i1min, i1max, ifact[8];
    static doublereal atemp, btemp;
    static integer isoid, itest, mmist;
    extern /* Subroutine */ int factor_(integer *, integer *, integer *, 
	    integer *, integer *);
    static integer ipower[8], ncount, istart, npower;

/*  TRANSPOSITION OF A RECTANGULAR MATRIX IN SITU. */
/*  BY NORMAN BRENNER, MIT, 1/72.  CF. ALG. 380, CACM, 5/70. */
/*  TRANSPOSITION OF THE N1 BY N2 MATRIX A AMOUNTS TO */
/*  REPLACING THE ELEMENT AT VECTOR POSITION I (0-ORIGIN) */
/*  WITH THE ELEMENT AT POSITION N1*I (MOD N1*N2-I). */
/*  EACH SUBCYCLE OF THIS PERMUTATION IS COMPLETED IN ORDER. */
/*  MOVED IS A LOGICAL WORK ARRAY OF LENGTH NWORK. */
/*  REALLY A(N1,N2), BUT N12 = N1*N2 */
    /* Parameter adjustments */
    --a;
    --moved;

    /* Function Body */
    if (*n1 < 2 || *n2 < 2) {
	return 0;
    }
    n = *n1;
    m = *n1 * *n2 - 1;
    if (*n1 != *n2) {
	goto L30;
    }
/*  SQUARE MATRICES ARE DONE SEPARATELY FOR SPEED */
    i1min = 2;
    i__1 = m;
    i__2 = n;
    for (i1max = n; i__2 < 0 ? i1max >= i__1 : i1max <= i__1; i1max += i__2) {
	i2 = i1min + n - 1;
	i__3 = i1max;
	for (i1 = i1min; i1 <= i__3; ++i1) {
	    atemp = a[i1];
	    a[i1] = a[i2];
	    a[i2] = atemp;
	    i2 += n;
/* L10: */
	}
	i1min = i1min + n + 1;
/* L20: */
    }
    return 0;
/*  MODULUS M IS FACTORED INTO PRIME POWERS.  EIGHT FACTORS */
/*  SUFFICE UP TO M = 2*3*5*7*11*13*17*19 = 9,767,520. */
L30:
    factor_(&m, ifact, ipower, nexp, &npower);
    i__2 = npower;
    for (ip = 1; ip <= i__2; ++ip) {
	iexp[ip - 1] = 0;
/* L40: */
    }
/*  GENERATE EVERY DIVISOR OF M LESS THAN M/2 */
    idiv = 1;
L50:
    if (idiv >= m / 2) {
	goto L190;
    }
/*  THE NUMBER OF ELEMENTS WHOSE INDEX IS DIVISIBLE BY IDIV */
/*  AND BY NO OTHER DIVISOR OF M IS THE EULER TOTIENT */
/*  FUNCTION, PHI(M/IDIV). */
    ncount = m / idiv;
    i__2 = npower;
    for (ip = 1; ip <= i__2; ++ip) {
	if (iexp[ip - 1] == nexp[ip - 1]) {
	    goto L60;
	}
	ncount = ncount / ifact[ip - 1] * (ifact[ip - 1] - 1);
L60:
	;
    }
    i__2 = *nwork;
    for (i__ = 1; i__ <= i__2; ++i__) {
	moved[i__] = FALSE_;
/* L70: */
    }
/*  THE STARTING POINT OF A SUBCYCLE IS DIVISIBLE ONLY BY IDIV */
/*  AND MUST NOT APPEAR IN ANY OTHER SUBCYCLE. */
    istart = idiv;
L80:
    mmist = m - istart;
    if (istart == idiv) {
	goto L120;
    }
    if (istart > *nwork) {
	goto L90;
    }
    if (moved[istart]) {
	goto L160;
    }
L90:
    isoid = istart / idiv;
    i__2 = npower;
    for (ip = 1; ip <= i__2; ++ip) {
	if (iexp[ip - 1] == nexp[ip - 1]) {
	    goto L100;
	}
	if (isoid % ifact[ip - 1] == 0) {
	    goto L160;
	}
L100:
	;
    }
    if (istart <= *nwork) {
	goto L120;
    }
    itest = istart;
L110:
    itest = n * itest % m;
    if (itest < istart || itest > mmist) {
	goto L160;
    }
    if (itest > istart && itest < mmist) {
	goto L110;
    }
L120:
    atemp = a[istart + 1];
    btemp = a[mmist + i__];
    ia1 = istart;
L130:
    ia2 = n * ia1 % m;
    mmia1 = m - ia1;
    mmia2 = m - ia2;
    if (ia1 <= *nwork) {
	moved[ia1] = TRUE_;
    }
    if (mmia1 <= *nwork) {
	moved[mmia1] = TRUE_;
    }
    ncount += -2;
/*  MOVE TWO ELEMENTS, THE SECOND FROM THE NEGATIVE */
/*  SUBCYCLE.  CHECK FIRST FOR SUBCYCLE CLOSURE. */
    if (ia2 == istart) {
	goto L140;
    }
    if (mmia2 == istart) {
	goto L150;
    }
    a[ia1 + 1] = a[ia2 + 1];
    a[mmia1 + 1] = a[mmia2 + 1];
    ia1 = ia2;
    goto L130;
L140:
    a[ia1 + 1] = atemp;
    a[mmia1 + 1] = btemp;
    goto L160;
L150:
    a[ia1 + 1] = btemp;
    a[mmia1 + 1] = atemp;
L160:
    istart += idiv;
    if (ncount > 0) {
	goto L80;
    }
    i__2 = npower;
    for (ip = 1; ip <= i__2; ++ip) {
	if (iexp[ip - 1] == nexp[ip - 1]) {
	    goto L170;
	}
	++iexp[ip - 1];
	idiv *= ifact[ip - 1];
	goto L50;
L170:
	iexp[ip - 1] = 0;
	idiv /= ipower[ip - 1];
/* L180: */
    }
L190:
    return 0;
} /* dxpose_ */

/* Subroutine */ int factor_(integer *n, integer *ifact, integer *ipower, 
	integer *nexp, integer *npower)
{
    static integer ip, idiv, ifcur, npart, iquot;

/*  FACTOR N INTO ITS PRIME POWERS, NPOWER IN NUMBER. */
/*  E.G., FOR N=1970=2**3 *5 *7**2, NPOWER=3, IFACT=3,5,7, */
/*  IPOWER=8,5,49, AND NEXP=3,1,2. */
    /* Parameter adjustments */
    --nexp;
    --ipower;
    --ifact;

    /* Function Body */
    ip = 0;
    ifcur = 0;
    npart = *n;
    idiv = 2;
L10:
    iquot = npart / idiv;
    if (npart - idiv * iquot != 0) {
	goto L60;
    } else {
	goto L20;
    }
L20:
    if (idiv - ifcur <= 0) {
	goto L40;
    } else {
	goto L30;
    }
L30:
    ++ip;
    ifact[ip] = idiv;
    ipower[ip] = idiv;
    ifcur = idiv;
    nexp[ip] = 1;
    goto L50;
L40:
    ipower[ip] = idiv * ipower[ip];
    ++nexp[ip];
L50:
    npart = iquot;
    goto L10;
L60:
    if (iquot - idiv <= 0) {
	goto L100;
    } else {
	goto L70;
    }
L70:
    if (idiv - 2 <= 0) {
	goto L80;
    } else {
	goto L90;
    }
L80:
    idiv = 3;
    goto L10;
L90:
    idiv += 2;
    goto L10;
L100:
    if (npart - 1 <= 0) {
	goto L140;
    } else {
	goto L110;
    }
L110:
    if (npart - ifcur <= 0) {
	goto L130;
    } else {
	goto L120;
    }
L120:
    ++ip;
    ifact[ip] = npart;
    ipower[ip] = npart;
    nexp[ip] = 1;
    goto L140;
L130:
    ipower[ip] = npart * ipower[ip];
    ++nexp[ip];
L140:
    *npower = ip;
    return 0;
} /* factor_ */

/* ACM Algorithm 380 */
/* Subroutine */ int dtrans_(doublereal *a, integer *m, integer *n, integer *
	mn, integer *move, integer *iwrk, integer *iok)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal b;
    static integer i__, j, k, i1, i2, j1, m2, n1, ia, ib, kmi, max__, ncount;

/* A IS A ONE-DIMENSIONAL ARRAY OF LENGTH MN=M*N, WHICH */
/* CONTAINS THE M BY N MATRIX TO BE TRANSPOSED (STORED */
/* COLUMNWISE).MOVE IS A ONE-DIMENSIONAL ARRAY OF LENGTH IWRK */
/* USED TO STORE INFORMATION TO SPEED UP THE PROCESS. THE */
/* VALUE IWRK=(M+N)/2 IS RECOMMENDED. IOK INDICATES THE */
/* SUCCESS OR FAILURE OF THE ROUTINE. */
/* NORMAL RETURN IOK=0 */
/* ERRORS           IOK=-1, MN NOT EQUAL TO M*N. */
/*                  IOK=-2, IWRK NEGATIVE OR ZERO. */
/*                  IOK.GT.0, (SHOULD NEVER OCCUR). IN THIS CASE */
/* WE SET IOK EQUAL TO THE FINAL VALUE OF I WHEN THE SEARCH */
/* IS COMPLETED BUT SOME LOOPS HAVE NOT BEEN MOVED. */
/* CHECK ARGUMENTS AND INITIALISE */
    /* Parameter adjustments */
    --a;
    --move;

    /* Function Body */
    if (*m < 2 || *n < 2) {
	goto L60;
    }
    if (*mn != *m * *n) {
	goto L92;
    }
    if (*iwrk < 1) {
	goto L93;
    }
    if (*m == *n) {
	goto L70;
    }
    ncount = 2;
    m2 = *m - 2;
    i__1 = *iwrk;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	move[i__] = 0;
    }
    if (m2 < 1) {
	goto L12;
    }
/* COUNT NUMBER,NCOUNT, OF SINGLE POINTS. */
    i__1 = m2;
    for (ia = 1; ia <= i__1; ++ia) {
	ib = ia * (*n - 1) / (*m - 1);
	if (ia * (*n - 1) != ib * (*m - 1)) {
	    goto L11;
	}
	++ncount;
	i__ = ia * *n + ib;
	if (i__ > *iwrk) {
	    goto L11;
	}
	move[i__] = 1;
L11:
	;
    }
/* SET INITIAL VALUES FOR SEARCH. */
L12:
    k = *mn - 1;
    kmi = k - 1;
    max__ = *mn;
    i__ = 1;
/* AT LEAST ONE LOOP MUST BE RE-ARRANGED. */
    goto L30;
/* SEARCH FOR LOOPS TO BE REARRANGED. */
L20:
    max__ = k - i__;
    ++i__;
    kmi = k - i__;
    if (i__ > max__) {
	goto L90;
    }
    if (i__ > *iwrk) {
	goto L21;
    }
    if (move[i__] < 1) {
	goto L30;
    }
    goto L20;
L21:
    if (i__ == *m * i__ - k * (i__ / *n)) {
	goto L20;
    }
    i1 = i__;
L22:
    i2 = *m * i1 - k * (i1 / *n);
    if (i2 <= i__ || i2 >= max__) {
	goto L23;
    }
    i1 = i2;
    goto L22;
L23:
    if (i2 != i__) {
	goto L20;
    }
/* REARRANGE ELEMENTS OF A LOOP. */
L30:
    i1 = i__;
L31:
    b = a[i1 + 1];
L32:
    i2 = *m * i1 - k * (i1 / *n);
    if (i1 <= *iwrk) {
	move[i1] = 2;
    }
/* L33: */
    ++ncount;
    if (i2 == i__ || i2 >= kmi) {
	goto L35;
    }
L34:
    a[i1 + 1] = a[i2 + 1];
    i1 = i2;
    goto L32;
L35:
    if (max__ == kmi || i2 == i__) {
	goto L41;
    }
    max__ = kmi;
    goto L34;
/* TEST FOR SYMMETRIC PAIR OF LOOPS. */
L41:
    a[i1 + 1] = b;
    if (ncount >= *mn) {
	goto L60;
    }
    if (i2 == max__ || max__ == kmi) {
	goto L20;
    }
    max__ = kmi;
    i1 = max__;
    goto L31;
/* NORMAL RETURN. */
L60:
    *iok = 0;
    return 0;
/* IF MATRIX IS SQUARE, EXCHANGE ELEMENTS A(I,J) AND A(J,I). */
L70:
    n1 = *n - 1;
    i__1 = n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j1 = i__ + 1;
	i__2 = *n;
	for (j = j1; j <= i__2; ++j) {
	    i1 = i__ + (j - 1) * *n;
	    i2 = j + (i__ - 1) * *m;
	    b = a[i1];
	    a[i1] = a[i2];
	    a[i2] = b;
/* L71: */
	}
    }
    goto L60;
/* ERROR RETURNS. */
L90:
    *iok = i__;
L91:
    return 0;
L92:
    *iok = -1;
    goto L91;
L93:
    *iok = -2;
    goto L91;
} /* dtrans_ */

/* ACM Algorithm 380. */
/* Subroutine */ int dtransr_(doublereal *a, integer *m, integer *n, integer *
	mn, integer *move, integer *iwrk, integer *iok)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal b, c__, d__;
    static integer i__, j, k, i1, i2, j1, n1, im, i1c, i2c, ir0, ir1, ir2, 
	    kmi, max__, ncount;

/* ***** */
/*  ALGORITHM 380 - REVISED */
/* ***** */
/*  A IS A ONE-DIMENSIONAL ARRAY OF LENGTH MN=M*N, WHICH */
/*  CONTAINS THE MXN MATRIX TO BE TRANSPOSED (STORED */
/*  COLUMWISE). MOVE IS A ONE-DIMENSIONAL ARRAY OF LENGTH IWRK */
/*  USED TO STORE INFORMATION TO SPEED UP THE PROCESS.  THE */
/*  VALUE IWRK=(M+N)/2 IS RECOMMENDED. IOK INDICATES THE */
/*  SUCCESS OR FAILURE OF THE ROUTINE. */
/*  NORMAL RETURN  IOK=0 */
/*  ERRORS         IOK=-1 ,MN NOT EQUAL TO M*N */
/*                 IOK=-2 ,IWRK NEGATIVE OR ZERO */
/*                 IOK.GT.0, (SHOULD NEVER OCCUR),IN THIS CASE */
/*  WE SET IOK EQUAL TO THE FINAL VALUE OF I WHEN THE SEARCH */
/*  IS COMPLETED BUT SOME LOOPS HAVE NOT BEEN MOVED */
/*  NOTE * MOVE(I) WILL STAY ZERO FOR FIXED POINTS */
/* CHECK ARGUMENTS AND INITIALIZE. */
    /* Parameter adjustments */
    --a;
    --move;

    /* Function Body */
    if (*m < 2 || *n < 2) {
	goto L120;
    }
    if (*mn != *m * *n) {
	goto L180;
    }
    if (*iwrk < 1) {
	goto L190;
    }
    if (*m == *n) {
	goto L130;
    }
    ncount = 2;
    k = *mn - 1;
    i__1 = *iwrk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	move[i__] = 0;
/* L10: */
    }
    if (*m < 3 || *n < 3) {
	goto L30;
    }
/* CALCULATE THE NUMBER OF FIXED POINTS, EUCLIDS ALGORITHM */
/* FOR GCD(M-1,N-1). */
    ir2 = *m - 1;
    ir1 = *n - 1;
L20:
    ir0 = ir2 % ir1;
    ir2 = ir1;
    ir1 = ir0;
    if (ir0 != 0) {
	goto L20;
    }
    ncount = ncount + ir2 - 1;
/* SET INITIAL VALUES FOR SEARCH */
L30:
    i__ = 1;
    im = *m;
/* AT LEAST ONE LOOP MUST BE RE-ARRANGED */
    goto L80;
/* SEARCH FOR LOOPS TO REARRANGE */
L40:
    max__ = k - i__;
    ++i__;
    if (i__ > max__) {
	goto L160;
    }
    im += *m;
    if (im > k) {
	im -= k;
    }
    i2 = im;
    if (i__ == i2) {
	goto L40;
    }
    if (i__ > *iwrk) {
	goto L60;
    }
    if (move[i__] == 0) {
	goto L80;
    }
    goto L40;
L50:
    i2 = *m * i1 - k * (i1 / *n);
L60:
    if (i2 <= i__ || i2 >= max__) {
	goto L70;
    }
    i1 = i2;
    goto L50;
L70:
    if (i2 != i__) {
	goto L40;
    }
/* REARRANGE THE ELEMENTS OF A LOOP AND ITS COMPANION LOOP */
L80:
    i1 = i__;
    kmi = k - i__;
    b = a[i1 + 1];
    i1c = kmi;
    c__ = a[i1c + 1];
L90:
    i2 = *m * i1 - k * (i1 / *n);
    i2c = k - i2;
    if (i1 <= *iwrk) {
	move[i1] = 2;
    }
    if (i1c <= *iwrk) {
	move[i1c] = 2;
    }
    ncount += 2;
    if (i2 == i__) {
	goto L110;
    }
    if (i2 == kmi) {
	goto L100;
    }
    a[i1 + 1] = a[i2 + 1];
    a[i1c + 1] = a[i2c + 1];
    i1 = i2;
    i1c = i2c;
    goto L90;
/* FINAL STORE AND TEST FOR FINISHED */
L100:
    d__ = b;
    b = c__;
    c__ = d__;
L110:
    a[i1 + 1] = b;
    a[i1c + 1] = c__;
    if (ncount < *mn) {
	goto L40;
    }
/* NORMAL RETURN */
L120:
    *iok = 0;
    return 0;
/* IF MATRIX IS SQUARE,EXCHANGE ELEMENTS A(I,J) AND A(J,I). */
L130:
    n1 = *n - 1;
    i__1 = n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j1 = i__ + 1;
	i__2 = *n;
	for (j = j1; j <= i__2; ++j) {
	    i1 = i__ + (j - 1) * *n;
	    i2 = j + (i__ - 1) * *m;
	    b = a[i1];
	    a[i1] = a[i2];
	    a[i2] = b;
/* L140: */
	}
/* L150: */
    }
    goto L120;
/* ERROR RETURNS. */
L160:
    *iok = i__;
L170:
    return 0;
L180:
    *iok = -1;
    goto L170;
L190:
    *iok = -2;
    goto L170;
} /* dtransr_ */

