C http://jin.ece.uiuc.edu/routines/routines.html
C 		****************************************
C 		*           DISK TO ACCOMPANY          *
C 		*   COMPUTATION OF SPECIAL FUNCTIONS   *
C 		*                                      *
C 		*   Shanjie Zhang and Jianming Jin     *
C 		*                                      *
C 		*   Copyright 1996 by John Wiley &     *
C 		*              Sons, Inc.              *
C 		*                                      *
C 		****************************************
C
C I. INTRODUCTION
C
C      As stated in the preface of our book "Computation of Special 
C Functions,"  the purpose of this book is to share with the reader  
C a set of computer programs (130 in total) which we have developed 
C during the past several years for computing a variety of  special  
C mathematical functions.  For your convenience, we attach to the
C book this diskette that contains all the computer programs  
C listed or mentioned in the book. 

C      In this diskette,  we place all the programs under directory 
C SMF\PROGRAMS. In order to illustrate the use of these programs 
C and facilitate your testing of the programs, we wrote a short 
C simple main program for each program so that you can readily test 
C them.    

C      All the programs are written in FORTRAN-77 and tested on PCs
C and workstations. Therefore, they should run on any computer with 
C implementation of the FORTRAN-77 standard.  

C      Although we have made a great effort to test these programs,  
C we would not be surprised  to find some errors in them.  We would 
C appreciate it if you can bring to our attention any errors you find.
C You can do this by either writing us directly at the location
C (e-mail: j-jin1@uiuc.edu) or writing to the publisher, whose address 
C appears on the back cover of the book.  However, we must note that
C all these programs are sold "as is," and we cannot guarantee to 
C correct the errors reported by readers on any fixed schedule.

C      All the programs and subroutines  contained in this diskette 
C are copyrighted.   However,  we give permission to the reader who
C purchases this book to incorporate any of these programs into his
C or her programs provided that the copyright is acknowledged. 

C
C       ==================================================
C       Purpose: This program computes the psi function
C                using subroutine PSI
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       Examples:
C                   x          Psi(x)
C                ------------------------
C                  .25      -4.227453533
C                  .50      -1.963510026
C                  .75      -1.085860880
C                 1.00       -.577215665
C                 1.25       -.227453533
C                 1.50        .036489974
C                 1.75        .247472454
C                 2.00        .422784335
C       ==================================================
C
        SUBROUTINE PSI(X,PS)
C
C       ======================================
C       Purpose: Compute the psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
C http://jin.ece.uiuc.edu/routines/routines.html
C 		****************************************
C 		*           DISK TO ACCOMPANY          *
C 		*   COMPUTATION OF SPECIAL FUNCTIONS   *
C 		*                                      *
C 		*   Shanjie Zhang and Jianming Jin     *
C 		*                                      *
C 		*   Copyright 1996 by John Wiley &     *
C 		*              Sons, Inc.              *
C 		*                                      *
C 		****************************************
C
C I. INTRODUCTION
C
C      As stated in the preface of our book "Computation of Special 
C Functions,"  the purpose of this book is to share with the reader  
C a set of computer programs (130 in total) which we have developed 
C during the past several years for computing a variety of  special  
C mathematical functions.  For your convenience, we attach to the
C book this diskette that contains all the computer programs  
C listed or mentioned in the book. 

C      In this diskette,  we place all the programs under directory 
C SMF\PROGRAMS. In order to illustrate the use of these programs 
C and facilitate your testing of the programs, we wrote a short 
C simple main program for each program so that you can readily test 
C them.    

C      All the programs are written in FORTRAN-77 and tested on PCs
C and workstations. Therefore, they should run on any computer with 
C implementation of the FORTRAN-77 standard.  

C      Although we have made a great effort to test these programs,  
C we would not be surprised  to find some errors in them.  We would 
C appreciate it if you can bring to our attention any errors you find.
C You can do this by either writing us directly at the location
C (e-mail: j-jin1@uiuc.edu) or writing to the publisher, whose address 
C appears on the back cover of the book.  However, we must note that
C all these programs are sold "as is," and we cannot guarantee to 
C correct the errors reported by readers on any fixed schedule.

C      All the programs and subroutines  contained in this diskette 
C are copyrighted.   However,  we give permission to the reader who
C purchases this book to incorporate any of these programs into his
C or her programs provided that the copyright is acknowledged. 

C
C       ===================================================
C       Purpose: This program computes the gamma function
C                â(x) for x > 0 using subroutine LGAMA
C       Examples:
C                  x           â(x)
C                -------------------------
C                 0.5     .1772453851D+01
C                 2.5     .1329340388D+01
C                 5.0     .2400000000D+02
C                 7.5     .1871254306D+04
C                10.0     .3628800000D+06
C       ===================================================
C

        SUBROUTINE LGAMA(KF,X,GL)
C
C       ==================================================
C       Purpose: Compute gamma function â(x) or ln[â(x)]
C       Input:   x  --- Argument of â(x) ( x > 0 )
C                KF --- Function code
C                       KF=1 for â(x); KF=0 for ln[â(x)]
C       Output:  GL --- â(x) or ln[â(x)]
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(10)
        DATA A/8.333333333333333D-02,-2.777777777777778D-03,
     &         7.936507936507937D-04,-5.952380952380952D-04,
     &         8.417508417508418D-04,-1.917526917526918D-03,
     &         6.410256410256410D-03,-2.955065359477124D-02,
     &         1.796443723688307D-01,-1.39243221690590D+00/
        X0=X
        IF (X.EQ.1.0.OR.X.EQ.2.0) THEN
           GL=0.0D0
           GO TO 20
        ELSE IF (X.LE.7.0) THEN
           N=INT(7-X)
           X0=X+N
        ENDIF
        X2=1.0D0/(X0*X0)
        XP=6.283185307179586477D0
        GL0=A(10)
        DO 10 K=9,1,-1
10         GL0=GL0*X2+A(K)
        GL=GL0/X0+0.5D0*DLOG(XP)+(X0-.5D0)*DLOG(X0)-X0
        IF (X.LE.7.0) THEN
           DO 15 K=1,N
              GL=GL-DLOG(X0-1.0D0)
15            X0=X0-1.0D0
        ENDIF
20      IF (KF.EQ.1) GL=DEXP(GL)
        RETURN
        END

C  This is William Cody's erf implementations

      SUBROUTINE CALERF(ARG,RESULT,JINT)
C------------------------------------------------------------------
C
C This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
C   for a real argument  x.  It contains three FUNCTION type
C   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
C   and one SUBROUTINE type subprogram, CALERF.  The calling
C   statements for the primary entries are:
C
C                   Y=ERF(X)     (or   Y=DERF(X)),
C
C                   Y=ERFC(X)    (or   Y=DERFC(X)),
C   and
C                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
C
C   The routine  CALERF  is intended for internal packet use only,
C   all computations within the packet being concentrated in this
C   routine.  The function subprograms invoke  CALERF  with the
C   statement
C
C          CALL CALERF(ARG,RESULT,JINT)
C
C   where the parameter usage is as follows
C
C      Function                     Parameters for CALERF
C       call              ARG                  Result          JINT
C
C     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
C     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
C     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
C
C   The main computation evaluates near-minimax approximations
C   from "Rational Chebyshev approximations for the error function"
C   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
C   transportable program uses rational functions that theoretically
C   approximate  erf(x)  and  erfc(x)  to at least 18 significant
C   decimal digits.  The accuracy achieved depends on the arithmetic
C   system, the compiler, the intrinsic functions, and proper
C   selection of the machine-dependent constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   XMIN   = the smallest positive floating-point number.
C   XINF   = the largest positive finite floating-point number.
C   XNEG   = the largest negative argument acceptable to ERFCX;
C            the negative of the solution to the equation
C            2*exp(x*x) = XINF.
C   XSMALL = argument below which erf(x) may be represented by
C            2*x/sqrt(pi)  and above which  x*x  will not underflow.
C            A conservative value is the largest machine number X
C            such that   1.0 + X = 1.0   to machine precision.
C   XBIG   = largest argument acceptable to ERFC;  solution to
C            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
C            W(x) = exp(-x*x)/[x*sqrt(pi)].
C   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
C            machine precision.  A conservative value is
C            1/[2*sqrt(XSMALL)]
C   XMAX   = largest acceptable argument to ERFCX; the minimum
C            of XINF and 1/[sqrt(pi)*XMIN].
C
C   Approximate values for some important machines are:
C
C                          XMIN       XINF        XNEG     XSMALL
C
C  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
C  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
C  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
C  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
C  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
C  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
C
C
C                          XBIG       XHUGE       XMAX
C
C  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
C  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
C  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
C  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
C  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
C  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns  ERFC = 0      for  ARG .GE. XBIG;
C
C                       ERFCX = XINF  for  ARG .LT. XNEG;
C      and
C                       ERFCX = 0     for  ARG .GE. XMAX.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, EXP
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: March 19, 1990
C
C------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL
      DOUBLE PRECISION
     1     A,ARG,B,C,D,DEL,FOUR,HALF,P,ONE,Q,RESULT,SIXTEN,SQRPI,
     2     TWO,THRESH,X,XBIG,XDEN,XHUGE,XINF,XMAX,XNEG,XNUM,XSMALL,
     3     Y,YSQ,ZERO
      DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)
C------------------------------------------------------------------
C  Mathematical constants
C------------------------------------------------------------------
CS    DATA FOUR,ONE,HALF,TWO,ZERO/4.0E0,1.0E0,0.5E0,2.0E0,0.0E0/,
CS   1     SQRPI/5.6418958354775628695E-1/,THRESH/0.46875E0/,
CS   2     SIXTEN/16.0E0/
      DATA FOUR,ONE,HALF,TWO,ZERO/4.0D0,1.0D0,0.5D0,2.0D0,0.0D0/,
     1     SQRPI/5.6418958354775628695D-1/,THRESH/0.46875D0/,
     2     SIXTEN/16.0D0/
C------------------------------------------------------------------
C  Machine-dependent constants
C------------------------------------------------------------------
CS    DATA XINF,XNEG,XSMALL/3.40E+38,-9.382E0,5.96E-8/,
CS   1     XBIG,XHUGE,XMAX/9.194E0,2.90E3,4.79E37/
      DATA XINF,XNEG,XSMALL/1.79D308,-26.628D0,1.11D-16/,
     1     XBIG,XHUGE,XMAX/26.543D0,6.71D7,2.53D307/
C------------------------------------------------------------------
C  Coefficients for approximation to  erf  in first interval
C------------------------------------------------------------------
CS    DATA A/3.16112374387056560E00,1.13864154151050156E02,
CS   1       3.77485237685302021E02,3.20937758913846947E03,
CS   2       1.85777706184603153E-1/
CS    DATA B/2.36012909523441209E01,2.44024637934444173E02,
CS   1       1.28261652607737228E03,2.84423683343917062E03/
      DATA A/3.16112374387056560D00,1.13864154151050156D02,
     1       3.77485237685302021D02,3.20937758913846947D03,
     2       1.85777706184603153D-1/
      DATA B/2.36012909523441209D01,2.44024637934444173D02,
     1       1.28261652607737228D03,2.84423683343917062D03/
C------------------------------------------------------------------
C  Coefficients for approximation to  erfc  in second interval
C------------------------------------------------------------------
CS    DATA C/5.64188496988670089E-1,8.88314979438837594E0,
CS   1       6.61191906371416295E01,2.98635138197400131E02,
CS   2       8.81952221241769090E02,1.71204761263407058E03,
CS   3       2.05107837782607147E03,1.23033935479799725E03,
CS   4       2.15311535474403846E-8/
CS    DATA D/1.57449261107098347E01,1.17693950891312499E02,
CS   1       5.37181101862009858E02,1.62138957456669019E03,
CS   2       3.29079923573345963E03,4.36261909014324716E03,
CS   3       3.43936767414372164E03,1.23033935480374942E03/
      DATA C/5.64188496988670089D-1,8.88314979438837594D0,
     1       6.61191906371416295D01,2.98635138197400131D02,
     2       8.81952221241769090D02,1.71204761263407058D03,
     3       2.05107837782607147D03,1.23033935479799725D03,
     4       2.15311535474403846D-8/
      DATA D/1.57449261107098347D01,1.17693950891312499D02,
     1       5.37181101862009858D02,1.62138957456669019D03,
     2       3.29079923573345963D03,4.36261909014324716D03,
     3       3.43936767414372164D03,1.23033935480374942D03/
C------------------------------------------------------------------
C  Coefficients for approximation to  erfc  in third interval
C------------------------------------------------------------------
CS    DATA P/3.05326634961232344E-1,3.60344899949804439E-1,
CS   1       1.25781726111229246E-1,1.60837851487422766E-2,
CS   2       6.58749161529837803E-4,1.63153871373020978E-2/
CS    DATA Q/2.56852019228982242E00,1.87295284992346047E00,
CS   1       5.27905102951428412E-1,6.05183413124413191E-2,
CS   2       2.33520497626869185E-3/
      DATA P/3.05326634961232344D-1,3.60344899949804439D-1,
     1       1.25781726111229246D-1,1.60837851487422766D-2,
     2       6.58749161529837803D-4,1.63153871373020978D-2/
      DATA Q/2.56852019228982242D00,1.87295284992346047D00,
     1       5.27905102951428412D-1,6.05183413124413191D-2,
     2       2.33520497626869185D-3/
C------------------------------------------------------------------
      X = ARG
      Y = ABS(X)
      IF (Y .LE. THRESH) THEN
C------------------------------------------------------------------
C  Evaluate  erf  for  |X| <= 0.46875
C------------------------------------------------------------------
            YSQ = ZERO
            IF (Y .GT. XSMALL) YSQ = Y * Y
            XNUM = A(5)*YSQ
            XDEN = YSQ
            DO 20 I = 1, 3
               XNUM = (XNUM + A(I)) * YSQ
               XDEN = (XDEN + B(I)) * YSQ
   20       CONTINUE
            RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
            IF (JINT .NE. 0) RESULT = ONE - RESULT
            IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
            GO TO 800
C------------------------------------------------------------------
C  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
C------------------------------------------------------------------
         ELSE IF (Y .LE. FOUR) THEN
            XNUM = C(9)*Y
            XDEN = Y
            DO 120 I = 1, 7
               XNUM = (XNUM + C(I)) * Y
               XDEN = (XDEN + D(I)) * Y
  120       CONTINUE
            RESULT = (XNUM + C(8)) / (XDEN + D(8))
            IF (JINT .NE. 2) THEN
               YSQ = AINT(Y*SIXTEN)/SIXTEN
               DEL = (Y-YSQ)*(Y+YSQ)
               RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
            END IF
C------------------------------------------------------------------
C  Evaluate  erfc  for |X| > 4.0
C------------------------------------------------------------------
         ELSE
            RESULT = ZERO
            IF (Y .GE. XBIG) THEN
               IF ((JINT .NE. 2) .OR. (Y .GE. XMAX)) GO TO 300
               IF (Y .GE. XHUGE) THEN
                  RESULT = SQRPI / Y
                  GO TO 300
               END IF
            END IF
            YSQ = ONE / (Y * Y)
            XNUM = P(6)*YSQ
            XDEN = YSQ
            DO 240 I = 1, 4
               XNUM = (XNUM + P(I)) * YSQ
               XDEN = (XDEN + Q(I)) * YSQ
  240       CONTINUE
            RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
            RESULT = (SQRPI -  RESULT) / Y
            IF (JINT .NE. 2) THEN
               YSQ = AINT(Y*SIXTEN)/SIXTEN
               DEL = (Y-YSQ)*(Y+YSQ)
               RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
            END IF
      END IF
C------------------------------------------------------------------
C  Fix up for negative argument, erf, etc.
C------------------------------------------------------------------
  300 IF (JINT .EQ. 0) THEN
            RESULT = (HALF - RESULT) + HALF
            IF (X .LT. ZERO) RESULT = -RESULT
         ELSE IF (JINT .EQ. 1) THEN
            IF (X .LT. ZERO) RESULT = TWO - RESULT
         ELSE
            IF (X .LT. ZERO) THEN
               IF (X .LT. XNEG) THEN
                     RESULT = XINF
                  ELSE
                     YSQ = AINT(X*SIXTEN)/SIXTEN
                     DEL = (X-YSQ)*(X+YSQ)
                     Y = EXP(YSQ*YSQ) * EXP(DEL)
                     RESULT = (Y+Y) - RESULT
               END IF
            END IF
      END IF
  800 RETURN
C---------- Last card of CALERF ----------
      END
CS    REAL FUNCTION ERF(X)
      DOUBLE PRECISION FUNCTION DERF(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for erf(x).
C   (see comments heading CALERF).
C
C   Author/date: W. J. Cody, January 8, 1985
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
C------------------------------------------------------------------
      JINT = 0
      CALL CALERF(X,RESULT,JINT)
CS    ERF = RESULT
      DERF = RESULT
      RETURN
C---------- Last card of DERF ----------
      END
CS    REAL FUNCTION ERFC(X)
      DOUBLE PRECISION FUNCTION DERFC(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for erfc(x).
C   (see comments heading CALERF).
C
C   Author/date: W. J. Cody, January 8, 1985
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
C------------------------------------------------------------------
      JINT = 1
      CALL CALERF(X,RESULT,JINT)
CS    ERFC = RESULT
      DERFC = RESULT
      RETURN
C---------- Last card of DERFC ----------
      END
CS    REAL FUNCTION ERFCX(X)
      DOUBLE PRECISION FUNCTION DERFCX(X)
C------------------------------------------------------------------
C
C This subprogram computes approximate values for exp(x*x) * erfc(x).
C   (see comments heading CALERF).
C
C   Author/date: W. J. Cody, March 30, 1987
C
C------------------------------------------------------------------
      INTEGER JINT
CS    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
C------------------------------------------------------------------
      JINT = 2
      CALL CALERF(X,RESULT,JINT)
CS    ERFCX = RESULT
      DERFCX = RESULT
      RETURN
C---------- Last card of DERFCX ----------
      END
C     ALGORITHM 488 COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN COMM. ACM, VOL. 17, NO. 12,
C     P. 704.
      FUNCTION DGRAND(N)                                                 GRA   10
C EXCEPT ON THE FIRST CALL GRAND RETURNS A
C PSEUDO-RANDOM NUMBER HAVING A GAUSSIAN (I.E.
C NORMAL) DISTRIBUTION WITH ZERO MEAN AND UNIT
C STANDARD DEVIATION.  THUS, THE DENSITY IS  F(X) =
C EXP(-0.5*X**2)/SQRT(2.0*PI). THE FIRST CALL
C INITIALIZES GRAND AND RETURNS ZERO.
C THE PARAMETER N IS DUMMY.
C GRAND CALLS A FUNCTION RAND, AND IT IS ASSUMED THAT
C SUCCESSIVE CALLS TO RAND(0) GIVE INDEPENDENT
C PSEUDO- RANDOM NUMBERS DISTRIBUTED UNIFORMLY ON (0,
C 1), POSSIBLY INCLUDING 0 (BUT NOT 1).
C THE METHOD USED WAS SUGGESTED BY VON NEUMANN, AND
C IMPROVED BY FORSYTHE, AHRENS, DIETER AND BRENT.
C ON THE AVERAGE THERE ARE 1.37746 CALLS OF RAND FOR
C EACH CALL OF GRAND.
C WARNING - DIMENSION AND DATA STATEMENTS BELOW ARE
C           MACHINE-DEPENDENT.
C DIMENSION OF D MUST BE AT LEAST THE NUMBER OF BITS
C IN THE FRACTION OF A FLOATING-POINT NUMBER.
C THUS, ON MOST MACHINES THE DATA STATEMENT BELOW
C CAN BE TRUNCATED.
C IF THE INTEGRAL OF SQRT(2.0/PI)*EXP(-0.5*X**2) FROM
C A(I) TO INFINITY IS 2**(-I), THEN D(I) = A(I) -
C A(I-1).
      DOUBLE PRECISION D
      DIMENSION D(60)
      DATA D(1), D(2), D(3), D(4), D(5), D(6), D(7),
     * D(8), D(9), D(10), D(11), D(12), D(13),
     * D(14), D(15), D(16), D(17), D(18), D(19),
     * D(20), D(21), D(22), D(23), D(24), D(25),
     * D(26), D(27), D(28), D(29), D(30), D(31),
     * D(32) /0.674489750,0.475859630,0.383771164,
     * 0.328611323,0.291142827,0.263684322,
     * 0.242508452,0.225567444,0.211634166,
     * 0.199924267,0.189910758,0.181225181,
     * 0.173601400,0.166841909,0.160796729,
     * 0.155349717,0.150409384,0.145902577,
     * 0.141770033,0.137963174,0.134441762,
     * 0.131172150,0.128125965,0.125279090,
     * 0.122610883,0.120103560,0.117741707,
     * 0.115511892,0.113402349,0.111402720,
     * 0.109503852,0.107697617/
      DATA D(33), D(34), D(35), D(36), D(37), D(38),
     * D(39), D(40), D(41), D(42), D(43), D(44),
     * D(45), D(46), D(47), D(48), D(49), D(50),
     * D(51), D(52), D(53), D(54), D(55), D(56),
     * D(57), D(58), D(59), D(60)
     * /0.105976772,0.104334841,0.102766012,
     * 0.101265052,0.099827234,0.098448282,
     * 0.097124309,0.095851778,0.094627461,
     * 0.093448407,0.092311909,0.091215482,
     * 0.090156838,0.089133867,0.088144619,
     * 0.087187293,0.086260215,0.085361834,
     * 0.084490706,0.083645487,0.082824924,
     * 0.082027847,0.081253162,0.080499844,
     * 0.079766932,0.079053527,0.078358781,
     * 0.077681899/
C END OF MACHINE-DEPENDENT STATEMENTS
C U MUST BE PRESERVED BETWEEN CALLS.
      DATA U /0.0/
      DOUBLE PRECISION A
C INITIALIZE DISPLACEMENT A AND COUNTER I.
      A = 0.0
      I = 0
C INCREMENT COUNTER AND DISPLACEMENT IF LEADING BIT
C OF U IS ONE.
   10 U = U + U
      IF (U.LT.1.0) GO TO 20
      U = U - 1.0
      I = I + 1
      A = A - D(I)
      GO TO 10
C FORM W UNIFORM ON 0 .LE. W .LT. D(I+1) FROM U.
   20 W = D(I+1)*U
C FORM V = 0.5*((W-A)**2 - A**2). NOTE THAT 0 .LE. V
C .LT. LOG(2).
      V = W*(0.5*W-A)
C GENERATE NEW UNIFORM U.
   30 U = RAND(0)
C ACCEPT W AS A RANDOM SAMPLE IF V .LE. U.
      IF (V.LE.U) GO TO 40
C GENERATE RANDOM V.
      V = RAND(0)
C LOOP IF U .GT. V.
      IF (U.GT.V) GO TO 30
C REJECT W AND FORM A NEW UNIFORM U FROM V AND U.
      U = (V-U)/(1.0-U)
      GO TO 20
C FORM NEW U (TO BE USED ON NEXT CALL) FROM U AND V.
   40 U = (U-V)/(1.0-V)
C USE FIRST BIT OF U FOR SIGN, RETURN NORMAL VARIATE.
      U = U + U
      IF (U.LT.1.0) GO TO 50
      U = U - 1.0
      DGRAND = W - A
      RETURN
   50 DGRAND = A - W
      RETURN
      END

C Algorithm 467 from ACM (can't get this working at the moment).
C See http://pdp-10.trailing-edge.com/red405a2/11/uetp/lib/467.for
	SUBROUTINE DXPOSE(A, N1, N2, N12, MOVED, NWORK)
C  TRANSPOSITION OF A RECTANGULAR MATRIX IN SITU.
C  BY NORMAN BRENNER, MIT, 1/72.  CF. ALG. 380, CACM, 5/70.
C  TRANSPOSITION OF THE N1 BY N2 MATRIX A AMOUNTS TO
C  REPLACING THE ELEMENT AT VECTOR POSITION I (0-ORIGIN)
C  WITH THE ELEMENT AT POSITION N1*I (MOD N1*N2-I).
C  EACH SUBCYCLE OF THIS PERMUTATION IS COMPLETED IN ORDER.
C  MOVED IS A LOGICAL WORK ARRAY OF LENGTH NWORK.
	LOGICAL MOVED
	DOUBLE PRECISION A
	DOUBLE PRECISION ATEMP
	DOUBLE PRECISION BTEMP
	DIMENSION A(N12), MOVED(NWORK)
C  REALLY A(N1,N2), BUT N12 = N1*N2
	DIMENSION IFACT(8), IPOWER(8), NEXP(8), IEXP(8)
	IF (N1.LT.2 .OR. N2.LT.2) RETURN
	N = N1
	M = N1*N2 - 1
	IF (N1.NE.N2) GO TO 30
C  SQUARE MATRICES ARE DONE SEPARATELY FOR SPEED
	I1MIN = 2
	DO 20 I1MAX=N,M,N
	I2 = I1MIN + N - 1
	DO 10 I1=I1MIN,I1MAX
	ATEMP = A(I1)
	A(I1) = A(I2)
	A(I2) = ATEMP
	I2 = I2 + N
10	CONTINUE
	I1MIN = I1MIN + N + 1
20	CONTINUE
	RETURN
C  MODULUS M IS FACTORED INTO PRIME POWERS.  EIGHT FACTORS
C  SUFFICE UP TO M = 2*3*5*7*11*13*17*19 = 9,767,520.
30	CALL FACTOR(M, IFACT, IPOWER, NEXP, NPOWER)
	DO 40 IP=1,NPOWER
	IEXP(IP) = 0
40	CONTINUE
C  GENERATE EVERY DIVISOR OF M LESS THAN M/2
	IDIV = 1
50	IF (IDIV.GE.M/2) GO TO 190
C  THE NUMBER OF ELEMENTS WHOSE INDEX IS DIVISIBLE BY IDIV
C  AND BY NO OTHER DIVISOR OF M IS THE EULER TOTIENT
C  FUNCTION, PHI(M/IDIV).
	NCOUNT = M/IDIV
	DO 60 IP=1,NPOWER
	IF (IEXP(IP).EQ.NEXP(IP)) GO TO 60
	NCOUNT = (NCOUNT/IFACT(IP))*(IFACT(IP)-1)
60	CONTINUE
	DO 70 I=1,NWORK
	MOVED(I) = .FALSE.
70	CONTINUE
C  THE STARTING POINT OF A SUBCYCLE IS DIVISIBLE ONLY BY IDIV
C  AND MUST NOT APPEAR IN ANY OTHER SUBCYCLE.
	ISTART = IDIV
80	MMIST = M - ISTART
	IF (ISTART.EQ.IDIV) GO TO 120
	IF (ISTART.GT.NWORK) GO TO 90
	IF (MOVED(ISTART)) GO TO 160
90	ISOID = ISTART/IDIV
	DO 100 IP=1,NPOWER
	IF (IEXP(IP).EQ.NEXP(IP)) GO TO 100
	IF (MOD(ISOID,IFACT(IP)).EQ.0) GO TO 160
100	CONTINUE
	IF (ISTART.LE.NWORK) GO TO 120
	ITEST = ISTART
110	ITEST = MOD(N*ITEST,M)
	IF (ITEST.LT.ISTART .OR. ITEST.GT.MMIST) GO TO 160
	IF (ITEST.GT.ISTART .AND. ITEST.LT.MMIST) GO TO 110
120	ATEMP = A(ISTART+1)
	BTEMP = A(MMIST+I)
	IA1 = ISTART
130	IA2 = MOD(N*IA1,M)
	MMIA1 = M - IA1
	MMIA2 = M - IA2
	IF (IA1.LE.NWORK) MOVED(IA1) = .TRUE.
	IF (MMIA1.LE.NWORK) MOVED(MMIA1) = .TRUE.
	NCOUNT = NCOUNT - 2
C  MOVE TWO ELEMENTS, THE SECOND FROM THE NEGATIVE
C  SUBCYCLE.  CHECK FIRST FOR SUBCYCLE CLOSURE.
	IF (IA2.EQ.ISTART) GO TO 140
	IF (MMIA2.EQ.ISTART) GO TO 150
	A(IA1+1) = A(IA2+1)
	A(MMIA1+1) = A(MMIA2+1)
	IA1 = IA2
	GO TO 130
140	A(IA1+1) = ATEMP
	A(MMIA1+1) = BTEMP
	GO TO 160
150	A(IA1+1) = BTEMP
	A(MMIA1+1) = ATEMP
160	ISTART = ISTART + IDIV
	IF (NCOUNT.GT.0) GO TO 80
	DO 180 IP=1,NPOWER
	IF (IEXP(IP).EQ.NEXP(IP)) GO TO 170
	IEXP(IP) = IEXP(IP) + 1
	IDIV = IDIV*IFACT(IP)
	GO TO 50
170	IEXP(IP) = 0
	IDIV = IDIV/IPOWER(IP)
180	CONTINUE
190	RETURN
	END

	SUBROUTINE FACTOR(N, IFACT, IPOWER, NEXP, NPOWER)
C  FACTOR N INTO ITS PRIME POWERS, NPOWER IN NUMBER.
C  E.G., FOR N=1970=2**3 *5 *7**2, NPOWER=3, IFACT=3,5,7,
C  IPOWER=8,5,49, AND NEXP=3,1,2.
	DIMENSION IFACT(8), IPOWER(8), NEXP(8)
	IP = 0
	IFCUR = 0
	NPART = N
	IDIV = 2
10	IQUOT = NPART/IDIV
	IF (NPART-IDIV*IQUOT) 60, 20, 60
20	IF (IDIV-IFCUR) 40, 40, 30
30	IP = IP + 1
	IFACT(IP) = IDIV
	IPOWER(IP) = IDIV
	IFCUR = IDIV
	NEXP(IP) = 1
	GO TO 50
40	IPOWER(IP) = IDIV*IPOWER(IP)
	NEXP(IP) = NEXP(IP) + 1
50	NPART = IQUOT
	GO TO 10
60	IF (IQUOT-IDIV) 100, 100, 70
70	IF (IDIV-2) 80, 80, 90
80	IDIV = 3
	GO TO 10
 90	IDIV = IDIV + 2
	GO TO 10
100	IF (NPART-1) 140, 140, 110
110	IF (NPART-IFCUR) 130, 130, 120
120	IP = IP + 1
	IFACT(IP) = NPART
	IPOWER(IP) = NPART
	NEXP(IP) = 1
	GO TO 140
130	IPOWER(IP) = NPART*IPOWER(IP)
	NEXP(IP) = NEXP(IP) + 1
140	NPOWER = IP
	RETURN
	END

C ACM Algorithm 380
      SUBROUTINE DTRANS(A,M,N,MN,MOVE,IWRK,IOK)
C A IS A ONE-DIMENSIONAL ARRAY OF LENGTH MN=M*N, WHICH
C CONTAINS THE M BY N MATRIX TO BE TRANSPOSED (STORED
C COLUMNWISE).MOVE IS A ONE-DIMENSIONAL ARRAY OF LENGTH IWRK
C USED TO STORE INFORMATION TO SPEED UP THE PROCESS. THE
C VALUE IWRK=(M+N)/2 IS RECOMMENDED. IOK INDICATES THE
C SUCCESS OR FAILURE OF THE ROUTINE.
C NORMAL RETURN IOK=0
C ERRORS           IOK=-1, MN NOT EQUAL TO M*N.
C                  IOK=-2, IWRK NEGATIVE OR ZERO.
C                  IOK.GT.0, (SHOULD NEVER OCCUR). IN THIS CASE
C WE SET IOK EQUAL TO THE FINAL VALUE OF I WHEN THE SEARCH
C IS COMPLETED BUT SOME LOOPS HAVE NOT BEEN MOVED.
      DOUBLE PRECISION A
      DIMENSION A(MN),MOVE(IWRK)
C CHECK ARGUMENTS AND INITIALISE
      IF(M.LT.2.OR.N.LT.2) GO TO 60
      IF(MN.NE.M*N) GO TO 92
      IF(IWRK.LT.1) GO TO 93
      IF(M.EQ.N) GO TO 70
      NCOUNT=2
      M2=M-2
      DO 10 I=1,IWRK
   10   MOVE(I)=0
      IF(M2.LT.1) GO TO 12
C COUNT NUMBER,NCOUNT, OF SINGLE POINTS.
      DO 11 IA=1,M2
        IB=IA*(N-1)/(M-1)
        IF(IA*(N-1).NE.IB*(M-1)) GO TO 11
        NCOUNT=NCOUNT+1
        I=IA*N+IB
        IF(I.GT.IWRK) GO TO 11
        MOVE(I)=1
   11   CONTINUE
C SET INITIAL VALUES FOR SEARCH.
   12 K=MN-1
      KMI=K-1
      MAX=MN
      I=1
C AT LEAST ONE LOOP MUST BE RE-ARRANGED.
      GO TO 30
C SEARCH FOR LOOPS TO BE REARRANGED.
   20        MAX=K-I
             I=I+1
             KMI=K-I
             IF(I.GT.MAX) GO TO 90
             IF(I.GT.IWRK) GO TO 21
             IF(MOVE(I).LT.1) GO TO 30
             GO TO 20
   21        IF(I.EQ.M*I-K*(I/N)) GO TO 20
             I1=I
   22        I2=M*I1-K*(I1/N)
             IF(I2.LE.I .OR. I2.GE.MAX)  GO TO 23
             I1=I2
             GO TO 22
   23        IF(I2.NE.I) GO TO 20
C REARRANGE ELEMENTS OF A LOOP.
   30             I1=I
   31             B=A(I1+1)
   32             I2=M*I1-K*(I1/N)
                  IF(I1.LE.IWRK) MOVE(I1)=2
   33             NCOUNT=NCOUNT+1
                  IF(I2.EQ.I .OR. I2.GE.KMI) GO TO 35
   34             A(I1+1)=A(I2+1)
                  I1=I2
                  GO TO 32
   35             IF(MAX.EQ.KMI .OR. I2.EQ.I) GO TO 41
                  MAX=KMI
                  GO TO 34
C TEST FOR SYMMETRIC PAIR OF LOOPS.
   41        A(I1+1)=B
             IF(NCOUNT.GE.MN) GO TO 60
             IF(I2.EQ.MAX .OR. MAX.EQ.KMI) GO TO 20
             MAX=KMI
             I1=MAX
             GO TO 31
C NORMAL RETURN.
   60 IOK=0
      RETURN
C IF MATRIX IS SQUARE, EXCHANGE ELEMENTS A(I,J) AND A(J,I).
   70 N1=N-1
      DO 71 I=1,N1
        J1=I+1
        DO 71 J=J1,N
          I1=I+(J-1)*N
          I2=J+(I-1)*M
          B=A(I1)
          A(I1)=A(I2)
          A(I2)=B
   71     CONTINUE
      GO TO 60
C ERROR RETURNS.
   90 IOK=I
   91 RETURN
   92 IOK=-1
      GO TO 91
   93 IOK=-2
      GO TO 91
      END

C ACM Algorithm 380.
      SUBROUTINE DTRANSR(A, M, N, MN, MOVE, IWRK, IOK)                    TRA   10
C *****
C  ALGORITHM 380 - REVISED
C *****
C  A IS A ONE-DIMENSIONAL ARRAY OF LENGTH MN=M*N, WHICH
C  CONTAINS THE MXN MATRIX TO BE TRANSPOSED (STORED
C  COLUMWISE). MOVE IS A ONE-DIMENSIONAL ARRAY OF LENGTH IWRK
C  USED TO STORE INFORMATION TO SPEED UP THE PROCESS.  THE
C  VALUE IWRK=(M+N)/2 IS RECOMMENDED. IOK INDICATES THE
C  SUCCESS OR FAILURE OF THE ROUTINE.
C  NORMAL RETURN  IOK=0
C  ERRORS         IOK=-1 ,MN NOT EQUAL TO M*N
C                 IOK=-2 ,IWRK NEGATIVE OR ZERO
C                 IOK.GT.0, (SHOULD NEVER OCCUR),IN THIS CASE
C  WE SET IOK EQUAL TO THE FINAL VALUE OF I WHEN THE SEARCH
C  IS COMPLETED BUT SOME LOOPS HAVE NOT BEEN MOVED
C  NOTE * MOVE(I) WILL STAY ZERO FOR FIXED POINTS
      DOUBLE PRECISION A
      DIMENSION A(MN), MOVE(IWRK)
C CHECK ARGUMENTS AND INITIALIZE.
      IF (M.LT.2 .OR. N.LT.2) GO TO 120
      IF (MN.NE.M*N) GO TO 180
      IF (IWRK.LT.1) GO TO 190
      IF (M.EQ.N) GO TO 130
      NCOUNT = 2
      K = MN - 1
      DO 10 I=1,IWRK
        MOVE(I) = 0
   10 CONTINUE
      IF (M.LT.3 .OR. N.LT.3) GO TO 30
C CALCULATE THE NUMBER OF FIXED POINTS, EUCLIDS ALGORITHM
C FOR GCD(M-1,N-1).
      IR2 = M - 1
      IR1 = N - 1
   20 IR0 = MOD(IR2,IR1)
      IR2 = IR1
      IR1 = IR0
      IF (IR0.NE.0) GO TO 20
      NCOUNT = NCOUNT + IR2 - 1
C SET INITIAL VALUES FOR SEARCH
   30 I = 1
      IM = M
C AT LEAST ONE LOOP MUST BE RE-ARRANGED
      GO TO 80
C SEARCH FOR LOOPS TO REARRANGE
   40 MAX = K - I
      I = I + 1
      IF (I.GT.MAX) GO TO 160
      IM = IM + M
      IF (IM.GT.K) IM = IM - K
      I2 = IM
      IF (I.EQ.I2) GO TO 40
      IF (I.GT.IWRK) GO TO 60
      IF (MOVE(I).EQ.0) GO TO 80
      GO TO 40
   50 I2 = M*I1 - K*(I1/N)
   60 IF (I2.LE.I .OR. I2.GE.MAX) GO TO 70
      I1 = I2
      GO TO 50
   70 IF (I2.NE.I) GO TO 40
C REARRANGE THE ELEMENTS OF A LOOP AND ITS COMPANION LOOP
   80 I1 = I
      KMI = K - I
      B = A(I1+1)
      I1C = KMI
      C = A(I1C+1)
   90 I2 = M*I1 - K*(I1/N)
      I2C = K - I2
      IF (I1.LE.IWRK) MOVE(I1) = 2
      IF (I1C.LE.IWRK) MOVE(I1C) = 2
      NCOUNT = NCOUNT + 2
      IF (I2.EQ.I) GO TO 110
      IF (I2.EQ.KMI) GO TO 100
      A(I1+1) = A(I2+1)
      A(I1C+1) = A(I2C+1)
      I1 = I2
      I1C = I2C
      GO TO 90
C FINAL STORE AND TEST FOR FINISHED
  100 D = B
      B = C
      C = D
  110 A(I1+1) = B
      A(I1C+1) = C
      IF (NCOUNT.LT.MN) GO TO 40
C NORMAL RETURN
  120 IOK = 0
      RETURN
C IF MATRIX IS SQUARE,EXCHANGE ELEMENTS A(I,J) AND A(J,I).
  130 N1 = N - 1
      DO 150 I=1,N1
        J1 = I + 1
        DO 140 J=J1,N
          I1 = I + (J-1)*N
          I2 = J + (I-1)*M
          B = A(I1)
          A(I1) = A(I2)
          A(I2) = B
  140   CONTINUE
  150 CONTINUE
      GO TO 120
C ERROR RETURNS.
  160 IOK = I
  170 RETURN
  180 IOK = -1
      GO TO 170
  190 IOK = -2
      GO TO 170
      END
