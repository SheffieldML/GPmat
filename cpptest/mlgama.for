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
