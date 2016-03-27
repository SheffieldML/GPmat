#include "fintrf.h"
C
#if 0
C     
C     
#endif
	SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
	IMPLICIT NONE
	MWPOINTER PLHS(*), PRHS(*)
	INTEGER*4 MXISCOMPLEX
	INTEGER NLHS, NRHS
	MWPOINTER M, N
	MWSIZE SIZE
	MWSIZE MXGETM, MXGETN

C       Check for proper number of arguments. 
	IF (NRHS .NE. 1) THEN
	   CALL MEXERRMSGTXT('One input required.')
	ENDIF
C	Check if input #1 is complex
	IF (MXISCOMPLEX(PRHS(1)) .NE. 1) THEN
	   CALL MEXERRMSGTXT('Input array is not complex')
	ENDIF
c$$$C       Check that input #1 is a column vector.
c$$$	N = MXGETN(PRHS(1))
c$$$	IF( N .NE. 1) THEN
c$$$	   CALL MEXERRMSGTXT('Input should be a column vector.')
c$$$	ENDIF
	
	M = MXGETM(PRHS(1))
	N = MXGETN(PRHS(1))
	SIZE = M*N
!	CALL INTERMEDIATE GATEWAY:
	CALL GATEWAY(PLHS,PRHS,SIZE)
	
	END 
!	==============================================================

	SUBROUTINE GATEWAY(PLHS, PRHS, SIZE)
	IMPLICIT NONE
	MWPOINTER PLHS(1), PRHS(1)
	MWSIZE SIZE
	INTEGER I
	MWPOINTER ROWS, COLS
	MWPOINTER MXCREATEDOUBLEMATRIX, MXGETPR, MXGETPI
	INTEGER MXGETM, MXGETN
	REAL*8 A(SIZE), B(SIZE), C(SIZE), D(SIZE)
	LOGICAL FLAG   
	CHARACTER*120 LINE

	CALL MXCOPYPTRTOREAL8(MXGETPR(PRHS(1)),A,SIZE)
	CALL MXCOPYPTRTOREAL8(MXGETPI(PRHS(1)),B,SIZE)
	ROWS = MXGETM(PRHS(1))
	COLS = MXGETN(PRHS(1))
	PLHS(1) = MXCREATEDOUBLEMATRIX(ROWS,COLS,1)

!       CALL COMPUTATIONAL ROUTINE:
	DO 20, I=1,SIZE
	   CALL WOFZ(A(I),B(I),C(I),D(I),FLAG)
 20	CONTINUE

	CALL MXCOPYREAL8TOPTR(C,MXGETPR(PLHS(1)),SIZE)
	CALL MXCOPYREAL8TOPTR(D,MXGETPI(PLHS(1)),SIZE)
	END
!	===============================================================
	SUBROUTINE WOFZ (XI, YI, U, V, FLAG)
      IMPLICIT REAL*8 (A-H, O-Z)
      LOGICAL A, B, FLAG
      PARAMETER (FACTOR   = 1.12837916709551257388D0,
     $          RMAXREAL = 0.5D+154,
     $          RMAXEXP  = 708.503061461606D0,
     $          RMAXGONI = 3.53711887601422D+15)
      FLAG = .FALSE.
      XABS = DABS(XI)
      YABS = DABS(YI)
      X    = XABS/6.3
      Y    = YABS/4.4
      IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) GOTO 100
      QRHO = X**2 + Y**2
      XABSQ = XABS**2
      XQUAD = XABSQ - YABS**2
      YQUAD = 2*XABS*YABS
      A     = QRHO.LT.0.085264D0
      IF (A) THEN
        QRHO  = (1-0.85*Y)*DSQRT(QRHO)
        N     = IDNINT(6 + 72*QRHO)
        J     = 2*N+1
        XSUM  = 1.0/J
        YSUM  = 0.0D0
        DO 10 I=N, 1, -1
          J    = J - 2
          XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I
          YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I
          XSUM = XAUX + 1.0/J
 10     CONTINUE
        U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0
        V1   =  FACTOR*(XSUM*XABS - YSUM*YABS)
        DAUX =  DEXP(-XQUAD)
        U2   =  DAUX*DCOS(YQUAD)
        V2   = -DAUX*DSIN(YQUAD)
        U    = U1*U2 - V1*V2
        V    = U1*V2 + V1*U2
      ELSE
        IF (QRHO.GT.1.0) THEN
          H    = 0.0D0
          KAPN = 0
          QRHO = DSQRT(QRHO)
          NU   = IDINT(3 + (1442/(26*QRHO+77)))
        ELSE
          QRHO = (1-Y)*DSQRT(1-QRHO)
          H    = 1.88*QRHO
          H2   = 2*H
          KAPN = IDNINT(7  + 34*QRHO)
          NU   = IDNINT(16 + 26*QRHO)
        ENDIF
        B = (H.GT.0.0)
        IF (B) QLAMBDA = H2**KAPN
        RX = 0.0
        RY = 0.0
        SX = 0.0
        SY = 0.0
        DO 11 N=NU, 0, -1
          NP1 = N + 1
          TX  = YABS + H + NP1*RX
          TY  = XABS - NP1*RY
          C   = 0.5/(TX**2 + TY**2)
          RX  = C*TX
          RY  = C*TY
          IF ((B).AND.(N.LE.KAPN)) THEN
            TX = QLAMBDA + SX
            SX = RX*TX - RY*SY
            SY = RY*TX + RX*SY
            QLAMBDA = QLAMBDA/H2
          ENDIF
 11     CONTINUE
        IF (H.EQ.0.0) THEN
          U = FACTOR*RX
          V = FACTOR*RY
        ELSE
          U = FACTOR*SX
          V = FACTOR*SY
        END IF
        IF (YABS.EQ.0.0) U = DEXP(-XABS**2)
      END IF
      IF (YI.LT.0.0) THEN
        IF (A) THEN
          U2    = 2*U2
          V2    = 2*V2
        ELSE
          XQUAD =  -XQUAD
          IF ((YQUAD.GT.RMAXGONI).OR.
     *        (XQUAD.GT.RMAXEXP)) GOTO 100
          W1 =  2*DEXP(XQUAD)
          U2  =  W1*DCOS(YQUAD)
          V2  = -W1*DSIN(YQUAD)
        END IF
        U = U2 - U
        V = V2 - V
        IF (XI.GT.0.0) V = -V
      ELSE
        IF (XI.LT.0.0) V = -V
      END IF
      RETURN
  100 FLAG = .TRUE.
      RETURN
      END
