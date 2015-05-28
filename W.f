	  SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
	  IMPLICIT NONE
	  INTEGER PLHS(*), PRHS(*)
	  INTEGER NLHS, NRHS
	  INTEGER MXCREATEFULL, MXGETPR
	  INTEGER MXGETM, MXGETN

	  integer M(1), N(1)

	  M(1) = MXGETM(PRHS(1))
	  N(1) = MXGETN(PRHS(1))
	
!	CALL INTERMEDIATE GATEWAY:
	  CALL GATEWAY(PLHS,PRHS,M(1))

	  END 
!	==============================================================

	  SUBROUTINE GATEWAY(PLHS, PRHS, T)
	  IMPLICIT NONE
      INTEGER PLHS(1), PRHS(1), T, I
	  INTEGER MXCREATEFULL, MXGETPR, MXGETPI
	  INTEGER MXGETM, MXGETN
	  DOUBLE PRECISION A(T), B(T), C(T), D(T)
      LOGICAL FLAG   

      call mxcopyptrtoreal8(mxgetpr(prhs(1)),A,T)
	  call mxcopyptrtoreal8(mxgetpi(prhs(1)),B,T)

	  plhs(1) = mxcreatefull(T,1,1)

!	  CALL COMPUTATIONAL ROUTINE:
      DO 20, I=1,T 
      CALL WOFZ(A(I),B(I),C(I),D(I),FLAG)
20    CONTINUE
 	  CALL MXCOPYreal8TOPTR(C,mxgetpr(plhs(1)),T)
      CALL MXCOPYreal8TOPTR(D,mxgetpi(plhs(1)),T)
	  END
!	===============================================================
      SUBROUTINE WOFZ (XI, YI, U, V, FLAG)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
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
