#include "fintrf.h"
C
#if 0
C     
C     
#endif	
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
	  IMPLICIT NONE
	  INTEGER PLHS(*), PRHS(*)
	  INTEGER MXGETM, MXGETN
      INTEGER NLHS, NRHS
	  INTEGER M, N


C      Check for proper number of arguments. 
      IF (NRHS .NE. 4) THEN
         CALL MEXERRMSGTXT('Four inputs required.')
      ENDIF

C	  Check that input #1 is a scalar.
      M = MXGETM(PRHS(1))
      N = MXGETN(PRHS(1))
      IF(M .NE. 1 .OR. N .NE. 1) THEN
         CALL MEXERRMSGTXT('Input #1 is not a scalar.')
      ENDIF

C	  Check that input #2 is a scalar.
      M = MXGETM(PRHS(2))
      N = MXGETN(PRHS(2))
      IF(M .NE. 1 .OR. N .NE. 1) THEN
         CALL MEXERRMSGTXT('Input #2 is not a scalar.')
      ENDIF

C	  Check that input #3 is a column vector.
      N = MXGETN(PRHS(3))
      IF( N .NE. 1) THEN
         CALL MEXERRMSGTXT('Input #3 is not a column vector.')
      ENDIF

C	  Check that input #4 is a column vector.
      N = MXGETN(PRHS(4))
      IF( N .NE. 1) THEN
         CALL MEXERRMSGTXT('Input #4 is not a column vector.')
      ENDIF

	  M = MXGETM(PRHS(3))
	  N = MXGETM(PRHS(4))
	
C	CALL INTERMEDIATE GATEWAY:
	  CALL GATEWAY(PLHS,PRHS,M,N)

	  END 
C ------ GATEWAY SUBROUTINE

	  SUBROUTINE GATEWAY(PLHS, PRHS, ROWS, COLS)
	  IMPLICIT NONE
      INTEGER PLHS(*), PRHS(*)
	  INTEGER MXCREATEDOUBLEMATRIX, MXGETPR, MXGETPI
	  INTEGER MXISCOMPLEX
      INTEGER ROWS, COLS, ISCOMPLEX, SIZE
	  REAL *8 SIG2, T1(ROWS), T2(COLS), ALPHA, OMEGA
	  REAL *8 RGRAD(ROWS, COLS), IGRAD(ROWS, COLS)
      COMPLEX*16 GAMMA

      SIZE = ROWS*COLS
      
      PLHS(1) = MXCREATEDOUBLEMATRIX(ROWS,COLS,1)
      CALL MXCOPYPTRTOREAL8(MXGETPR(PRHS(1)), ALPHA ,1)  

C ------  Convert MATLAB variables to FORTRAN variables  
      
      ISCOMPLEX = MXISCOMPLEX(PRHS(1))       
      IF(ISCOMPLEX.EQ.1) THEN
         CALL MXCOPYPTRTOREAL8(MXGETPI(PRHS(1)),OMEGA,1)          
      ELSE 
         OMEGA = 0
      END IF        
      CALL MXCOPYPTRTOREAL8(MXGETPR(PRHS(2)), SIG2 ,1)
      CALL MXCOPYPTRTOREAL8(MXGETPR(PRHS(3)), T1 ,ROWS)
      CALL MXCOPYPTRTOREAL8(MXGETPR(PRHS(4)), T2 ,COLS)
      GAMMA = DCMPLX(ALPHA,OMEGA)

C ------ Call computational subroutine     

      CALL COMPUP(GAMMA,SIG2, T1, T2, ROWS, COLS, RGRAD, IGRAD )   

C ------  Convert FORTRAN variables to MATLAB variables  

      CALL MXCOPYREAL8TOPTR(RGRAD,MXGETPR(PLHS(1)),SIZE)
      CALL MXCOPYREAL8TOPTR(IGRAD,MXGETPI(PLHS(1)),SIZE)

      RETURN 
      END

C ----- Computational subroutine
      SUBROUTINE COMPUP(GAMMA, SIG2, T1, T2, ROWS, COLS, 
     +   RGRAD, IGRAD)
	  INTEGER I,J, ROWS, COLS 
      REAL *8 RGRAD(ROWS, COLS), IGRAD(ROWS, COLS), PI, SQRTPI
	  REAL *8 SIG2, T1(ROWS), T2(COLS) 
      REAL *8 DEV, REALZ1, REALZ2(COLS) 
	  COMPLEX *16 GAMMA, WOFZ1, WOFZ2(COLS)
      COMPLEX *16 Z1, Z2(COLS), GRAD
	  
      PI = 3.141592653589793D0
      SQRTPI = DSQRT(PI)
      DEV = DSQRT(SIG2)
      DO 20, I=1,ROWS 
        DO 10, J=1,COLS             
             IF(I.EQ.1) THEN
                 Z2(J) = T2(J)/DEV + DEV*GAMMA/2
                 REALZ2(J) = DREAL(Z2(J))
                 IF(REALZ2(J).GE.0) THEN    
                    CALL WOFZ(DCMPLX(0,1)*Z2(J), WOFZ2(J))
                 ELSE
                    CALL WOFZ(DCMPLX(0,-1)*Z2(J), WOFZ2(J))
                 END IF   
             END IF  
             Z1 = (T1(I)-T2(J))/DEV - DEV*GAMMA/2   
             REALZ1 = DREAL(Z1)
             IF(REALZ1.GE.0) THEN    
                CALL WOFZ(DCMPLX(0,1)*Z1, WOFZ1)
             ELSE
                CALL WOFZ(DCMPLX(0,-1)*Z1, WOFZ1)
             END IF
C		real(z1) >= 0 and real(z2)>=0
             IF(REALZ1. GE. 0. AND .REALZ2(J). GE .0  ) THEN
                GRAD = DEV*(GAMMA**2)*CDEXP(SIG2*(GAMMA**2)/4
     +                  - GAMMA*(T1(I)-T2(J))) 
     +                  - 2*CDEXP(-((T1(I)-T2(J))**2/SIG2) 
     +                  + CDLOG(((T1(I)-T2(J))**2)*WOFZ1/(DEV**3) 
     +                  + ((T1(I)-T2(J))/SIG2+GAMMA/2)
     +                  * (1/SQRTPI	-Z1*WOFZ1))) 
     +                  - 2*CDEXP(-((T2(J)**2)/SIG2)-GAMMA*T1(I) 
     +                  + CDLOG((T2(J)**2)*WOFZ2(J)/(DEV**3) 
     +                  + (T2(J)/SIG2-GAMMA/2)*
     +                  ((1/SQRTPI)-Z2(J)*WOFZ2(J))))
              END IF
C		real(z1) < 0 and real(z2)>=0
		   IF(REALZ1. LT. 0. AND .REALZ2(J). GE .0  ) THEN
            GRAD =   2*CDEXP(-(T1(I)-T2(J))**2/SIG2 
     +                  + CDLOG(((T1(I)-T2(J))**2)*WOFZ1/(DEV**3) 
     +                  - ((T1(I)-T2(J))/SIG2+GAMMA/2)
     +                  * (1/SQRTPI	+Z1*WOFZ1))) 
     +                  - 2*CDEXP(-(T2(J)**2)/SIG2-GAMMA*T1(I) 
     +                  + CDLOG((T2(J)**2)*WOFZ2(J)/(DEV**3) 
     +                  + (T2(J)/SIG2-GAMMA/2)*
     +                  (1/SQRTPI-Z2(J)*WOFZ2(J))))
             END IF
C		real(z1) >= 0 and real(z2)<0
		   IF(REALZ1. GE. 0. AND .REALZ2(J). LT .0  ) THEN
              GRAD =  - 2*CDEXP(-(T1(I)-T2(J))**2/SIG2 
     +                  + CDLOG(((T1(I)-T2(J))**2)*WOFZ1/(DEV**3) 
     +                  + ((T1(I)-T2(J))/SIG2+GAMMA/2)
     +                  * (1/SQRTPI	-Z1*WOFZ1))) 
     +                  + 2*CDEXP(-(T2(J)**2)/SIG2-GAMMA*T1(I) 
     +                  + CDLOG((T2(J)**2)*WOFZ2(J)/(DEV**3) 
     +                  - (T2(J)/SIG2-GAMMA/2)* 
     +                  (1/SQRTPI+Z2(J)*WOFZ2(J))))
             END IF
C		real(z1) < 0 and real(z2)<0             
		   IF(REALZ1. LT. 0. AND .REALZ2(J). LT .0  ) THEN
            GRAD = - DEV*(GAMMA**2)*CDEXP(SIG2*(GAMMA**2)/4
     +                  - GAMMA*(T1(I)-T2(J))) 
     +                  + 2*CDEXP(-(T1(I)-T2(J))**2/SIG2 
     +                  + CDLOG(((T1(I)-T2(J))**2)*WOFZ1/(DEV**3) 
     +                  - ((T1(I)-T2(J))/SIG2+GAMMA/2)
     +                  * (1/SQRTPI	+Z1*WOFZ1))) 
     +                  + 2*CDEXP(-((T2(J)**2)/SIG2-GAMMA*T1(I)) 
     +                  + CDLOG((T2(J)**2)*WOFZ2(J)/(DEV**3) 
     +                  - (T2(J)/SIG2-GAMMA/2)*
     +                  ((1/SQRTPI)+Z2(J)*WOFZ2(J))))
           END IF
           RGRAD(I,J) = DREAL(GRAD) 
           IGRAD(I,J) = DIMAG(GRAD)
10      CONTINUE
20    CONTINUE
      RETURN  
      END

C-----  Complex error function due to Poppe 
      SUBROUTINE WOFZ (Z, W)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION XI, YI, U, V
      COMPLEX*16  Z, W
      LOGICAL A, B, FLAG
      PARAMETER (FACTOR   = 1.12837916709551257388D0,
     +          RMAXREAL = 0.5D+154,
     +          RMAXEXP  = 708.503061461606D0,
     +          RMAXGONI = 3.53711887601422D+15)
      XI = DREAL(Z)
      YI = DIMAG(Z) 
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
     +        (XQUAD.GT.RMAXEXP)) GOTO 100
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
      W = DCMPLX(U, V)
      RETURN
  100 FLAG = .TRUE.
      RETURN
      END
