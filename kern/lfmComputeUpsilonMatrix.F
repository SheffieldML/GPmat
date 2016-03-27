#include "fintrf.h"
C
#if 0
C     
C     
#endif
	  SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
	  IMPLICIT NONE
	  MWPOINTER PLHS(*), PRHS(*)
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
	  MWPOINTER PLHS(*), PRHS(*)
	  MWPOINTER MXCREATEDOUBLEMATRIX, MXGETPR, MXGETPI
	  INTEGER MXISCOMPLEX
      INTEGER ROWS, COLS, ISCOMPLEX, SIZE
	  DOUBLE PRECISION SIG2, T1(ROWS), T2(COLS), ALPHA, OMEGA
	  DOUBLE PRECISION RUPSI(ROWS, COLS), IUPSI(ROWS, COLS)
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

      CALL COMPUP(GAMMA,SIG2, T1, T2, ROWS, COLS, RUPSI, IUPSI )   

C ------  Convert FORTRAN variables to MATLAB variables  

      CALL MXCOPYREAL8TOPTR(RUPSI,MXGETPR(PLHS(1)),SIZE)
      CALL MXCOPYREAL8TOPTR(IUPSI,MXGETPI(PLHS(1)),SIZE)

      RETURN 
      END

C ----- Computational subroutine
      SUBROUTINE COMPUP(GAMMA, SIG2, T1, T2, ROWS, COLS, 
     +   RUPSI, IUPSI)        
	  INTEGER I,J, ROWS, COLS 
      DOUBLE PRECISION RUPSI(ROWS, COLS), IUPSI(ROWS, COLS)
	  DOUBLE PRECISION SIG2, T1(ROWS), T2(COLS)
      DOUBLE PRECISION DEV, REALZ1, REALZ2(COLS) 
	  COMPLEX *16 GAMMA, WOFZ1, WOFZ2(COLS)
      COMPLEX *16 Z1, Z2, UPSI
	  

      DEV = DSQRT(SIG2)
      DO 20, I=1,ROWS 
        DO 10, J=1,COLS             
             IF(I.EQ.1) THEN
                 Z2 = T2(J)/DEV + DEV*GAMMA/2
                 REALZ2(J) = DREAL(Z2)
                 IF(REALZ2(J).GE.0) THEN    
                    CALL WOFZ(DCMPLX(0,1)*Z2, WOFZ2(J))
                 ELSE
                    CALL WOFZ(DCMPLX(0,-1)*Z2, WOFZ2(J))
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
                UPSI = 2*EXP(SIG2*(GAMMA**2)/4 - 
     +              GAMMA*(T1(I)-T2(J))) - EXP(-(T1(I)-T2(J))**2/SIG2
     +         + LOG(WOFZ1)) - EXP(-T2(J)**2/SIG2 - GAMMA*T1(I) 
     +         + LOG(WOFZ2(J)))
             END IF
C		real(z1) < 0 and real(z2)>=0
		   IF(REALZ1. LT. 0. AND .REALZ2(J). GE .0  ) THEN
                UPSI =  EXP(-(T1(I)-T2(J))**2/SIG2
     +         + LOG(WOFZ1)) - EXP(-T2(J)**2/SIG2 - GAMMA*T1(I) 
     +         + LOG(WOFZ2(J)))
             END IF
C		real(z1) >= 0 and real(z2)<0
		   IF(REALZ1. GE. 0. AND .REALZ2(J). LT .0  ) THEN
                UPSI = -  EXP(-(T1(I)-T2(J))**2/SIG2
     +         + LOG(WOFZ1)) + EXP(-T2(J)**2/SIG2 - GAMMA*T1(I) 
     +         + LOG(WOFZ2(J)))
             END IF
C		real(z1) < 0 and real(z2)<0             
		   IF(REALZ1. LT. 0. AND .REALZ2(J). LT .0  ) THEN
                UPSI = - 2*EXP(SIG2*(GAMMA**2)/4 - 
     +              GAMMA*(T1(I)-T2(J))) + EXP(-(T1(I)-T2(J))**2/SIG2
     +         + LOG(WOFZ1)) + EXP(-T2(J)**2/SIG2 - GAMMA*T1(I) 
     +         + LOG(WOFZ2(J)))
             END IF
           RUPSI(I,J) = DREAL(UPSI) 
           IUPSI(I,J) = DIMAG(UPSI)
C            RUPSI(I,J) = T1(1)
C            IUPSI(I,J) = T2(1)
10      CONTINUE
20    CONTINUE
      RETURN  
      END
      
C-----  Complex error function due to Hui et al.

      SUBROUTINE WOFZ(Z, W)       
      DOUBLE PRECISION A0, A1, A2, A3, A4, A5, A6
      DOUBLE PRECISION B0, B1, B2, B3, B4, B5, B6
      COMPLEX*16  Z, T, W
      DOUBLE PRECISION X,Y
      A0 = 122.607931777104326
      A1 = 214.382388694706425 
      A2 = 181.928533092181549 
      A3 = 93.155580458138441 
      A4 = 30.180142196210589 
      A5 = 5.912626209773153 
      A6 = 0.564189583562615 
      B0 = 122.607931773875350
      B1 = 352.730625110963558
      B2 = 457.334478783897737
      B3 = 348.703917719495792
      B4 = 170.354001821091472
      B5 = 53.992906912940207
      B6 = 10.479857114260399
      X = DREAL(Z)
      Y = DIMAG(Z)  
      T = DCMPLX(Y, -X)
      W  = ((((((A6*T + A5)*T + A4)*T + A3)*T + A2)*T + A1)*T+A0)
     +    /(((((((T + B6)*T + B5)*T + B4)*T + B3)*T + B2)*T + B1)*T+B0)
      RETURN
      END
