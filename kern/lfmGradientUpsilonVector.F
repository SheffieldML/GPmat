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
      IF (NRHS .NE. 3) THEN
         CALL MEXERRMSGTXT('Three inputs required.')
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
         CALL MEXERRMSGTXT('Input #4 is not a column vector.')
      ENDIF

C	  
	  M = MXGETM(PRHS(3))
	  
	
C	 CALL INTERMEDIATE GATEWAY:
	  CALL GATEWAY(PLHS,PRHS,M)

	  END 
C ------ GATEWAY SUBROUTINE

	  SUBROUTINE GATEWAY(PLHS, PRHS, ROWS)
	  IMPLICIT NONE
      MWPOINTER PLHS(*), PRHS(*)
	  MWPOINTER MXCREATEDOUBLEMATRIX, MXGETPR, MXGETPI
	  INTEGER MXISCOMPLEX
      INTEGER ROWS,ISCOMPLEX
	  DOUBLE PRECISION SIG2, T1(ROWS), ALPHA, OMEGA
	  DOUBLE PRECISION RGUPSI(ROWS), IGUPSI(ROWS)
      COMPLEX*16 GAMMA

      PLHS(1) = MXCREATEDOUBLEMATRIX(ROWS,1,1)

C ------ Convert MATLAB variables to FORTRAN variables  

C ------ Convert gamma 
      CALL MXCOPYPTRTOREAL8(MXGETPR(PRHS(1)), ALPHA ,1)        
      ISCOMPLEX = MXISCOMPLEX(PRHS(1))       
      IF(ISCOMPLEX.EQ.1) THEN
         CALL MXCOPYPTRTOREAL8(MXGETPI(PRHS(1)),OMEGA,1)          
      ELSE 
         OMEGA = 0
      END IF
      GAMMA = DCMPLX(ALPHA,OMEGA)
C ------ Convert sigma2        
      CALL MXCOPYPTRTOREAL8(MXGETPR(PRHS(2)), SIG2 ,1)
C ------ Convert t1        
      CALL MXCOPYPTRTOREAL8(MXGETPR(PRHS(3)), T1 ,ROWS)

C ------ Call computational subroutine     

      CALL COMPGUP(GAMMA,SIG2,T1, ROWS, RGUPSI, IGUPSI)   

C ------  Convert FORTRAN variables to MATLAB variables  

      CALL MXCOPYREAL8TOPTR(RGUPSI,MXGETPR(PLHS(1)),ROWS)
      CALL MXCOPYREAL8TOPTR(IGUPSI,MXGETPI(PLHS(1)),ROWS)

      RETURN 
      END

C ----- Computational subroutine
      SUBROUTINE COMPGUP(GAMMA, SIG2, T1, ROWS,RGUPSI, IGUPSI)
      PARAMETER(PI = 3.141592653589793D0)   
	  INTEGER I, ROWS 
      DOUBLE PRECISION RGUPSI(ROWS), IGUPSI(ROWS), SQRTPI
	  DOUBLE PRECISION SIG2, T1(ROWS), SIG, REALZ1, REALZ2 
	  COMPLEX *16 GAMMA, WOFZ1, WOFZ2
      COMPLEX *16 Z1, Z2, GUPSI
	  
      SQRTPI = SQRT(PI)
      SIG = DSQRT(SIG2)
      Z2 =  SIG*GAMMA/2
      REALZ2 = DREAL(Z2)
      IF(REALZ2.GE.0) THEN    
         CALL WOFZ(DCMPLX(0,1)*Z2, WOFZ2)
      ELSE
         CALL WOFZ(DCMPLX(0,-1)*Z2, WOFZ2)
      END IF  
      DO 20, I=1,ROWS 
             Z1 = (T1(I))/SIG - SIG*GAMMA/2   
             REALZ1 = DREAL(Z1)
             IF(REALZ1.GE.0) THEN    
                CALL WOFZ(DCMPLX(0,1)*Z1, WOFZ1)
             ELSE
                CALL WOFZ(DCMPLX(0,-1)*Z1, WOFZ1)
             END IF
C		real(z1) >= 0 and real(z2)>=0
             IF(REALZ1. GE. 0. AND .REALZ2. GE .0  ) THEN
                GUPSI = EXP(SIG2*(GAMMA**2)/4 - GAMMA*(T1(I))
     +                + LOG(SIG2*GAMMA - 2*(T1(I))))
     +                - EXP(-(T1(I))**2/SIG2 + LOG(SIG)
     +                + LOG(1/SQRTPI - Z1*WOFZ1))     
     +                + EXP(- GAMMA*T1(I)
     +                + LOG(T1(I)*WOFZ2 
     +                + SIG*(1/SQRTPI - Z2*WOFZ2)))
             END IF
C		real(z1) < 0 and real(z2)>=0
		     IF(REALZ1. LT. 0. AND .REALZ2. GE .0  ) THEN
                GUPSI = - EXP(-(T1(I))**2/SIG2 + LOG(SIG)
     +                + LOG(1/SQRTPI + Z1*WOFZ1))     
     +                + EXP( - GAMMA*T1(I)
     +                + LOG(T1(I)*WOFZ2 
     +                + SIG*(1/SQRTPI - Z2*WOFZ2)))
             END IF
C		real(z1) >= 0 and real(z2)<0
		   IF(REALZ1. GE. 0. AND .REALZ2. LT .0  ) THEN
                GUPSI = - EXP(-(T1(I))**2/SIG2 + LOG(SIG)
     +                + LOG(1/SQRTPI - Z1*WOFZ1))     
     +                - EXP(- GAMMA*T1(I)
     +                + LOG(T1(I)*WOFZ2 
     +                - SIG*(1/SQRTPI + Z2*WOFZ2)))
             END IF
C		real(z1) < 0 and real(z2)<0             
		   IF(REALZ1. LT. 0. AND .REALZ2. LT .0  ) THEN
                GUPSI = - EXP(SIG2*(GAMMA**2)/4 - GAMMA*(T1(I))
     +                + LOG(SIG2*GAMMA - 2*(T1(I))))
     +                - EXP(-(T1(I))**2/SIG2 + LOG(SIG)
     +                + LOG(1/SQRTPI + Z1*WOFZ1))
     +                - EXP(- GAMMA*T1(I)
     +                + LOG(T1(I)*WOFZ2 
     +                - SIG*(1/SQRTPI + Z2*WOFZ2)))
           END IF
           RGUPSI(I) = DREAL(GUPSI) 
           IGUPSI(I) = DIMAG(GUPSI)
20    CONTINUE
      RETURN  
      END
      
C-----  Complex error function due to Hui et al.

      SUBROUTINE WOFZ(Z, W)       
      DOUBLE PRECISION A0, A1, A2, A3, A4, A5, A6
      DOUBLE PRECISION B0, B1, B2, B3, B4, B5, B6
      COMPLEX*16  Z, T, W
      DOUBLE PRECISION X,Y
      A0 = 122.607931777104326D0
      A1 = 214.382388694706425D0 
      A2 = 181.928533092181549D0 
      A3 = 93.155580458138441D0 
      A4 = 30.180142196210589D0 
      A5 = 5.912626209773153D0 
      A6 = 0.564189583562615D0 
      B0 = 122.607931773875350D0
      B1 = 352.730625110963558D0
      B2 = 457.334478783897737D0
      B3 = 348.703917719495792D0
      B4 = 170.354001821091472D0
      B5 = 53.992906912940207D0
      B6 = 10.479857114260399D0
      X = DREAL(Z)
      Y = DIMAG(Z)  
      T = DCMPLX(Y, -X)
      W  = ((((((A6*T + A5)*T + A4)*T + A3)*T + A2)*T + A1)*T+A0)
     +    /(((((((T + B6)*T + B5)*T + B4)*T + B3)*T + B2)*T + B1)*T+B0)
      RETURN
      END
