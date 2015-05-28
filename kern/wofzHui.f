	  SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
	  IMPLICIT NONE
	  INTEGER PLHS(*), PRHS(*)
	  INTEGER MXGETM, MXGETN
      INTEGER NLHS, NRHS
	  INTEGER M(1), N(1)

	  M(1) = MXGETM(PRHS(1))
	  N(1) = MXGETN(PRHS(1))
	
!	CALL INTERMEDIATE GATEWAY:
	  CALL GATEWAY(PLHS,PRHS,M(1),N(1))

	  END 
!	==============================================================

	  SUBROUTINE GATEWAY(PLHS, PRHS, T, P)
	  IMPLICIT NONE
      INTEGER PLHS(1), PRHS(1)
	  INTEGER MXCREATEDOUBLEMATRIX, MXGETPR, MXGETPI
	  INTEGER MXGETM, MXGETN, MXISCOMPLEX
      INTEGER T, P, I, J, ISCOMPLEX
	  DOUBLE PRECISION A(T,P), B(T,P), C(T,P), D(T,P)

      call mxCopyPtrToReal8(mxGetPr(prhs(1)),A,T*P)
!   Checks if there are imaginary numbers in data
      ISCOMPLEX= MXISCOMPLEX(PRHS(1))               
      IF (ISCOMPLEX.EQ.1) THEN   
	  call mxCopyPtrToReal8(mxGetPi(prhs(1)),B,T*P)
      END IF               

	  plhs(1) = mxCreateDoubleMatrix(T,P,1)

!	  CALL COMPUTATIONAL ROUTINE:
      DO 20, I=1,T 
      DO 10, J=1,P
      IF (ISCOMPLEX.EQ.0) THEN   
	  B(I,J) = 0
      END IF   
      CALL WOFZ(A(I,J),B(I,J),C(I,J),D(I,J))
10    CONTINUE
20    CONTINUE
 	  CALL mxCopyReal8ToPtr(C,mxgetpr(plhs(1)),T*P)
      CALL mxCopyReal8ToPtr(D,mxgetpi(plhs(1)),T*P)
	  END
!	===============================================================
      SUBROUTINE WOFZ (X, Y, U, V)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION A0, A1, A2, A3, A4, A5, A6
      DOUBLE PRECISION B0, B1, B2, B3, B4, B5, B6
      COMPLEX *16  T, W
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
      T  = DCMPLX(Y,-X)
      W  = ((((((A6*T + A5)*T + A4)*T + A3)*T + A2)*T + A1)*T+A0)
     $    /(((((((T + B6)*T + B5)*T + B4)*T + B3)*T + B2)*T + B1)*T+B0)
      U = DREAL(W)
      V = DIMAG(W)
      RETURN
      END
