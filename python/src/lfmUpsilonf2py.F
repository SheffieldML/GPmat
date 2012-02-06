C ----- Computational subroutine to compute Upsilon
C------ when it is a vector
      SUBROUTINE COMPUPVEC(GAMMA, SIG2, T1, ROWS,
     +   UPSI)
      INTEGER I,  ROWS
      REAL *8 SIG2, T1(ROWS), DEV, REALZ1, REALZ2
      COMPLEX *16 GAMMA, WOFZ1, WOFZ2
      COMPLEX *16 Z1, Z2, UPSI(ROWS)
cf2py intent(in) :: GAMMA, SIG2, T1, ROWS
cf2py intent(out) :: UPSI
      DEV = DSQRT(SIG2)
      Z2 = DEV*GAMMA/2
      REALZ2 = DREAL(Z2)

      IF(REALZ2 .GE.0) THEN
         CALL WOFZ(DCMPLX(0,1)*Z2, WOFZ2)
      ELSE
         CALL WOFZ(DCMPLX(0,-1)*Z2, WOFZ2)
      END IF

      DO 20, I=1,ROWS
             Z1 = T1(I)/DEV - DEV*GAMMA/2
             REALZ1 = DREAL(Z1)
             IF(REALZ1.GE.0) THEN
                CALL WOFZ(DCMPLX(0,1)*Z1, WOFZ1)
             ELSE
                CALL WOFZ(DCMPLX(0,-1)*Z1, WOFZ1)
             END IF
C		real(z1) >= 0 and real(z2)>=0
             IF(REALZ1. GE. 0. AND .REALZ2. GE .0  ) THEN
                UPSI(I) = 2*CDEXP(SIG2*(GAMMA**2)/4 -
     +              GAMMA*(T1(I))) - CDEXP(-(T1(I))**2/SIG2
     +         + CDLOG(WOFZ1)) - CDEXP( - GAMMA*T1(I)
     +         + CDLOG(WOFZ2))
             END IF
C		real(z1) < 0 and real(z2)>=0
		   IF(REALZ1. LT. 0. AND .REALZ2. GE .0  ) THEN
                UPSI(I) =  CDEXP(-(T1(I))**2/SIG2
     +         + CDLOG(WOFZ1)) - CDEXP(- GAMMA*T1(I)
     +         + CDLOG(WOFZ2))
             END IF
C		real(z1) >= 0 and real(z2)<0
		   IF(REALZ1. GE. 0. AND .REALZ2. LT .0  ) THEN
                UPSI(I) = -  CDEXP(-(T1(I))**2/SIG2
     +         + CDLOG(WOFZ1)) + CDEXP(- GAMMA*T1(I)
     +         + CDLOG(WOFZ2))
             END IF
C		real(z1) < 0 and real(z2)<0
		   IF(REALZ1. LT. 0. AND .REALZ2. LT .0  ) THEN
                UPSI(I) = - 2*CDEXP(SIG2*(GAMMA**2)/4 -
     +              GAMMA*(T1(I))) + CDEXP(-(T1(I))**2/SIG2
     +         + CDLOG(WOFZ1)) + CDEXP( - GAMMA*T1(I)
     +         + CDLOG(WOFZ2))
             END IF
20    CONTINUE
      RETURN
      END

C ----- Computational subroutine to compute Upsilon
C------ when it is a matrix
      SUBROUTINE COMPUPMAT(GAMMA, SIG2, T1, T2, ROWS, COLS,
     +   UPSI)
      INTEGER I,J, ROWS, COLS
      REAL *8 SIG2, T1(ROWS), T2(COLS)
      REAL *8 DEV, REALZ1, REALZ2(COLS)
      COMPLEX *16 GAMMA, WOFZ1, WOFZ2(COLS)
      COMPLEX *16 Z1, Z2, UPSI(ROWS, COLS)
cf2py intent(in) :: GAMMA, SIG2, T1, T2, ROWS, COLS
cf2py intent(out) :: UPSI
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
                UPSI(I,J) = 2*CDEXP(SIG2*(GAMMA**2)/4 -
     +              GAMMA*(T1(I)-T2(J))) - CDEXP(-(T1(I)-T2(J))**2/SIG2
     +         + CDLOG(WOFZ1)) - CDEXP(-T2(J)**2/SIG2 - GAMMA*T1(I)
     +         + CDLOG(WOFZ2(J)))
             END IF
C		real(z1) < 0 and real(z2)>=0
		   IF(REALZ1. LT. 0. AND .REALZ2(J). GE .0  ) THEN
                UPSI(I,J) =  CDEXP(-(T1(I)-T2(J))**2/SIG2
     +         + CDLOG(WOFZ1)) - CDEXP(-T2(J)**2/SIG2 - GAMMA*T1(I)
     +         + CDLOG(WOFZ2(J)))
             END IF
C		real(z1) >= 0 and real(z2)<0
		   IF(REALZ1. GE. 0. AND .REALZ2(J). LT .0  ) THEN
                UPSI(I,J) = -  CDEXP(-(T1(I)-T2(J))**2/SIG2
     +         + CDLOG(WOFZ1)) + CDEXP(-T2(J)**2/SIG2 - GAMMA*T1(I)
     +         + CDLOG(WOFZ2(J)))
             END IF
C		real(z1) < 0 and real(z2)<0
		   IF(REALZ1. LT. 0. AND .REALZ2(J). LT .0  ) THEN
                UPSI(I,J) = - 2*CDEXP(SIG2*(GAMMA**2)/4 -
     +              GAMMA*(T1(I)-T2(J))) + CDEXP(-(T1(I)-T2(J))**2/SIG2
     +         + CDLOG(WOFZ1)) + CDEXP(-T2(J)**2/SIG2 - GAMMA*T1(I)
     +         + CDLOG(WOFZ2(J)))
             END IF
10      CONTINUE
20    CONTINUE
      RETURN
      END

C ----- Computational subroutine to compute the gradient of Upsilon
C------ when it is a vector
      SUBROUTINE GRADUPVEC(GAMMA, SIG2, T1, ROWS, GUPSI)
      INTEGER I, ROWS
      REAL *8 PI, SQRTPI
      REAL *8 SIG2, T1(ROWS), SIG, REALZ1, REALZ2
      COMPLEX *16 GAMMA, WOFZ1, WOFZ2
      COMPLEX *16 Z1, Z2, GUPSI(ROWS)
cf2py intent(in) :: GAMMA, SIG2, T1, ROWS
cf2py intent(out) :: GUPSI
      PI = 3.141592653589793D0
      SQRTPI = DSQRT(PI)
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
                GUPSI(I) = CDEXP(SIG2*(GAMMA**2)/4 - GAMMA*(T1(I))
     +                + CDLOG(SIG2*GAMMA - 2*(T1(I))))
     +                - CDEXP(-(T1(I))**2/SIG2 + DLOG(SIG)
     +                + CDLOG(1/SQRTPI - Z1*WOFZ1))
     +                + CDEXP(- GAMMA*T1(I)
     +                + CDLOG(T1(I)*WOFZ2
     +                + SIG*(1/SQRTPI - Z2*WOFZ2)))
             END IF
C		real(z1) < 0 and real(z2)>=0
		     IF(REALZ1. LT. 0. AND .REALZ2. GE .0  ) THEN
                GUPSI(I) = - CDEXP(-(T1(I))**2/SIG2 + DLOG(SIG)
     +                + CDLOG(1/SQRTPI + Z1*WOFZ1))
     +                + CDEXP( - GAMMA*T1(I)
     +                + CDLOG(T1(I)*WOFZ2
     +                + SIG*(1/SQRTPI - Z2*WOFZ2)))
             END IF
C		real(z1) >= 0 and real(z2)<0
		   IF(REALZ1. GE. 0. AND .REALZ2. LT .0  ) THEN
                GUPSI(I) = - CDEXP(-(T1(I))**2/SIG2 + DLOG(SIG)
     +                + CDLOG(1/SQRTPI - Z1*WOFZ1))
     +                - CDEXP(- GAMMA*T1(I)
     +                + CDLOG(T1(I)*WOFZ2
     +                - SIG*(1/SQRTPI + Z2*WOFZ2)))
             END IF
C		real(z1) < 0 and real(z2)<0
		   IF(REALZ1. LT. 0. AND .REALZ2. LT .0  ) THEN
                GUPSI(I) = - CDEXP(SIG2*(GAMMA**2)/4 - GAMMA*(T1(I))
     +                + CDLOG(SIG2*GAMMA - 2*(T1(I))))
     +                - CDEXP(-(T1(I))**2/SIG2 + DLOG(SIG)
     +                + CDLOG(1/SQRTPI + Z1*WOFZ1))
     +                - CDEXP(- GAMMA*T1(I)
     +                + CDLOG(T1(I)*WOFZ2
     +                - SIG*(1/SQRTPI + Z2*WOFZ2)))
           END IF
20    CONTINUE
      RETURN
      END

C ----- Computational subroutine to compute the gradient of Upsilon
C------ when it is a matrix
      SUBROUTINE GRADUPMAT(GAMMA, SIG2, T1, T2, ROWS, COLS,
     +   GUPSI)
      INTEGER I,J, ROWS, COLS
      REAL *8 PI, SQRTPI
      REAL *8 SIG2, T1(ROWS), T2(COLS)
      REAL *8 SIG, REALZ1, REALZ2(COLS)
      COMPLEX *16 GAMMA, WOFZ1, WOFZ2(COLS)
      COMPLEX *16 Z1, Z2(COLS), GUPSI(ROWS, COLS)
cf2py intent(in) :: GAMMA, SIG2, T1, T2, ROWS, COLS
cf2py intent(out) :: GUPSI
      PI = 3.141592653589793D0
      SIG = DSQRT(SIG2)
      SQRTPI = DSQRT(PI)
      DO 20, I=1,ROWS
        DO 10, J=1,COLS
             IF(I.EQ.1) THEN
                 Z2(J) = T2(J)/SIG + SIG*GAMMA/2
                 REALZ2(J) = DREAL(Z2(J))
                 IF(REALZ2(J).GE.0) THEN
                    CALL WOFZ(DCMPLX(0,1)*Z2(J), WOFZ2(J))
                 ELSE
                    CALL WOFZ(DCMPLX(0,-1)*Z2(J), WOFZ2(J))
                 END IF
             END IF
             Z1 = (T1(I)-T2(J))/SIG - SIG*GAMMA/2
             REALZ1 = DREAL(Z1)
             IF(REALZ1.GE.0) THEN
                CALL WOFZ(DCMPLX(0,1)*Z1, WOFZ1)
             ELSE
                CALL WOFZ(DCMPLX(0,-1)*Z1, WOFZ1)
             END IF
C		real(z1) >= 0 and real(z2)>=0
             IF(REALZ1. GE. 0. AND .REALZ2(J). GE .0  ) THEN
                GUPSI(I,J) = CDEXP(SIG2*(GAMMA**2)/4
     +                - GAMMA*(T1(I)-T2(J))
     +                + CDLOG(SIG2*GAMMA - 2*(T1(I)-T2(J))))
     +                - CDEXP(-(T1(I)-T2(J))**2/SIG2 + DLOG(SIG)
     +                + CDLOG(1/SQRTPI - Z1*WOFZ1))
     +                + CDEXP(-T2(J)**2/SIG2 - GAMMA*T1(I)
     +                + CDLOG(T1(I)*WOFZ2(J)
     +                + SIG*(1/SQRTPI - Z2(J)*WOFZ2(J))))
             END IF
C		real(z1) < 0 and real(z2)>=0
		     IF(REALZ1. LT. 0. AND .REALZ2(J). GE .0  ) THEN
                GUPSI(I,J) = - CDEXP(-(T1(I)-T2(J))**2/SIG2 + DLOG(SIG)
     +                + CDLOG(1/SQRTPI + Z1*WOFZ1))
     +                + CDEXP(-T2(J)**2/SIG2 - GAMMA*T1(I)
     +                + CDLOG(T1(I)*WOFZ2(J)
     +                + SIG*(1/SQRTPI - Z2(J)*WOFZ2(J))))
             END IF
C		real(z1) >= 0 and real(z2)<0
		   IF(REALZ1. GE. 0. AND .REALZ2(J). LT .0  ) THEN
                GUPSI(I,J) = - CDEXP(-(T1(I)-T2(J))**2/SIG2 + DLOG(SIG)
     +                + CDLOG(1/SQRTPI - Z1*WOFZ1))
     +                - CDEXP(-T2(J)**2/SIG2 - GAMMA*T1(I)
     +                + CDLOG(T1(I)*WOFZ2(J)
     +                - SIG*(1/SQRTPI + Z2(J)*WOFZ2(J))))
             END IF
C		real(z1) < 0 and real(z2)<0
		   IF(REALZ1. LT. 0. AND .REALZ2(J). LT .0  ) THEN
                GUPSI(I,J) = - CDEXP(SIG2*(GAMMA**2)/4
     +                - GAMMA*(T1(I)-T2(J))
     +                + CDLOG(SIG2*GAMMA - 2*(T1(I)-T2(J))) )
     +                - CDEXP(-(T1(I)-T2(J))**2/SIG2 + DLOG(SIG)
     +                + CDLOG(1/SQRTPI + Z1*WOFZ1))
     +                - CDEXP(-T2(J)**2/SIG2 - GAMMA*T1(I)
     +                + CDLOG(T1(I)*WOFZ2(J)
     +                - SIG*(1/SQRTPI + Z2(J)*WOFZ2(J))))
           END IF
10      CONTINUE
20    CONTINUE
      RETURN
      END

C ----- Special computational subroutine to compute the gradient of Upsilon
C------ wtr to Sigma, when Upsilon is a vector
      SUBROUTINE GRADSIGUPVEC(GAMMA, SIG2, T1, ROWS,
     +   GRAD)
      INTEGER I, ROWS
      REAL *8 PI, SQRTPI
      REAL *8 SIG2, T1(ROWS)
      REAL *8 DEV, REALZ1, REALZ2
      COMPLEX *16 GAMMA, WOFZ1, WOFZ2
      COMPLEX *16 Z1, Z2, GRAD(ROWS)
cf2py intent(in) :: GAMMA, SIG2, T1, ROWS
cf2py intent(out) :: GRAD
      PI = 3.141592653589793D0
      SQRTPI = DSQRT(PI)
      DEV    = DSQRT(SIG2)
      Z2     =  DEV*GAMMA/2
      REALZ2 = DREAL(Z2)
      IF(REALZ2.GE.0) THEN
         CALL WOFZ(DCMPLX(0,1)*Z2, WOFZ2)
      ELSE
         CALL WOFZ(DCMPLX(0,-1)*Z2, WOFZ2)
      END IF
      DO 20, I=1,ROWS
             Z1 = (T1(I))/DEV - DEV*GAMMA/2
             REALZ1 = DREAL(Z1)
             IF(REALZ1.GE.0) THEN
                CALL WOFZ(DCMPLX(0,1)*Z1, WOFZ1)
             ELSE
                CALL WOFZ(DCMPLX(0,-1)*Z1, WOFZ1)
             END IF
C		real(z1) >= 0 and real(z2)>=0
             IF(REALZ1. GE. 0. AND .REALZ2. GE .0  ) THEN
                GRAD(I) = DEV*(GAMMA**2)*CDEXP(SIG2*(GAMMA**2)/4
     +                  - GAMMA*(T1(I)))
     +                  - 2*CDEXP(-((T1(I))**2/SIG2)
     +                  + CDLOG(((T1(I))**2)*WOFZ1/(DEV**3)
     +                  + ((T1(I))/SIG2+GAMMA/2)
     +                  * (1/SQRTPI-Z1*WOFZ1)))
     +                  - 2*CDEXP(-GAMMA*T1(I) + CDLOG((-GAMMA/2)*
     +                  ((1/SQRTPI)-Z2*WOFZ2)))
              END IF
C		real(z1) < 0 and real(z2)>=0
              IF(REALZ1. LT. 0. AND .REALZ2. GE .0  ) THEN
                GRAD(I) =   2*CDEXP(-(T1(I))**2/SIG2
     +                  + CDLOG(((T1(I))**2)*WOFZ1/(DEV**3)
     +                  - ((T1(I))/SIG2+GAMMA/2)
     +                  * (1/SQRTPI+Z1*WOFZ1)))
     +                  - 2*CDEXP(-GAMMA*T1(I)
     +                  + CDLOG((-GAMMA/2)*(1/SQRTPI-Z2*WOFZ2)))
              END IF
C		real(z1) >= 0 and real(z2)<0
              IF(REALZ1. GE. 0. AND .REALZ2. LT .0  ) THEN
                GRAD(I) =  - 2*CDEXP(-(T1(I))**2/SIG2
     +                  + CDLOG(((T1(I))**2)*WOFZ1/(DEV**3)
     +                  + ((T1(I))/SIG2+GAMMA/2)
     +                  * (1/SQRTPI-Z1*WOFZ1)))
     +                  + 2*CDEXP(-GAMMA*T1(I)
     +                  + CDLOG(GAMMA/2*(1/SQRTPI+Z2*WOFZ2)))
              END IF
C		real(z1) < 0 and real(z2)<0
              IF(REALZ1. LT. 0. AND .REALZ2. LT .0  ) THEN
                GRAD(I) = - DEV*(GAMMA**2)*CDEXP(SIG2*(GAMMA**2)/4
     +                  - GAMMA*(T1(I)))
     +                  + 2*CDEXP(-(T1(I))**2/SIG2
     +                  + CDLOG(((T1(I))**2)*WOFZ1/(DEV**3)
     +                  - ((T1(I))/SIG2+GAMMA/2)
     +                  * (1/SQRTPI+Z1*WOFZ1)))
     +                  + 2*CDEXP(-GAMMA*T1(I)
     +                  + CDLOG(GAMMA/2*((1/SQRTPI)+Z2*WOFZ2)))
              END IF
20    CONTINUE
      RETURN
      END

C ----- Special computational subroutine to compute the gradient of Upsilon
C------ wtr to Sigma, when Upsilon is a matrix
      SUBROUTINE GRADSIGUPMAT(GAMMA, SIG2, T1, T2, ROWS, COLS,
     +   GRAD)
      INTEGER I,J, ROWS, COLS
      REAL *8 PI, SQRTPI
      REAL *8 SIG2, T1(ROWS), T2(COLS)
      REAL *8 DEV, REALZ1, REALZ2(COLS)
      COMPLEX *16 GAMMA, WOFZ1, WOFZ2(COLS)
      COMPLEX *16 Z1, Z2(COLS), GRAD(ROWS, COLS)
cf2py intent(in) :: GAMMA, SIG2, T1, T2, ROWS, COLS
cf2py intent(out) :: GRAD
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
                GRAD(I,J) = DEV*(GAMMA**2)*CDEXP(SIG2*(GAMMA**2)/4
     +                  - GAMMA*(T1(I)-T2(J)))
     +                  - 2*CDEXP(-((T1(I)-T2(J))**2/SIG2)
     +                  + CDLOG(((T1(I)-T2(J))**2)*WOFZ1/(DEV**3)
     +                  + ((T1(I)-T2(J))/SIG2+GAMMA/2)
     +                  * (1/SQRTPI-Z1*WOFZ1)))
     +                  - 2*CDEXP(-((T2(J)**2)/SIG2)-GAMMA*T1(I)
     +                  + CDLOG((T2(J)**2)*WOFZ2(J)/(DEV**3)
     +                  + (T2(J)/SIG2-GAMMA/2)*
     +                  ((1/SQRTPI)-Z2(J)*WOFZ2(J))))
              END IF
C		real(z1) < 0 and real(z2)>=0
              IF(REALZ1. LT. 0. AND .REALZ2(J). GE .0  ) THEN
                 GRAD(I,J) =   2*CDEXP(-(T1(I)-T2(J))**2/SIG2
     +                  + CDLOG(((T1(I)-T2(J))**2)*WOFZ1/(DEV**3)
     +                  - ((T1(I)-T2(J))/SIG2+GAMMA/2)
     +                  * (1/SQRTPI+Z1*WOFZ1)))
     +                  - 2*CDEXP(-(T2(J)**2)/SIG2-GAMMA*T1(I)
     +                  + CDLOG((T2(J)**2)*WOFZ2(J)/(DEV**3)
     +                  + (T2(J)/SIG2-GAMMA/2)*
     +                  (1/SQRTPI-Z2(J)*WOFZ2(J))))
              END IF
C		real(z1) >= 0 and real(z2)<0
              IF(REALZ1. GE. 0. AND .REALZ2(J). LT .0  ) THEN
                 GRAD(I,J) =  - 2*CDEXP(-(T1(I)-T2(J))**2/SIG2
     +                  + CDLOG(((T1(I)-T2(J))**2)*WOFZ1/(DEV**3)
     +                  + ((T1(I)-T2(J))/SIG2+GAMMA/2)
     +                  * (1/SQRTPI-Z1*WOFZ1)))
     +                  + 2*CDEXP(-(T2(J)**2)/SIG2-GAMMA*T1(I)
     +                  + CDLOG((T2(J)**2)*WOFZ2(J)/(DEV**3)
     +                  - (T2(J)/SIG2-GAMMA/2)*
     +                  (1/SQRTPI+Z2(J)*WOFZ2(J))))
              END IF
C		real(z1) < 0 and real(z2)<0
              IF(REALZ1. LT. 0. AND .REALZ2(J). LT .0  ) THEN
                 GRAD(I,J) = - DEV*(GAMMA**2)*CDEXP(SIG2*(GAMMA**2)/4
     +                  - GAMMA*(T1(I)-T2(J)))
     +                  + 2*CDEXP(-(T1(I)-T2(J))**2/SIG2
     +                  + CDLOG(((T1(I)-T2(J))**2)*WOFZ1/(DEV**3)
     +                  - ((T1(I)-T2(J))/SIG2+GAMMA/2)
     +                  * (1/SQRTPI+Z1*WOFZ1)))
     +                  + 2*CDEXP(-((T2(J)**2)/SIG2-GAMMA*T1(I))
     +                  + CDLOG((T2(J)**2)*WOFZ2(J)/(DEV**3)
     +                  - (T2(J)/SIG2-GAMMA/2)*
     +                  ((1/SQRTPI)+Z2(J)*WOFZ2(J))))
              END IF
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
cf2py intent(in) :: Z
cf2py intent(out) :: W
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
100   FLAG = .TRUE.
      RETURN
      END
