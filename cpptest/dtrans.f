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
