C     ALGORITHM 599, COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.9, NO. 2,
C     JUN., 1983, P. 255-257.
C**********************************************************************CSUN   10
C**********************************************************************CSUN   20
C**********************************************************************CSUN   30
C                                                                      CSUN   40
C                                                                      CSUN   50
C                                                                      CSUN   60
C     F O R T R A N  SOFTWARE PACKAGE FOR RANDOM NUMBER GENERATION     CSUN   70
C                                                                      CSUN   80
C                                                                      CSUN   90
C                                                                      CSUN  100
C**********************************************************************CSUN  110
C**********************************************************************CSUN  120
C**********************************************************************CSUN  130
C                                                                       SUN  140
C                                                                       SUN  150
C                                                                       SUN  160
C     CONTENTS:                                                         SUN  170
C                                                                       SUN  180
C     1) SUNIF  -  0,1 -UNIFORM DISTRIBUTION                            SUN  190
C                                                                       SUN  200
C     2) SEXPO  - (STANDARD-) EXPONENTIAL DISTRIBUTION                  SUN  210
C                                                                       SUN  220
C     3) SNORM  - (STANDARD-) NORMAL DISTRIBUTION                       SUN  230
C                                                                       SUN  240
C     4) SGAMMA - (STANDARD-) GAMMA DISTRIBUTION                        SUN  250
C                                                                       SUN  260
C     5) KPOISS - POISSON DISTRIBUTION                                  SUN  270
C                                                                       SUN  280
C                                                                       SUN  290
C     THIS PACKAGE CONSTITUTES A FORTRAN-77 DOCUMENTATION OF A SET OF   SUN  300
C     ASSEMBLER FUNCTIONS FOR SAMPLING FROM THE ABOVE DISTRIBUTIONS.    SUN  310
C     ALL ROUTINES MAKE AMPLE USE OF BINARY REPRESENTATIONS OF NUMBERS, SUN  320
C     THEY ARE AMONG THE MOST ACCURATE AND FAST SAMPLING FUNCTIONS      SUN  330
C     KNOWN. THE FORTRAN PROGRAMS BELOW YIELD THE SAME RANDOM NUMBER    SUN  340
C     SEQUENCES AS THE ONES FROM OUR ASSEMBLER PACKAGE, BUT THEY ARE    SUN  350
C     OF COURSE MUCH SLOWER (BY FACTORS 5-8 ON OUR SIEMENS 7760         SUN  360
C     COMPUTER.)                                                        SUN  370
C     THE SET OF ROUTINES WILL ALSO BE ACCEPTABLE TO FORTRAN IV         SUN  380
C     COMPILERS WHICH ALLOW DATA STATEMENTS FOR ARRAYS WITHOUT          SUN  390
C     IMPLICIT DO-LOOPS.                                                SUN  400
C                                                                       SUN  410
C                                                                       SUN  420
C     REMARKS:                                                          SUN  430
C                                                                       SUN  440
C     -  NO CARE IS TAKEN TO ENSURE THAT THE PARAMETER VALUES LIE       SUN  450
C        IN THE ALLOWED RANGE (E.G. A/MU > 0.0 FOR SGAMMA/KPOISS).      SUN  460
C                                                                       SUN  470
C     -  THE PARAMETER 'IR' MUST BE SET TO SOME  4*K+1 > 0  BEFORE      SUN  480
C        THE FIRST CALL OF ANY OF THE GENERATORS. THEREAFTER IR         SUN  490
C        MUST NOT BE ALTERED UNTIL A NEW INITIALIZATION IS DESIRED.     SUN  500
C                                                                       SUN  510
C     -  THE PACKAGE PROVIDES RANDOM DEVIATES OF 6-7 DIGITS ACCURACY.   SUN  520
C        ON MORE ACCURATE COMPUTERS THE CONSTANTS IN SEXPO, SNORM,      SUN  530
C        SGAMMA AND KPOISS OUGHT TO BE ADJUSTED ACCORDING TO LOCAL      SUN  540
C        COMMENTS OR WITH THE AID OF THE TABLES IN THE LITERATURE       SUN  550
C        QUOTED AT THE BEGINNING OF EACH FUNCTION.                      SUN  560
C                                                                       SUN  570
C                                                                       SUN  580
C**********************************************************************CSUN  590
C**********************************************************************CSUN  600
C                                                                      CSUN  610
C                                                                      CSUN  620
C       0 , 1   - U N I F O R M  DISTRIBUTION                          CSUN  630
C                                                                      CSUN  640
C                                                                      CSUN  650
C**********************************************************************CSUN  660
C**********************************************************************CSUN  670
C                                                                      CSUN  680
C     FOR DETAILS SEE:                                                 CSUN  690
C                                                                      CSUN  700
C               AHRENS, J.H., DIETER, U. AND GRUBE, A.                 CSUN  710
C               PSEUDO-RANDOM NUMBERS:  A NEW PROPOSAL                 CSUN  720
C                     FOR THE CHOICE OF MULTIPLICATORS                 CSUN  730
C               COMPUTING, 6 (1970), 121 - 138                         CSUN  740
C                                                                      CSUN  750
C**********************************************************************CSUN  760
C                                                                       SUN  770
      REAL FUNCTION SUNIF(IR)                                           SUN  780
      DOUBLE PRECISION R,FACTOR,TWO28                                   SUN  790
C                                                                       SUN  800
C     FACTOR - INTEGER OF THE FORM 8*K+5 AS CLOSE AS POSSIBLE           SUN  810
C              TO  2**26 * (SQRT(5)-1)/2     (GOLDEN SECTION)           SUN  820
C     TWO28  = 2**28  (I.E. 28 SIGNIFICANT BITS FOR DEVIATES)           SUN  830
C                                                                       SUN  840
      DATA FACTOR /41475557.0D0/, TWO28 /268435456.0D0/                 SUN  850
C                                                                       SUN  860
C     RETURNS SAMPLE U FROM THE  0,1 -UNIFORM DISTRIBUTION              SUN  870
C     BY A MULTIPLICATIVE CONGRUENTIAL GENERATOR OF THE FORM            SUN  880
C        R := R * FACTOR (MOD 1) .                                      SUN  890
C     IN THE FIRST CALL R IS INITIALIZED TO                             SUN  900
C        R := IR / 2**28 ,                                              SUN  910
C     WHERE IR MUST BE OF THE FORM  IR = 4*K+1.                         SUN  920
C     THEN R ASSUMES ALL VALUES  0 < (4*K+1)/2**28 < 1 DURING           SUN  930
C     A FULL PERIOD 2**26 OF SUNIF.                                     SUN  940
C     THE PARAMETER IR IS USED ONLY IN THE FIRST CALL FOR               SUN  950
C     INITIALIZATION OF SUNIF. THEREAFTER (WHEN NEGATIVE)               SUN  960
C     IR BECOMES A DUMMY VARIABLE.                                      SUN  970
C                                                                       SUN  980
      IF (IR .GE. 0) GO TO 1                                            SUN  990
C                                                                       SUN 1000
C     STANDARD CASE:  SAMPLING                                          SUN 1010
C                                                                       SUN 1020
      R=DMOD(R*FACTOR,1.0D0)                                            SUN 1030
      SUNIF=SNGL(R)                                                     SUN 1040
      RETURN                                                            SUN 1050
C                                                                       SUN 1060
C     FIRST CALL: INITIALIZATION                                        SUN 1070
C                                                                       SUN 1080
1     R=DBLE(FLOAT(IR))/TWO28                                           SUN 1090
      R=DMOD(R*FACTOR,1.0D0)                                            SUN 1100
      SUNIF=SNGL(R)                                                     SUN 1110
      IR=-1                                                             SUN 1120
      RETURN                                                            SUN 1130
      END                                                               SUN 1140
C                                                                       SEX   10
C**********************************************************************CSEX   20
C**********************************************************************CSEX   30
C                                                                      CSEX   40
C                                                                      CSEX   50
C     (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                CSEX   60
C                                                                      CSEX   70
C                                                                      CSEX   80
C**********************************************************************CSEX   90
C**********************************************************************CSEX  100
C                                                                      CSEX  110
C     FOR DETAILS SEE:                                                 CSEX  120
C                                                                      CSEX  130
C               AHRENS, J.H. AND DIETER, U.                            CSEX  140
C               COMPUTER METHODS FOR SAMPLING FROM THE                 CSEX  150
C               EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  CSEX  160
C               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               CSEX  170
C                                                                      CSEX  180
C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       CSEX  190
C     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       CSEX  200
C                                                                      CSEX  210
C**********************************************************************CSEX  220
C                                                                       SEX  230
      REAL FUNCTION SEXPO(IR)                                           SEX  240
      DIMENSION Q(8)                                                    SEX  250
      EQUIVALENCE (Q(1),Q1)                                             SEX  260
C                                                                       SEX  270
C     Q(N) = SUM(ALOG(2.0)**K/K])    K=1,..,N ,      THE HIGHEST N      SEX  280
C     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION      SEX  290
C                                                                       SEX  300
      DATA Q/.6931472,.9333737,.9888778,.9984959,                       SEX  310
     ,.9998293,.9999833,.9999986,.9999999/                              SEX  320
C                                                                       SEX  330
   1  A=0.0                                                             SEX  340
      U=SUNIF(IR)                                                       SEX  350
      GO TO 2                                                           SEX  360
   3  A=A+Q1                                                            SEX  370
   2  U=U+U                                                             SEX  380
      IF (U.LE.1.0) GO TO 3                                             SEX  390
   4  U=U-1.0                                                           SEX  400
      IF (U.GT.Q1) GO TO 6                                              SEX  410
   5  SEXPO=A+U                                                         SEX  420
      RETURN                                                            SEX  430
   6  I=1                                                               SEX  440
      USTAR=SUNIF(IR)                                                   SEX  450
      UMIN=USTAR                                                        SEX  460
   7  USTAR=SUNIF(IR)                                                   SEX  470
      IF (USTAR.LT.UMIN) UMIN=USTAR                                     SEX  480
   8  I=I+1                                                             SEX  490
      IF (U.GT.Q(I)) GO TO 7                                            SEX  500
   9  SEXPO=A+UMIN*Q1                                                   SEX  510
      RETURN                                                            SEX  520
      END                                                               SEX  530
C                                                                       SNO   10
C**********************************************************************CSNO   20
C**********************************************************************CSNO   30
C                                                                      CSNO   40
C                                                                      CSNO   50
C     (STANDARD-)  N O R M A L  DISTRIBUTION                           CSNO   60
C                                                                      CSNO   70
C                                                                      CSNO   80
C**********************************************************************CSNO   90
C**********************************************************************CSNO  100
C                                                                      CSNO  110
C     FOR DETAILS SEE:                                                 CSNO  120
C                                                                      CSNO  130
C               AHRENS, J.H. AND DIETER, U.                            CSNO  140
C               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             CSNO  150
C               SAMPLING FROM THE NORMAL DISTRIBUTION.                 CSNO  160
C               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          CSNO  170
C                                                                      CSNO  180
C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  CSNO  190
C     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  CSNO  200
C                                                                      CSNO  210
C**********************************************************************CSNO  220
C                                                                       SNO  230
      REAL FUNCTION SNORM(IR)                                           SNO  240
      DIMENSION A(32),D(31),T(31),H(31)                                 SNO  250
C                                                                       SNO  260
C     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND             SNO  270
C     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE                  SNO  280
C                                                                       SNO  290
      DATA A/0.0,.3917609E-1,.7841241E-1,.1177699,.1573107,             SNO  300
     ,.1970991,.2372021,.2776904,.3186394,.3601299,.4022501,            SNO  310
     ,.4450965,.4887764,.5334097,.5791322,.6260990,.6744898,            SNO  320
     ,.7245144,.7764218,.8305109,.8871466,.9467818,1.009990,            SNO  330
     ,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,            SNO  340
     ,1.675940,1.862732,2.153875/                                       SNO  350
      DATA D/5*0.0,.2636843,.2425085,.2255674,.2116342,.1999243,        SNO  360
     ,.1899108,.1812252,.1736014,.1668419,.1607967,.1553497,            SNO  370
     ,.1504094,.1459026,.1417700,.1379632,.1344418,.1311722,            SNO  380
     ,.1281260,.1252791,.1226109,.1201036,.1177417,.1155119,            SNO  390
     ,.1134023,.1114027,.1095039/                                       SNO  400
      DATA T/.7673828E-3,.2306870E-2,.3860618E-2,.5438454E-2,           SNO  410
     ,.7050699E-2,.8708396E-2,.1042357E-1,.1220953E-1,.1408125E-1,      SNO  420
     ,.1605579E-1,.1815290E-1,.2039573E-1,.2281177E-1,.2543407E-1,      SNO  430
     ,.2830296E-1,.3146822E-1,.3499233E-1,.3895483E-1,.4345878E-1,      SNO  440
     ,.4864035E-1,.5468334E-1,.6184222E-1,.7047983E-1,.8113195E-1,      SNO  450
     ,.9462444E-1,.1123001,.1364980,.1716886,.2276241,.3304980,         SNO  460
     ,.5847031/                                                         SNO  470
      DATA H/.3920617E-1,.3932705E-1,.3950999E-1,.3975703E-1,           SNO  480
     ,.4007093E-1,.4045533E-1,.4091481E-1,.4145507E-1,.4208311E-1,      SNO  490
     ,.4280748E-1,.4363863E-1,.4458932E-1,.4567523E-1,.4691571E-1,      SNO  500
     ,.4833487E-1,.4996298E-1,.5183859E-1,.5401138E-1,.5654656E-1,      SNO  510
     ,.5953130E-1,.6308489E-1,.6737503E-1,.7264544E-1,.7926471E-1,      SNO  520
     ,.8781922E-1,.9930398E-1,.1155599,.1404344,.1836142,.2790016,      SNO  530
     ,.7010474/                                                         SNO  540
C                                                                       SNO  550
   1  U=SUNIF(IR)                                                       SNO  560
      S=0.0                                                             SNO  570
      IF (U.GE.0.5) S=1.0                                               SNO  580
      U=U+U-S                                                           SNO  590
   2  U=32.0*U                                                          SNO  600
      I=INT(U)                                                          SNO  610
      IF (I.EQ.0) GO TO 9                                               SNO  620
C                                                                       SNO  630
C                                START CENTER                           SNO  640
C                                                                       SNO  650
   3  USTAR=U-FLOAT(I)                                                  SNO  660
      AA=A(I)                                                           SNO  670
   4  IF (USTAR.LE.T(I)) GO TO 5                                        SNO  680
      W=(USTAR-T(I))*H(I)                                               SNO  690
C                                                                       SNO  700
C                                EXIT   (BOTH CASES)                    SNO  710
C                                                                       SNO  720
  17  Y=AA+W                                                            SNO  730
      SNORM=Y                                                           SNO  740
      IF (S.EQ.1.0) SNORM=-Y                                            SNO  750
      RETURN                                                            SNO  760
C                                                                       SNO  770
C                                CENTER CONTINUED                       SNO  780
C                                                                       SNO  790
   5  U=SUNIF(IR)                                                       SNO  800
      W=U*(A(I+1)-AA)                                                   SNO  810
      TT=(0.5*W+AA)*W                                                   SNO  820
      GO TO 6                                                           SNO  830
   8  TT=U                                                              SNO  840
      USTAR=SUNIF(IR)                                                   SNO  850
   6  IF (USTAR.GT.TT) GO TO 17                                         SNO  860
   7  U=SUNIF(IR)                                                       SNO  870
      IF (USTAR.GE.U) GO TO 8                                           SNO  880
      USTAR=SUNIF(IR)                                                   SNO  890
      GO TO 4                                                           SNO  900
C                                                                       SNO  910
C                                START TAIL                             SNO  920
C                                                                       SNO  930
   9  I=6                                                               SNO  940
      AA=A(32)                                                          SNO  950
      GO TO 10                                                          SNO  960
  11  AA=AA+D(I)                                                        SNO  970
      I=I+1                                                             SNO  980
  10  U=U+U                                                             SNO  990
      IF (U.LT.1.0) GO TO 11                                            SNO 1000
  12  U=U-1.0                                                           SNO 1010
  13  W=U*D(I)                                                          SNO 1020
      TT=(0.5*W+AA)*W                                                   SNO 1030
      GO TO 14                                                          SNO 1040
  16  TT=U                                                              SNO 1050
  14  USTAR=SUNIF(IR)                                                   SNO 1060
      IF (USTAR.GT.TT) GO TO 17                                         SNO 1070
  15  U=SUNIF(IR)                                                       SNO 1080
      IF (USTAR.GE.U) GO TO 16                                          SNO 1090
      U=SUNIF(IR)                                                       SNO 1100
      GO TO 13                                                          SNO 1110
      END                                                               SNO 1120
C                                                                       SGA   10
C**********************************************************************CSGA   20
C**********************************************************************CSGA   30
C                                                                      CSGA   40
C                                                                      CSGA   50
C     (STANDARD-)  G A M M A  DISTRIBUTION                             CSGA   60
C                                                                      CSGA   70
C                                                                      CSGA   80
C**********************************************************************CSGA   90
C**********************************************************************CSGA  100
C                                                                      CSGA  110
C               PARAMETER  A >= 1.0  ]                                 CSGA  120
C                                                                      CSGA  130
C**********************************************************************CSGA  140
C                                                                      CSGA  150
C     FOR DETAILS SEE:                                                 CSGA  160
C                                                                      CSGA  170
C               AHRENS, J.H. AND DIETER, U.                            CSGA  180
C               GENERATING GAMMA VARIATES BY A                         CSGA  190
C               MODIFIED REJECTION TECHNIQUE.                          CSGA  200
C               COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  CSGA  210
C                                                                      CSGA  220
C     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     CSGA  230
C                                 (STRAIGHTFORWARD IMPLEMENTATION)     CSGA  240
C                                                                      CSGA  250
C**********************************************************************CSGA  260
C                                                                      CSGA  270
C               PARAMETER  0.0 < A < 1.0  ]                            CSGA  280
C                                                                      CSGA  290
C**********************************************************************CSGA  300
C                                                                      CSGA  310
C     FOR DETAILS SEE:                                                 CSGA  320
C                                                                      CSGA  330
C               AHRENS, J.H. AND DIETER, U.                            CSGA  340
C               COMPUTER METHODS FOR SAMPLING FROM GAMMA,              CSGA  350
C               BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              CSGA  360
C               COMPUTING, 12 (1974), 223 - 246.                       CSGA  370
C                                                                      CSGA  380
C     (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)    CSGA  390
C                                                                      CSGA  400
C**********************************************************************CSGA  410
C                                                                       SGA  420
      REAL FUNCTION SGAMMA(IR,A)                                        SGA  430
C                                                                       SGA  440
C     INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR         SGA  450
C             A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION    SGA  460
C     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION           SGA  470
C                                                                       SGA  480
C     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))                    SGA  490
C     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)             SGA  500
C     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)                 SGA  510
C                                                                       SGA  520
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7 /.04166669,.02083148,                   SGA  530
     ,.00801191,.00144121,-.00007388,.00024511,.00024240/               SGA  540
      DATA A1,A2,A3,A4,A5,A6,A7 /.3333333,-.2500030,                    SGA  550
     ,.2000062,-.1662921,.1423657,-.1367177,.1233795/                   SGA  560
      DATA E1,E2,E3,E4,E5 /1.,.4999897,.1668290,.0407753,.0102930/      SGA  570
C                                                                       SGA  580
C     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"                  SGA  590
C     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380                SGA  600
C                                                                       SGA  610
      DATA AA /0.0/, AAA /0.0/, SQRT32 /5.656854/                       SGA  620
C                                                                       SGA  630
      IF (A .EQ. AA) GO TO 1                                            SGA  640
      IF (A .LT. 1.0) GO TO 12                                          SGA  650
C                                                                       SGA  660
C     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED               SGA  670
C                                                                       SGA  680
      AA=A                                                              SGA  690
      S2=A-0.5                                                          SGA  700
      S=SQRT(S2)                                                        SGA  710
      D=SQRT32-12.0*S                                                   SGA  720
C                                                                       SGA  730
C     STEP  2:  T=STANDARD NORMAL DEVIATE,                              SGA  740
C               X=(S,1/2)-NORMAL DEVIATE.                               SGA  750
C               IMMEDIATE ACCEPTANCE (I)                                SGA  760
C                                                                       SGA  770
   1  T=SNORM(IR)                                                       SGA  780
      X=S+0.5*T                                                         SGA  790
      SGAMMA=X*X                                                        SGA  800
      IF (T .GE. 0.0) RETURN                                            SGA  810
C                                                                       SGA  820
C     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)          SGA  830
C                                                                       SGA  840
      U=SUNIF(IR)                                                       SGA  850
      IF (D*U .LE. T*T*T) RETURN                                        SGA  860
C                                                                       SGA  870
C     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY                SGA  880
C                                                                       SGA  890
      IF (A .EQ. AAA) GO TO 4                                           SGA  900
      AAA=A                                                             SGA  910
      R=1.0/A                                                           SGA  920
      Q0=((((((Q7*R+Q6)*R+Q5)*R+Q4)*R+Q3)*R+Q2)*R+Q1)*R                 SGA  930
C                                                                       SGA  940
C               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A          SGA  950
C               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND          SGA  960
C               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS             SGA  970
C                                                                       SGA  980
      IF (A .LE. 3.686) GO TO 3                                         SGA  990
      IF (A .LE. 13.022) GO TO 2                                        SGA 1000
C                                                                       SGA 1010
C               CASE 3:  A .GT. 13.022                                  SGA 1020
C                                                                       SGA 1030
      B=1.77                                                            SGA 1040
      SI=.75                                                            SGA 1050
      C=.1515/S                                                         SGA 1060
      GO TO 4                                                           SGA 1070
C                                                                       SGA 1080
C               CASE 2:  3.686 .LT. A .LE. 13.022                       SGA 1090
C                                                                       SGA 1100
   2  B=1.654+.0076*S2                                                  SGA 1110
      SI=1.68/S+.275                                                    SGA 1120
      C=.062/S+.024                                                     SGA 1130
      GO TO 4                                                           SGA 1140
C                                                                       SGA 1150
C               CASE 1:  A .LE. 3.686                                   SGA 1160
C                                                                       SGA 1170
   3  B=.463+S-.178*S2                                                  SGA 1180
      SI=1.235                                                          SGA 1190
      C=.195/S-.079+.016*S                                              SGA 1200
C                                                                       SGA 1210
C     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE                      SGA 1220
C                                                                       SGA 1230
   4  IF (X .LE. 0.0) GO TO 7                                           SGA 1240
C                                                                       SGA 1250
C     STEP  6:  CALCULATION OF V AND QUOTIENT Q                         SGA 1260
C                                                                       SGA 1270
      V=T/(S+S)                                                         SGA 1280
      IF (ABS(V) .LE. 0.25) GO TO 5                                     SGA 1290
      Q=Q0-S*T+0.25*T*T+(S2+S2)*ALOG(1.0+V)                             SGA 1300
      GO TO 6                                                           SGA 1310
   5  Q=Q0+0.5*T*T*((((((A7*V+A6)*V+A5)*V+A4)*V+A3)*V+A2)*V+A1)*V       SGA 1320
C                                                                       SGA 1330
C     STEP  7:  QUOTIENT ACCEPTANCE (Q)                                 SGA 1340
C                                                                       SGA 1350
   6  IF (ALOG(1.0-U) .LE. Q) RETURN                                    SGA 1360
C                                                                       SGA 1370
C     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE                          SGA 1380
C               U= 0,1 -UNIFORM DEVIATE                                 SGA 1390
C               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE            SGA 1400
C                                                                       SGA 1410
   7  E=SEXPO(IR)                                                       SGA 1420
      U=SUNIF(IR)                                                       SGA 1430
      U=U+U-1.0                                                         SGA 1440
      T=B+SIGN(SI*E,U)                                                  SGA 1450
C                                                                       SGA 1460
C     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719           SGA 1470
C                                                                       SGA 1480
      IF (T .LT. (-.7187449)) GO TO 7                                   SGA 1490
C                                                                       SGA 1500
C     STEP 10:  CALCULATION OF V AND QUOTIENT Q                         SGA 1510
C                                                                       SGA 1520
      V=T/(S+S)                                                         SGA 1530
      IF (ABS(V) .LE. 0.25) GO TO 8                                     SGA 1540
      Q=Q0-S*T+0.25*T*T+(S2+S2)*ALOG(1.0+V)                             SGA 1550
      GO TO 9                                                           SGA 1560
   8  Q=Q0+0.5*T*T*((((((A7*V+A6)*V+A5)*V+A4)*V+A3)*V+A2)*V+A1)*V       SGA 1570
C                                                                       SGA 1580
C     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)     SGA 1590
C                                                                       SGA 1600
   9  IF (Q .LE. 0.0) GO TO 7                                           SGA 1610
      IF (Q .LE. 0.5) GO TO 10                                          SGA 1620
      W=EXP(Q)-1.0                                                      SGA 1630
      GO TO 11                                                          SGA 1640
  10  W=((((E5*Q+E4)*Q+E3)*Q+E2)*Q+E1)*Q                                SGA 1650
C                                                                       SGA 1660
C               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8                SGA 1670
C                                                                       SGA 1680
  11  IF (C*ABS(U) .GT. W*EXP(E-0.5*T*T)) GO TO 7                       SGA 1690
      X=S+0.5*T                                                         SGA 1700
      SGAMMA=X*X                                                        SGA 1710
      RETURN                                                            SGA 1720
C                                                                       SGA 1730
C     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))    SGA 1740
C                                                                       SGA 1750
  12  AA=0.0                                                            SGA 1760
      B=1.0+.3678794*A                                                  SGA 1770
  13  P=B*SUNIF(IR)                                                     SGA 1780
      IF (P .GE. 1.0) GO TO 14                                          SGA 1790
      SGAMMA=EXP(ALOG(P)/A)                                             SGA 1800
      IF (SEXPO(IR) .LT. SGAMMA) GO TO 13                               SGA 1810
      RETURN                                                            SGA 1820
  14  SGAMMA=-ALOG((B-P)/A)                                             SGA 1830
      IF (SEXPO(IR) .LT. (1.0-A)*ALOG(SGAMMA)) GO TO 13                 SGA 1840
      RETURN                                                            SGA 1850
      END                                                               SGA 1860
C                                                                       KPO   10
C**********************************************************************CKPO   20
C**********************************************************************CKPO   30
C                                                                      CKPO   40
C                                                                      CKPO   50
C     P O I S S O N  DISTRIBUTION                                      CKPO   60
C                                                                      CKPO   70
C                                                                      CKPO   80
C**********************************************************************CKPO   90
C**********************************************************************CKPO  100
C                                                                      CKPO  110
C     FOR DETAILS SEE:                                                 CKPO  120
C                                                                      CKPO  130
C               AHRENS, J.H. AND DIETER, U.                            CKPO  140
C               COMPUTER GENERATION OF POISSON DEVIATES                CKPO  150
C               FROM MODIFIED NORMAL DISTRIBUTIONS.                    CKPO  160
C               ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179. CKPO  170
C                                                                      CKPO  180
C     (SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)  CKPO  190
C                                                                      CKPO  200
C**********************************************************************CKPO  210
C                                                                       KPO  220
      INTEGER FUNCTION KPOISS(IR,MU)                                    KPO  230
C                                                                       KPO  240
C     INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR         KPO  250
C             MU=MEAN MU OF THE POISSON DISTRIBUTION                    KPO  260
C     OUTPUT: KPOISS=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION          KPO  270
C                                                                       KPO  280
      REAL MU, MUPREV, MUOLD                                            KPO  290
C                                                                       KPO  300
C     MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR B.    KPO  310
C     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT            KPO  320
C     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL            KPO  330
C                                                                       KPO  340
      DIMENSION FACT(10), PP(35)                                        KPO  350
      DATA MUPREV,MUOLD /0.,0./                                         KPO  360
      DATA A0,A1,A2,A3,A4,A5,A6,A7 /-.5,.3333333,-.2500068,             KPO  370
     ,.2000118,-.1661269,.1421878,-.1384794,.1250060/                   KPO  380
      DATA FACT /1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880./        KPO  390
C                                                                       KPO  400
C     SEPARATION OF CASES A AND B                                       KPO  410
C                                                                       KPO  420
      IF (MU .EQ. MUPREV) GO TO 1                                       KPO  430
      IF (MU .LT. 10.0) GO TO 12                                        KPO  440
C                                                                       KPO  450
C     C A S E  A. (RECALCULATION OF S,D,L IF MU HAS CHANGED)            KPO  460
C                                                                       KPO  470
      MUPREV=MU                                                         KPO  480
      S=SQRT(MU)                                                        KPO  490
      D=6.0*MU*MU                                                       KPO  500
C                                                                       KPO  510
C             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL   KPO  520
C             PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)   KPO  530
C             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .             KPO  540
C                                                                       KPO  550
      L=IFIX(MU-1.1484)                                                 KPO  560
C                                                                       KPO  570
C     STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE     KPO  580
C                                                                       KPO  590
   1  G=MU+S*SNORM(IR)                                                  KPO  600
      IF (G .LT. 0.0) GO TO 2                                           KPO  610
      KPOISS=IFIX(G)                                                    KPO  620
C                                                                       KPO  630
C     STEP I. IMMEDIATE ACCEPTANCE IF KPOISS IS LARGE ENOUGH            KPO  640
C                                                                       KPO  650
      IF (KPOISS .GE. L) RETURN                                         KPO  660
C                                                                       KPO  670
C     STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U         KPO  680
C                                                                       KPO  690
      FK=FLOAT(KPOISS)                                                  KPO  700
      DIFMUK=MU-FK                                                      KPO  710
      U=SUNIF(IR)                                                       KPO  720
      IF (D*U .GE. DIFMUK*DIFMUK*DIFMUK) RETURN                         KPO  730
C                                                                       KPO  740
C     STEP P. PREPARATIONS FOR STEPS Q AND H.                           KPO  750
C             (RECALCULATIONS OF PARAMETERS IF NECESSARY)               KPO  760
C             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7. KPO  770
C             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE KPO  780
C             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.   KPO  790
C             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION. KPO  800
C                                                                       KPO  810
   2  IF (MU .EQ. MUOLD) GO TO 3                                        KPO  820
      MUOLD=MU                                                          KPO  830
      OMEGA=.3989423/S                                                  KPO  840
      B1=.4166667E-1/MU                                                 KPO  850
      B2=.3*B1*B1                                                       KPO  860
      C3=.1428571*B1*B2                                                 KPO  870
      C2=B2-15.*C3                                                      KPO  880
      C1=B1-6.*B2+45.*C3                                                KPO  890
      C0=1.-B1+3.*B2-15.*C3                                             KPO  900
      C=.1069/MU                                                        KPO  910
   3  IF (G .LT. 0.0) GO TO 5                                           KPO  920
C                                                                       KPO  930
C             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)     KPO  940
C                                                                       KPO  950
      KFLAG=0                                                           KPO  960
      GO TO 7                                                           KPO  970
C                                                                       KPO  980
C     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)                           KPO  990
C                                                                       KPO 1000
   4  IF (FY-U*FY .LE. PY*EXP(PX-FX)) RETURN                            KPO 1010
C                                                                       KPO 1020
C     STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL   KPO 1030
C             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'             KPO 1040
C             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)           KPO 1050
C                                                                       KPO 1060
   5  E=SEXPO(IR)                                                       KPO 1070
      U=SUNIF(IR)                                                       KPO 1080
      U=U+U-1.0                                                         KPO 1090
      T=1.8+SIGN(E,U)                                                   KPO 1100
      IF (T .LE. (-.6744)) GO TO 5                                      KPO 1110
      KPOISS=IFIX(MU+S*T)                                               KPO 1120
      FK=FLOAT(KPOISS)                                                  KPO 1130
      DIFMUK=MU-FK                                                      KPO 1140
C                                                                       KPO 1150
C             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)     KPO 1160
C                                                                       KPO 1170
      KFLAG=1                                                           KPO 1180
      GO TO 7                                                           KPO 1190
C                                                                       KPO 1200
C     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)               KPO 1210
C                                                                       KPO 1220
   6  IF (C*ABS(U) .GT. PY*EXP(PX+E)-FY*EXP(FX+E)) GO TO 5              KPO 1230
      RETURN                                                            KPO 1240
C                                                                       KPO 1250
C     STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.               KPO 1260
C             CASE KPOISS .LT. 10 USES FACTORIALS FROM TABLE FACT       KPO 1270
C                                                                       KPO 1280
   7  IF (KPOISS .GE. 10) GO TO 8                                       KPO 1290
      PX=-MU                                                            KPO 1300
      PY=MU**KPOISS/FACT(KPOISS+1)                                      KPO 1310
      GO TO 11                                                          KPO 1320
C                                                                       KPO 1330
C             CASE KPOISS .GE. 10 USES POLYNOMIAL APPROXIMATION         KPO 1340
C             A0-A7 FOR ACCURACY WHEN ADVISABLE                         KPO 1350
C             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)                KPO 1360
C                                                                       KPO 1370
   8  DEL=.8333333E-1/FK                                                KPO 1380
      DEL=DEL-4.8*DEL*DEL*DEL                                           KPO 1390
      V=DIFMUK/FK                                                       KPO 1400
      IF (ABS(V) .LE. 0.25) GO TO 9                                     KPO 1410
      PX=FK*ALOG(1.0+V)-DIFMUK-DEL                                      KPO 1420
      GO TO 10                                                          KPO 1430
   9  PX=FK*V*V*(((((((A7*V+A6)*V+A5)*V+A4)*V+A3)*V+A2)*V+A1)*V+A0)-DEL KPO 1440
  10  PY=.3989423/SQRT(FK)                                              KPO 1450
  11  X=(0.5-DIFMUK)/S                                                  KPO 1460
      XX=X*X                                                            KPO 1470
      FX=-0.5*XX                                                        KPO 1480
      FY=OMEGA*(((C3*XX+C2)*XX+C1)*XX+C0)                               KPO 1490
      IF (KFLAG) 4,4,6                                                  KPO 1500
C                                                                       KPO 1510
C     C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)       KPO 1520
C                                                                       KPO 1530
  12  MUPREV=0.0                                                        KPO 1540
      IF (MU .EQ. MUOLD) GO TO 13                                       KPO 1550
      MUOLD=MU                                                          KPO 1560
      M=MAX0(1,IFIX(MU))                                                KPO 1570
      L=0                                                               KPO 1580
      P=EXP(-MU)                                                        KPO 1590
      Q=P                                                               KPO 1600
      P0=P                                                              KPO 1610
C                                                                       KPO 1620
C     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD                       KPO 1630
C                                                                       KPO 1640
  13  U=SUNIF(IR)                                                       KPO 1650
      KPOISS=0                                                          KPO 1660
      IF (U .LE. P0) RETURN                                             KPO 1670
C                                                                       KPO 1680
C     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE               KPO 1690
C             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES              KPO 1700
C             (0.458=PP(9) FOR MU=10)                                   KPO 1710
C                                                                       KPO 1720
      IF (L .EQ. 0) GO TO 15                                            KPO 1730
      J=1                                                               KPO 1740
      IF (U .GT. 0.458) J=MIN0(L,M)                                     KPO 1750
      DO 14 K=J,L                                                       KPO 1760
      IF (U .LE. PP(K)) GO TO 18                                        KPO 1770
  14  CONTINUE                                                          KPO 1780
      IF (L .EQ. 35) GO TO 13                                           KPO 1790
C                                                                       KPO 1800
C     STEP C. CREATION OF NEW POISSON PROBABILITIES P                   KPO 1810
C             AND THEIR CUMULATIVES Q=PP(K)                             KPO 1820
C                                                                       KPO 1830
  15  L=L+1                                                             KPO 1840
      DO 16 K=L,35                                                      KPO 1850
      P=P*MU/FLOAT(K)                                                   KPO 1860
      Q=Q+P                                                             KPO 1870
      PP(K)=Q                                                           KPO 1880
      IF (U .LE. Q) GO TO 17                                            KPO 1890
  16  CONTINUE                                                          KPO 1900
      L=35                                                              KPO 1910
      GO TO 13                                                          KPO 1920
  17  L=K                                                               KPO 1930
  18  KPOISS=K                                                          KPO 1940
      RETURN                                                            KPO 1950
      END                                                               KPO 1960
C**********************************************************************CMAN   10
C**********************************************************************CMAN   20
C**********************************************************************CMAN   30
C                                                                      CMAN   40
C                                                                      CMAN   50
C                                                                      CMAN   60
C     AUTOMATIC TEST DRIVER                                            CMAN   70
C     FOR                                                              CMAN   80
C     RANDOM NUMBER PACKAGE AHRENS/DIETER/KOHRT                        CMAN   90
C     (THIS DRIVER MAY CAUSE MANY EXP() UNDERFLOWS.)                    MAN  100
C                                                                      CMAN  110
C                                                                      CMAN  120
C                                                                      CMAN  130
C**********************************************************************CMAN  140
C**********************************************************************CMAN  150
C**********************************************************************CMAN  160
C     PROGRAM TEST                                                      MAN  170
      REAL MU                                                           MAN  180
      DIMENSION VPAR4(22),VPAR5(24),SAMPLE(10000),II(5)                 MAN  190
C                                                                       MAN  200
C     VPAR4 - VECTOR OF PARAMETER VALUES FOR CASE 4 : SGAMMA            MAN  210
C                                                                       MAN  220
      DATA VPAR4 /.0001,.25,.5,.75,.9999,1.,1.5,2.,                     MAN  230
     ,3.,4.,5.,7.,10.,15.,20.,30.,50.,100.,1000.,                       MAN  240
     ,10000.,100000.,1000000./                                          MAN  250
C                                                                       MAN  260
C     VPAR5 - VECTOR OF PARAMETER VALUES FOR CASE 5 : KPOISS            MAN  270
C                                                                       MAN  280
      DATA VPAR5 /.0001,1.,2.,5.,9.99,10.,                              MAN  290
     ,12.,15.,20.,25.,30.,40.,50.,75.,100.,150.,                        MAN  300
     ,200.,500.,1000.,2000.,5000.,1.E4,1.E5,1.E6/                       MAN  310
C                                                                       MAN  320
C     II - NUMBER OF RUNS FOR EACH DISTRIBUTION                         MAN  330
C                                                                       MAN  340
      DATA II /1,1,1,22,24/                                             MAN  350
C                                                                       MAN  360
C     FORMAT STATEMENTS                                                 MAN  370
C                                                                       MAN  380
   1  FORMAT(' ',/,' ',/,' LISTING OF TRIAL RUNS',                      MAN  390
     ,' FOR RANDOM NUMBER PACKAGE AHRENS/DIETER/KOHRT',                 MAN  400
C    ,/,' ===============================',                             MAN  410
     ,'====================================')                           MAN  420
   2  FORMAT('   FIRST 100 SAMPLES:',/,                                 MAN  430
     ,'   ..................',/,' ',/,(5E15.6))                         MAN  440
   3  FORMAT(' ',/,' ',/,'   TEST DATA:',                               MAN  450
     ,'     (]BASED ON 10000 SAMPLES])',                                MAN  460
     ,/,'   ..........',/,20X,'MEAN',11X,                               MAN  470
     ,'STD.DEV.',7X,'SKEWNESS',                                         MAN  480
     ,7X,'EXCESS',/,' ',/,                                              MAN  490
     ,'   TRUE VALUES:',4E15.6,/,                                       MAN  500
     ,'   SAMPLE DATA:',4E15.6,/,' ')                                   MAN  510
   4  FORMAT(' ',/,' 1.)   0,1 -UNIFORM DISTRIBUTION:',                 MAN  520
     ,/,' ********************************')                            MAN  530
   5  FORMAT(' ',/,' 2.)  (STANDARD-) EXPONENTIAL DISTRIBUTION:',       MAN  540
     ,/,' ******************************************')                  MAN  550
   6  FORMAT(' ',/,' 3.)  (STANDARD-) NORMAL DISTRIBUTION:',            MAN  560
     ,/,' ******************************************')                  MAN  570
   7  FORMAT(' ',/,' 4.)  (STANDARD-) GAMMA-(A) DISTRIBUTION:',         MAN  580
     ,/,' ****************************************')                    MAN  590
   8  FORMAT(' ',/,' 5.)  POISSON-(MU) DISTRIBUTION:',                  MAN  600
     ,/,' *******************************',/,' ',                       MAN  610
     ,/,'   (INTEGER SAMPLES ARE DISPLAYED AS REALS])',                 MAN  620
     ,/,' ',/,' ')                                                      MAN  630
   9  FORMAT(43X,'    GAMMA-(A):  A =',E13.6,                           MAN  640
     ,/,43X,'    ----------------------------')                         MAN  650
  10  FORMAT(43X,'POISSON-(MU):  MU =',E13.6,                           MAN  660
     ,/,43X,'--------------------------------')                         MAN  670
  11  FORMAT(' ',/,' ')                                                 MAN  680
C                                                                       MAN  690
C     DEFINE OUTPUT UNIT NUMBER                                         MAN  700
C                                                                       MAN  710
      NOUT=10                                                           MAN  720
C                                                                       MAN  730
C     OUTPUT: MAIN HEADING                                              MAN  740
C                                                                       MAN  750
      WRITE (NOUT,1)                                                    MAN  760
C                                                                       MAN  770
C     TRIAL RUNS FOR 5 DIFFERENT CASES:                                 MAN  780
C                                                                       MAN  790
C       NDIS=1 :   0,1 -UNIFORM DISTRIBUTION                            MAN  800
C       NDIS=2 :  (STANDARD-) EXPONENTIAL DISTRIBUTION                  MAN  810
C       NDIS=3 :  (STANDARD-) NORMAL DISTRIBUTION                       MAN  820
C       NDIS=4 :  (STANDARD-) GAMMA-(A) DISTRIBUTION                    MAN  830
C       NDIS=5 :  POISSON-(MU) DISTRIBUTION                             MAN  840
C                                                                       MAN  850
      DO 27 NDIS=1,5                                                    MAN  860
      NRUN=II(NDIS)                                                     MAN  870
C                                                                       MAN  880
C     OUTPUT: CASE HEADING                                              MAN  890
C                                                                       MAN  900
      WRITE (NOUT,11)                                                   MAN  910
      IF (NDIS .EQ. 1) WRITE (NOUT,4)                                   MAN  920
      IF (NDIS .EQ. 2) WRITE (NOUT,5)                                   MAN  930
      IF (NDIS .EQ. 3) WRITE (NOUT,6)                                   MAN  940
      IF (NDIS .EQ. 4) WRITE (NOUT,7)                                   MAN  950
      IF (NDIS .EQ. 5) WRITE (NOUT,8)                                   MAN  960
C                                                                       MAN  970
C     EACH CASE: ONE RUN FOR EVERY PARAMETER VALUE                      MAN  980
C                                                                       MAN  990
      DO 26 NPAR=1,NRUN                                                 MAN 1000
C                                                                       MAN 1010
C     CASE 4 AND 5: SET PARAMETER VALUES ACCORDING TO DATA VECTOR       MAN 1020
C                                                                       MAN 1030
      IF (NDIS .EQ. 4) A=VPAR4(NPAR)                                    MAN 1040
      IF (NDIS .EQ. 5) MU=VPAR5(NPAR)                                   MAN 1050
C                                                                       MAN 1060
C     SEED FOR UNIFORM RANDOM NUMBER GENERATOR IS INITIALIZED TO 4*0+1  MAN 1070
C                                                                       MAN 1080
      IR=1                                                              MAN 1090
C                                                                       MAN 1100
C     EACH CASE SEPARATELY: SAMPLING AND TEST DATA                      MAN 1110
C        T2 - STANDARD DEVIATION (=SQRT(VARIANCE)),                     MAN 1120
C        T1 - MEAN,   T3 - SKEWNESS,   T4 - EXCESS.                     MAN 1130
C                                                                       MAN 1140
      GO TO (12,14,16,18,20) , NDIS                                     MAN 1150
C                                                                       MAN 1160
C     CASE 1 :   0,1 -UNIFORM DISTRIBUTION                              MAN 1170
C                                                                       MAN 1180
  12  DO 13 I=1,10000                                                   MAN 1190
  13  SAMPLE(I)=SUNIF(IR)                                               MAN 1200
      T1=0.5                                                            MAN 1210
      T2=1.0/SQRT(12.0)                                                 MAN 1220
      T3=0.0                                                            MAN 1230
      T4=-1.2                                                           MAN 1240
      GO TO 23                                                          MAN 1250
C                                                                       MAN 1260
C     (STANDARD-) EXPONENTIAL DISTRIBUTION                              MAN 1270
C                                                                       MAN 1280
  14  DO 15 I=1,10000                                                   MAN 1290
  15  SAMPLE(I)=SEXPO(IR)                                               MAN 1300
      T1=1.0                                                            MAN 1310
      T2=1.0                                                            MAN 1320
      T3=2.0                                                            MAN 1330
      T4=6.0                                                            MAN 1340
      GO TO 23                                                          MAN 1350
C                                                                       MAN 1360
C     (STANDARD-) NORMAL DISTRIBUTION                                   MAN 1370
C                                                                       MAN 1380
  16  DO 17 I=1,10000                                                   MAN 1390
  17  SAMPLE(I)=SNORM(IR)                                               MAN 1400
      T1=0.0                                                            MAN 1410
      T2=1.0                                                            MAN 1420
      T3=0.0                                                            MAN 1430
      T4=0.0                                                            MAN 1440
      GO TO 23                                                          MAN 1450
C                                                                       MAN 1460
C     (STANDARD-) GAMMA-(A) DISTRIBUTION                                MAN 1470
C                                                                       MAN 1480
  18  DO 19 I=1,10000                                                   MAN 1490
  19  SAMPLE(I)=SGAMMA(IR,A)                                            MAN 1500
      T1=A                                                              MAN 1510
      T2=SQRT(A)                                                        MAN 1520
      T3=2.0/T2                                                         MAN 1530
      T4=6.0/A                                                          MAN 1540
      GO TO 22                                                          MAN 1550
C                                                                       MAN 1560
C     POISSON-(MU) DISTRIBUTION                                         MAN 1570
C                                                                       MAN 1580
  20  DO 21 I=1,10000                                                   MAN 1590
  21  SAMPLE(I)=FLOAT(KPOISS(IR,MU))                                    MAN 1600
      T1=MU                                                             MAN 1610
      T2=SQRT(MU)                                                       MAN 1620
      T3=1.0/T2                                                         MAN 1630
      T4=1.0/MU                                                         MAN 1640
C                                                                       MAN 1650
C     CASE 4 AND 5:  OUTPUT: PARAMETER VALUE                            MAN 1660
C                                                                       MAN 1670
  22  IF (NPAR .NE. 1) WRITE (NOUT,11)                                  MAN 1680
      IF (NDIS .EQ. 4) WRITE (NOUT, 9) A                                MAN 1690
      IF (NDIS .EQ. 5) WRITE (NOUT,10) MU                               MAN 1700
C                                                                       MAN 1710
  23  CONTINUE                                                          MAN 1720
C                                                                       MAN 1730
C     OUTPUT : FIRST 100 RANDOM DEVIATES FOR EACH RUN                   MAN 1740
C              (INTEGER SAMPLES ARE DISPLAYED AS REALS])                MAN 1750
C                                                                       MAN 1760
      IF (NDIS .LE. 3) WRITE (NOUT,11)                                  MAN 1770
      WRITE (NOUT,2) (SAMPLE(I),I=1,100)                                MAN 1780
C                                                                       MAN 1790
C     EVALUATION OF SAMPLE MEAN:    E1 - (1/N)*SUM(SAMPLE(I))           MAN 1800
C                                                                       MAN 1810
      S1=0.0                                                            MAN 1820
      DO 24 I=1,10000                                                   MAN 1830
  24  S1=S1+SAMPLE(I)                                                   MAN 1840
      E1=S1/10000.0                                                     MAN 1850
C                                                                       MAN 1860
C     EVALUATION OF FURTHER SAMPLE ESTIMATES :                          MAN 1870
C                                                                       MAN 1880
C     SK       - (1/N)*SUM((SAMPLE(I)-E1)**K)                           MAN 1890
C                SAMPLE CENTRAL MOMENTS  (K=2,3,4)                      MAN 1900
C                WITH RESPECT TO SAMPLE MEAN                            MAN 1910
C                                                                       MAN 1920
C     E2       - SQRT(S2)                                               MAN 1930
C                SAMPLE STANDARD DEVIATION                              MAN 1940
C                (=SQRT(SAMPLE VARIANCE))                               MAN 1950
C                                                                       MAN 1960
C     E3       - S3/S2**(3/2)                                           MAN 1970
C                SAMPLE SKEWNESS                                        MAN 1980
C                                                                       MAN 1990
C     E4       - S4/S2**2-3                                             MAN 2000
C                SAMPLE EXCESS                                          MAN 2010
C                                                                       MAN 2020
      S2=0.0                                                            MAN 2030
      S3=0.0                                                            MAN 2040
      S4=0.0                                                            MAN 2050
C                                                                       MAN 2060
      DO 25 I=1,10000                                                   MAN 2070
C                                                                       MAN 2080
      X1=SAMPLE(I)-E1                                                   MAN 2090
      X2=X1*X1                                                          MAN 2100
      X3=X2*X1                                                          MAN 2110
      X4=X2*X2                                                          MAN 2120
C                                                                       MAN 2130
      S2=S2+X2                                                          MAN 2140
      S3=S3+X3                                                          MAN 2150
      S4=S4+X4                                                          MAN 2160
C                                                                       MAN 2170
  25  CONTINUE                                                          MAN 2180
C                                                                       MAN 2190
      S2=S2/10000.0                                                     MAN 2200
      S3=S3/10000.0                                                     MAN 2210
      S4=S4/10000.0                                                     MAN 2220
C                                                                       MAN 2230
      E2=SQRT(S2)                                                       MAN 2240
      E3=S3/SQRT(S2*S2*S2)                                              MAN 2250
      E4=S4/(S2*S2)-3.0                                                 MAN 2260
C                                                                       MAN 2270
C     END OF EVALUATION                                                 MAN 2280
C     OUTPUT: CHARACTERISTIC DATA                                       MAN 2290
C                                                                       MAN 2300
      WRITE (NOUT,3) T1,T2,T3,T4,E1,E2,E3,E4                            MAN 2310
C                                                                       MAN 2320
C     END OF PARAMETER LOOP                                             MAN 2330
C                                                                       MAN 2340
  26  CONTINUE                                                          MAN 2350
C                                                                       MAN 2360
C     END OF DISTRIBUTION LOOP                                          MAN 2370
C                                                                       MAN 2380
  27  CONTINUE                                                          MAN 2390
C                                                                       MAN 2400
C     END OF PROGRAM                                                    MAN 2410
C                                                                       MAN 2420
      STOP                                                              MAN 2430
      END                                                               MAN 2440
