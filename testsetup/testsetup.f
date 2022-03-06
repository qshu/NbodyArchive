      PROGRAM TESTSETUP
*
*
*       Generation of initial coordinates & velocities.
*       -----------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  A(8)
      REAL*4  RAN2
      REAL*8 CMR(3),CMRDOT(3)
*
      PRINT*,' Enter N,NRAND:'
      READ*,N,NRAND
*
      IDUM1 = -NRAND
*       Choose between uniform density and Plummer model.
      KDUM = IDUM1
      TWOPI = 8.0D0*ATAN(1.D0)
*
      DO 1 I = 1,N
 1    BODY(I) = 1.D0/DBLE(N)
*
      ZMASS =1.D0
*
*       Initialize centre of mass terms.
      DO 25 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   25 CONTINUE
*
*       Generate initial conditions from Plummer model (A & A 37, 183).
      DO 40 I = 1,N
   30     A(1) = RAN2(KDUM)
          IF (A(1).LT.1.0D-10) GO TO 30
          RI = (A(1)**(-0.6666667) - 1.0)**(-0.5)
*       Reject distant particles.
          IF (RI.GT.10.0) GO TO 30
*
          A(2) = RAN2(KDUM)
          A(3) = RAN2(KDUM)
          X(3,I) = (1.0 - 2.0*A(2))*RI
          X(1,I) = SQRT(RI**2 - X(3,I)**2)*COS(TWOPI*A(3))
          X(2,I) = SQRT(RI**2 - X(3,I)**2)*SIN(TWOPI*A(3))
   32     A(4) = RAN2(KDUM)
          A(5) = RAN2(KDUM)
          A(6) = A(4)**2*(1.0 - A(4)**2)**3.5
          IF (0.1*A(5).GT.A(6)) GO TO 32
*
          A(8) = A(4)*SQRT(2.0)/(1.0 + RI**2)**0.25
          A(6) = RAN2(KDUM)
          A(7) = RAN2(KDUM)
          XDOT(3,I) = (1.0 - 2.0*A(6))*A(8)
          XDOT(1,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*COS(TWOPI*A(7))
          XDOT(2,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*SIN(TWOPI*A(7))
*
*       Accumulate c.m. terms.
          DO 35 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   35     CONTINUE
   40 CONTINUE
*
*       Scale coordinates & velocities to analytical expectation values.
      SX = 1.5*TWOPI/16.0
      SV = SQRT(ZMASS/SX)
      DO 50 I = 1,N
          DO 45 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
              X(K,I) = SX*X(K,I)
              XDOT(K,I) = SV*XDOT(K,I)
   45     CONTINUE
   50 CONTINUE
*
*       Save random number sequence in COMMON for future use.
*

      PRINT*,' Going to call ENERGY '
      CALL FLUSH(6)
*
      CALL ENERGY(ZKIN,POT)
*
      PRINT*,' Old Method ZKIN,POT,Q=',ZKIN,POT,ZKIN/POT
*
      DO 999 I=1,N
      WRITE(10,1111)BODY(I),(X(K,I),K=1,3),(XDOT(K,I),K=1,3),PHIDBL(I)
 1111 FORMAT(1X,1P,8(1X,D15.7))
 999  CONTINUE
*
*     test of different procedure
*
      IDUM1 = -NRAND
      KDUM = IDUM1
*
      RMAX = 10.D0
*
*       Generate initial conditions using rejection technique.
      DO 140 I = 1,N
*
 130  RR = RMAX*RAN2(KDUM)
      RR2 = RR*RR
*
      RHOR = 3.D0/(1.D0+RR2)**2.5D0
      FXR = RHOR*RR2
*
      XRAN = RAN2(KDUM)
      IF(XRAN.LT.FXR)THEN
*
      RI = RR
*
          A(2) = RAN2(KDUM)
          A(3) = RAN2(KDUM)
          X(3,I) = (1.0 - 2.0*A(2))*RI
          X(1,I) = SQRT(RI**2 - X(3,I)**2)*COS(TWOPI*A(3))
          X(2,I) = SQRT(RI**2 - X(3,I)**2)*SIN(TWOPI*A(3))
  132     A(4) = RAN2(KDUM)
          A(5) = RAN2(KDUM)
          A(6) = A(4)**2*(1.0 - A(4)**2)**3.5
          IF (0.1*A(5).GT.A(6)) GO TO 132
*
          A(8) = A(4)*SQRT(2.0)/(1.0 + RI**2)**0.25
          A(6) = RAN2(KDUM)
          A(7) = RAN2(KDUM)
          XDOT(3,I) = (1.0 - 2.0*A(6))*A(8)
          XDOT(1,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*COS(TWOPI*A(7))
          XDOT(2,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*SIN(TWOPI*A(7))
*       
*       Accumulate c.m. terms.
          DO 135 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
  135     CONTINUE
*
          ELSE
*
          IREJ = IREJ + 1
          GOTO 130
*
          END IF
*
  140 CONTINUE
*
*       Scale coordinates & velocities to analytical expectation values.
      SX = 1.5*TWOPI/16.0
      SV = SQRT(ZMASS/SX)
      DO 150 I = 1,N
          DO 145 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
              X(K,I) = SX*X(K,I)
              XDOT(K,I) = SV*XDOT(K,I)
  145     CONTINUE
  150 CONTINUE
*
*
      PRINT*,' Going to call ENERGY '
      CALL FLUSH(6)
*
      CALL ENERGY(ZKIN,POT)
*
      PRINT*,' New Method ZKIN,POT,Q=',ZKIN,POT,ZKIN/POT
*
      DO 998 I=1,N
      WRITE(11,1111)BODY(I),(X(K,I),K=1,3),(XDOT(K,I),K=1,3),PHIDBL(I)
 998  CONTINUE
*
      STOP
*
      END
