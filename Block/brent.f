      SUBROUTINE BRENT(TDOT,IP,DTU)
*
*       Brent's method for root finding.
*       --------------------------------
*
*       Coded by Denis & Marina Ryabov.
*       -------------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      REAL*8  TDOT(11)
      REAL*8  M
      REAL*8  COEFF(11)
      PARAMETER (C3=1.0/3.0D0,C6=1.0/6.0D0,C7=1.0/7.0D0,C9=1.0/9.0D0,
     &           C11=1.0/11.0D0)
      DATA COEFF/1D0, 0.5D0, C3, 0.25D0, 0.2D0, C6, C7, 0.125D0, C9,
     &           0.1D0, C11/
*
*       Solving t(DTU)=STEP equation
*
      STEPI1 = STEP(2*IP-1)
*
*       Include improved tolerance near small pericentre (R/a < 1.0D-02).
      IF (R(IP)*ABS(H(IP))/BODY(N+IP).LT.5.0D-03) THEN
          TOL = 1.0D-10
      ELSE
          TOL = 1.0D-08
      END IF
*
*       Include possible slow-down factor.
      IMOD = KSLOW(IP)
      IF (IMOD.GT.1) THEN
          ZMOD = FLOAT(ISLOW(IMOD))
          STEPI1 = STEPI1/ZMOD
          TOL = TOL/ZMOD  ! Not sure it is necessary
      END IF
*
*       Bracketing
      B = DTU
      FB = 0.0
      DO K = 11,1,-1
          FB = (FB + TDOT(K))*B*COEFF(K)
      END DO
      FB = FB - STEPI1

      IF (FB.EQ.0.0) GO TO 30  ! A wonder!!!

      IF (FB.GT.0.0) THEN
          A = 0.0
          FA = -STEPI1
      ELSE
*           Expand interval
          A = B
          FA = FB
          DO WHILE (FB.LT.0.0)
              B = B * 1.5
              FB = 0.0
              DO K = 11,1,-1
                  FB = (FB + TDOT(K))*B*COEFF(K)
              END DO
              FB = FB - STEPI1
          END DO
      END IF
*
      IF (ABS(FB).GT.1.0D-20) GO TO 88
      WRITE (6,86) TIME, FB
   86 FORMAT (' SMALL     T FB ',F10.5,1P,E10.2)
      STOP
   88 CONTINUE
*       Brent's method
      FC = FB
      DO 10 I = 1, 50
          IF (FB*FC.GT.0.0) THEN
              C = A
              FC = FA
              D = B-A
              EE = D
          END IF
          IF (ABS(FC).LT.ABS(FB)) THEN
              A = B
              B = C
              C = A
              FA = FB
              FB = FC
              FC = FA
          END IF
          M = 0.5*(C-B)
          TOL1 = 4.45D-16*ABS(B)
          IF (B.GT.0.0.AND.(ABS(FB).LT.TOL.OR.ABS(M).LE.TOL1)) GO TO 30
          IF (ABS(EE).LT.TOL1.OR.ABS(FA).LE.ABS(FB)) THEN
              D = M
              EE = M
          ELSE
              S = FB/FA
              IF (A.EQ.C) THEN
                  P = 2.0*M*S
                  Q = 1.0-S
              ELSE
                  Q = FA/FC
                  RR = FB/FC
                  P = S*(2.0*M*Q*(Q-RR)-(B-A)*(RR-1.0))
                  Q = (Q-1.0)*(RR-1.0)*(S-1.0)
              END IF
              IF (P.GT.0.0) THEN
                  Q = -Q
              ELSE
                  P = -P
              END IF
              IF (2.0*P.LT.MIN(3.0*M*Q-ABS(TOL1*Q),ABS(EE*Q))) THEN
                  EE = D
                  D = P/Q
              ELSE
                  D = M
                  EE = M
              END IF
          END IF
          A = B
          FA = FB
          IF (ABS(D).GT.TOL1) THEN
              B = B + D
          ELSE
              B = B + SIGN(TOL1,M)
          END IF
          FB = 0.0
          DO K = 11,1,-1
              FB = (FB + TDOT(K))*B*COEFF(K)
          END DO
          FB = FB - STEPI1
   10 CONTINUE
*
      WRITE (6,20)  NSTEPU, FB, B, DTAU(IP)
   20 FORMAT (' DANGER BRENT!    # F DTU DTAU ',I11,1P,4E10.2)
      STOP
*
   30 CONTINUE
      DTU = B
      RETURN
*
      END
