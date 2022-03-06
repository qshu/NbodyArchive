      SUBROUTINE CHSTAB(ITERM)
*
*
*       Chain stability test.
*       ---------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      REAL*8  M,MB,MB1,R2(NMX,NMX),XCM(3),VCM(3),XX(3,3),VV(3,3)
      INTEGER  IJ(NMX)
*
*
*       Sort particle separations (I1 & I2 form closest pair).
      CALL R2SORT(IJ,R2)
      I1 = IJ(1)
      I2 = IJ(2)
      I3 = IJ(3)
      MB = M(I1) + M(I2)
      MB1 = MB + M(I3)
*
*       Form output diagnostics.
      VREL2 = 0.0D0
      VREL21 = 0.0D0
      RDOT = 0.0D0
      RDOT3 = 0.0D0
      RI2 = 0.0D0
      DO 10 K = 1,3
          J1 = 3*(I1-1) + K
          J2 = 3*(I2-1) + K
          J3 = 3*(I3-1) + K
          VREL2 = VREL2 + (V(J1) - V(J2))**2
          RDOT = RDOT + (X(J1) - X(J2))*(V(J1) - V(J2))
          XCM(K) = (M(I1)*X(J1) + M(I2)*X(J2))/MB
          VCM(K) = (M(I1)*V(J1) + M(I2)*V(J2))/MB
          RI2 = RI2 + (X(J3) - XCM(K))**2
          VREL21 = VREL21 + (V(J3) - VCM(K))**2
          RDOT3 = RDOT3 + (X(J3) - XCM(K))*(V(J3) - VCM(K))
          XX(K,1) = X(J1)
          XX(K,2) = X(J2)
          XX(K,3) = X(J3)
          VV(K,1) = V(J1)
          VV(K,2) = V(J2)
          VV(K,3) = V(J3)
   10 CONTINUE
*
*       Evaluate orbital elements for inner and outer motion.
      RB = SQRT(R2(I1,I2))
      R3 = SQRT(RI2)
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
      ECC = SQRT((1.0D0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
      SEMI1 = 2.0/R3 - VREL21/MB1
      SEMI1 = 1.0/SEMI1
      ECC1 = SQRT((1.0D0 - R3/SEMI1)**2 + RDOT3**2/(SEMI1*MB1))
*
*       Form hierarchical stability ratio (Eggleton & Kiseleva 1995).
      QL = MB/M(I3)
      Q1 = MAX(M(I2)/M(I1),M(I1)/M(I2))
      Q3 = QL**0.33333
      Q13 = Q1**0.33333
      AR = 1.0 + 3.7/Q3 - 2.2/(1.0 + Q3) + 1.4/Q13*(Q3 - 1.0)/(Q3 + 1.0)
*
      EK = AR*SEMI*(1.0D0 + ECC)
      PMIN = SEMI1*(1.0D0 - ECC1)
*
*       Replace the EK criterion by the MA analytical stability formula.
      Q0 = M(I3)/MB
      IF (ECC1.LT.1.0) THEN
          XFAC = (1.0 + Q0)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
      ELSE
          XFAC = 1.0 + Q0
      END IF
      FE = 1.0
      PCRIT = 2.8*FE*XFAC**0.4*SEMI
*
*       Obtain the inclination.
      CALL INCLIN(XX,VV,XCM,VCM,ALPHA)
*
*       Include fudge factor for inclination effect.
      YFAC = 1.0 - 0.3*ALPHA/180.0
*       Add 1% for perturbation to avoid repeated switching.
      PCRIT = 1.01*YFAC*PCRIT
*
*       Check hierarchical stability condition for bound close pair.
      ITERM = 0
      IF (PMIN.GT.PCRIT.AND.SEMI.GT.0.0.AND.SEMI1.GT.0.0) THEN
          ITERM = -1
          WRITE (6,20)  ECC, ECC1, SEMI, SEMI1, PMIN, PCRIT, EK, ALPHA
   20     FORMAT (' NEW HIARCH    E =',F6.3,'  E1 =',F6.3,
     &                     '  A =',1P,E8.1,'  A1 =',E8.1,'  PM =',E9.2,
     &                     '  PC =',E9.2,'  EK =',E9.2,'  IN =',0P,F6.1)
          RI = SQRT(CM(1)**2 + CM(2)**2 + CM(3)**2)
          EMAX = 0.0
          K = INAME(I3)
          WRITE (81,30)  TIMEC, RI, NAMEC(K), QL, Q1, ECC, ECC1,
     &                   SEMI, SEMI1, PCRIT/PMIN, ALPHA, EMAX
   30     FORMAT (F8.1,F5.1,I6,2F6.2,2F6.3,1P,2E10.2,0P,F5.2,F6.1,F6.3)
          CALL FLUSH(81)
*       Include termination test for wide triple system (exclude ECC1 > 0.9).
      ELSE IF (PMIN.GT.2.5*SEMI*(1.0 + ECC).AND.SEMI.GT.0.0.AND.
     &         ECC1.LT.0.9) THEN
*       Wait for favourable configuration (R > SEMI and well separated).
          IF (RB.GT.SEMI.AND.R3.GT.5.0*RB) THEN
              APO = SEMI*(1.0 + ECC)
              WRITE (6,40)  ECC, ECC1, ALPHA, RB, R3, PCRIT, PMIN, APO
   40         FORMAT (' WIDE CHAIN    E E1 IN RB R3 PC PM APO ',
     &                                2F7.3,F7.1,1P,5E10.2)
              ITERM = -1
          END IF
      END IF
*
      RETURN
*
      END
