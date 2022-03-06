      SUBROUTINE BRAKE(IPAIR)
*
*
*       Magnetic braking & gravitational radiation.
*       -------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8 M1,M2
*     REAL*4 RL,Q
*     EXTERNAL RL
      SAVE NDIAG
      DATA NDIAG /0/
*
*
*       Set indices of KS components & c.m.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      I = N + IPAIR
      M1 = BODY(I1)*ZMBAR
      M2 = BODY(I2)*ZMBAR
*
*       Quit if active ROCHE contains WD (TIME > TEV from MDOT).
      IF ((KSTAR(I).EQ.21.OR.KSTAR(I).EQ.23).AND.
     &    MAX(KSTAR(I1),KSTAR(I2)).GE.10) THEN
          GO TO 40
      END IF
*
*       Identify the most magnetically active star (#J2).
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
          IF ((KSTAR(I1).EQ.1.AND.M1.GT.1.3).OR.
     &        (KSTAR(I1).EQ.0.AND.M1.LT.0.3)) THEN
              J2 = I2
          ELSE
              J2 = I1
          END IF
      ELSE
          IF ((KSTAR(I2).EQ.1.AND.M2.GT.1.3).OR.
     &        (KSTAR(I2).EQ.0.AND.M2.LT.0.3)) THEN
              J2 = I1
          ELSE
              J2 = I2
          END IF
      END IF
*
*       Set mass and radius of star #J2 in solar units.
      M2 = BODY(J2)*ZMBAR
      R2 = RADIUS(J2)*SU
*
*       Obtain semi-major axis and eccentricity.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
*
*       Define binary period, solar radius & GMS in cgs units.
      TB = 3.147D+07*YRS*SEMI*SQRT(SEMI/BODY(I))
      RSUN = 6.96D+10
      GMS = GM*ZMBAR*BODY(I)
*
*       Form time derivative of angular momentum (Regos & Tout 1995).
      ZJDOT = 3.8D-30*1.989D+33*M2*RSUN**4*R2**3*(TWOPI/TB)**3
*
*       Determine new semi-major axis from angular momentum relation.
      ZMU = BODY(I1)*BODY(I2)/BODY(I)
      ACM = 3.08D+18*RBAR*SEMI
      ADOT = 2.0*SQRT(ACM/GMS)*ZJDOT/(1.989D+33*ZMBAR*ZMU)
*
*       Define old primary and evaluate time scale for spin-down (yrs).
      TBR = ACM/(3.147D+07*ADOT)
*
*       Evaluate alternative expression derived by J/JDOT (factor 2 longer).
*     J1 = I1
*     IF (J1.EQ.J2) J1 = I2
*     TBRAKE = 4.4D+06*(BODY(J1)/BODY(I))/(BODY(I)*ZMBAR)*(ACM/RSUN)**5/
*    &                                                            R2**3
*
*       Include gravitational radiation (64*GM**3/5*c**5; Ap.J. 418, 147).
*     FAC = 64.0*((6.67D-08*1.989D+33)/3.0D+10)**3/(5.0*9.0D+20)
      AGDOT = 1.23D+27*BODY(I1)*BODY(I2)*BODY(I)*(ZMBAR/ACM)**3
      TGR = ACM/(3.147D+07*AGDOT)
*
*       Suppress magnetic braking for massive MS/low-mass or evolved star.
      IF (((M2.GT.1.3.OR.M2.LT.0.3).AND.KSTAR(J2).LE.2).OR.
     &      KSTAR(J2).GE.10) THEN
          ADOT = 0.0
          TBR = 1.0D+10
      END IF
*
*       Combine effects but include possible cutoff of magnetic braking.
      ADOT = ADOT + AGDOT
*
*       Convert from cgs to scaled units and update semi-major axis.
      ADOT = ADOT/(1.0D+05*VSTAR)
      SEMI1 = SEMI - ADOT*STEP(I1)
*
*       Include safety test on new semi-major axis.
      IF (ABS(SEMI1).LT.RADIUS(J2).OR.SEMI1.LT.0.0) THEN
          SEMI1 = RADIUS(J2)
      END IF
*
*       Transform original orbit to pericentre if R > A.
      IF (R(IPAIR).GT.SEMI) THEN
          CALL KSAPO(IPAIR)
      END IF
*
*       Form square of regularized velocity.
      V20 = 0.0
      DO 10 K = 1,4
          V20 = V20 + UDOT(K,IPAIR)**2
   10 CONTINUE
*
*       Update binding energy and collision energy.
      HI = H(IPAIR)
      H(IPAIR) = -0.5*BODY(I)/SEMI1
      ECOLL = ECOLL + ZMU*(HI - H(IPAIR))
*
*       Specify KS coordinate & velocity scaling factors at new pericentre.
      C2 = SQRT(SEMI1/SEMI)
      V2 = 0.5*(BODY(I) + H(IPAIR)*SEMI1*(1.0 - SQRT(ECC2)))
      C1 = SQRT(V2/V20)
*
*       Re-scale KS variables to new energy with constant eccentricity.
      R(IPAIR) = 0.0D0
      TDOT2(IPAIR) = 0.0D0
      DO 20 K = 1,4
          U(K,IPAIR) = C2*U(K,IPAIR)
          UDOT(K,IPAIR) = C1*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
*         TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0*U(K,IPAIR)*UDOT(K,IPAIR)
   20 CONTINUE
*
*       Transform back to apocentre for standard unperturbed motion.
      CALL KSAPO(IPAIR)
*
*       Include optional diagnostic output (also 2nd Roche stage).
      IF (KZ(28).GT.1.OR.KSTAR(I).EQ.23) THEN
          NDIAG = NDIAG + 1
          IF (NDIAG.EQ.1.OR.MOD(NDIAG,100).EQ.0) THEN
              RCOLL = RADIUS(I1) + RADIUS(I2)
              WRITE (6,25)  TIME, IPAIR, M2, R2, TBR, TGR, R(IPAIR),
     &                      RCOLL
   25         FORMAT (' BRAKE    T KS M2 R2 TBR TGR R RCOLL ',
     &                           F10.4,I4,F6.2,F7.3,1P,2E9.1,2E10.2)
          END IF
      END IF
*
*       Determine indices for primary & secondary star (donor & accretor).
*     J1 = I1
*     J2 = I2
*     Q = BODY(J1)/BODY(J2)
*     RL1 = RL(Q)*SEMI1
*       Evaluate Roche radius for the second star.
*     Q = 1.0/Q
*     RL2 = RL(Q)*SEMI1
*
*       Compare scaled Roche radii when choosing the primary.
*     IF (RADIUS(J1)/RL1.LT.RADIUS(J2)/RL2) THEN
*         J1 = I2
*         J2 = I1
*         RL1 = RL2
*     END IF
*
*       Check Roche mass transfer at every stage (AM CVn formation).
*     K1 = MAX(KSTAR(J1),KSTAR(J2))
*     IF (K1.EQ.11.OR.K1.EQ.12) THEN
*       See whether the Roche condition is satisfied (MS or largest WD).
*         IF (RADIUS(J1).GT.RL1) THEN
*             CALL ROCHE(IPAIR)
*         END IF
*         IF (IPHASE.LT.0) GO TO 40
*     END IF
*
*       Perform collision check (compact + evolved star).
      IF (R(IPAIR).LT.RADIUS(I1) + RADIUS(I2)) THEN
          RCOLL = RADIUS(I1) + RADIUS(I2)
          TK = DAYS*SEMI*SQRT(SEMI/BODY(I))
          WRITE (6,30)  TIME, NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2),
     &                  M2, R2, TBR, TGR, R(IPAIR), RCOLL, TK
   30     FORMAT (' END BRAKE    T NAM K* M2 R2 TBR TGR R RCOLL P ',
     &                           F10.4,2I6,2I4,F6.2,F7.3,1P,5E9.1)
*
*       Predict c.m. and combine the two stars inelastically.
          CALL XVPRED(I,0)
          KSPAIR = IPAIR
          CALL CMBODY(R(IPAIR),2)
      END IF
*
   40 RETURN
*
      END

