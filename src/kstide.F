      SUBROUTINE KSTIDE(IPAIR,QPERI)
!VAMPINST include
*
*
*       Tidal interaction of KS pair.
*       -----------------------------
*
      INCLUDE 'common6.h'
      REAL*8  DE(2)
      INTEGER IS(2)
      DATA  ECCM,ECCM2  /0.002,0.00000399/
*
*
*       Skip procedure if |QP - 4*R| < 0.01*QP or both stars zero size.
!VAMPINST subroutine_start
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      RX = MAX(RADIUS(I1),RADIUS(I2))
      SEMI = -0.5D0*BODY(N+IPAIR)/H(IPAIR)
      ECC = 1.0 - QPERI/SEMI
      ZF = 4.0
      IF (ECC.GT.0.95) ZF = 50.0
      IF (ABS(QPERI - ZF*RX).LT.0.01*QPERI) GO TO 50
      IF (RADIUS(I1) + RADIUS(I2).LE.0.0D0) GO TO 50
*
      IF (QPERI.LT.RADIUS(I1) + RADIUS(I2)) THEN
          if(rank.eq.0)
     &    WRITE (6,1)  NAME(I1), QPERI, RX
    1     FORMAT (' COLLISION TEST!    NAM QP R*  ',I5,1P,2E10.2)
      END IF
      IS(1) = KSTAR(I1)
      IS(2) = KSTAR(I2)
*
*       Obtain kinetic energy loss due to tidal interaction (DE > 0 here).
      CALL TIDES(QPERI,BODY(I1),BODY(I2),RADIUS(I1),RADIUS(I2),IS,DE)
*
*       Set c.m. index & reduced mass.
      I = N + IPAIR
      ZMU = BODY(I1)*BODY(I2)/BODY(I)
*
*       Determine pericentre variables U & UDOT by backwards reflection.
      CALL KSPERI(IPAIR)
*
*       Form semi-major axis & eccentricity (TDOT2 = 0 at pericentre).
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      ECC = 1.0 - R(IPAIR)/SEMI
      PERI = SEMI*(1.0D0 - ECC)
*
*       Restore circularization index if needed (exit from CHAIN).
      IF (ECC.LE.ECCM.AND.KSTAR(I).LT.20) THEN
          KSTAR(I) = 20
          GO TO 50
      END IF
*
*       Include safety check on energy loss to prevent new SEMI < R.
      DH = -(DE(1) + DE(2))/ZMU
      IF (H(IPAIR) + DH.LT.-0.5*BODY(I)/R(IPAIR)) THEN
          DH = -0.5*BODY(I)/R(IPAIR) - H(IPAIR)
          DE(1) = -ZMU*DH
          DE(2) = 0.0
      END IF
*
*       Obtain the new eccentricity from angular momentum conservation.
      AM0 = SEMI*(1.0D0 - ECC**2)
*       Suppress the old instantaneous PT procedure.
*     ECC2 = ECC**2 + 2.0D0*AM0*DH/BODY(I)
*     ECC2 = MAX(ECC2,ECCM2)
*
*       Adopt sequential circularization instead of standard PT.
      ECC2 = ECCM2
      ACIRC = AM0/(1.0 - ECC2)
      SEMI1 = AM0/(1.0 - ECC2)
      ECC1 = SQRT(ECC2)
*
*       Check Roche conditions for decreased semi-major axis.
      Q1 = BODY(I1)/BODY(I2)
      RL1 = 0.49*SEMI1/(0.6 + LOG(1.0 + Q1**0.333)/Q1**0.667)
      Q2 = 1.0/Q1
      RL2 = 0.49*SEMI1/(0.6 + LOG(1.0 + Q2**0.333)/Q2**0.667)
*
*       Find the dominant Roche radius (not used yet).
      IF (RADIUS(I1)/RL1.GE.RADIUS(I2)/RL2) THEN
          RL = RADIUS(I1)
      ELSE
          RL = RADIUS(I2)
          RL1 = RL2
      END IF
*
*       Accept circularized orbit if ACIRC < 4*R1 (use maximum radius).
      R1 = RX
      IF (ACIRC.LT.ZF*R1) THEN
          SEMI1 = ACIRC
      ELSE
*       Obtain E1 from A1*(1 - E1**2) = AM0 using A1*(1 - E1) = 4*R1.
          ECC1 = AM0/(ZF*R1) - 1.0
          ECC1 = MAX(ECC1,ECCM)
          ECC1 = MAX(ECC1,0.9*ECC)
*       Set new semi-major axis from angular momentum conservation.
          SEMI1 = AM0/(1.0 - ECC1**2)
      END IF
*
*       Form the corresponding energy change.
      DH = 0.5*BODY(I)*(1.0/SEMI - 1.0/SEMI1)
*     DH = (ECC2 - ECC**2)*BODY(I)/(2.0D0*AM0)
      DE(1) = -ZMU*DH
      DE(2) = 0.0
*
*       Skip possible hyperbolic case.
      IF (H(IPAIR) + DH.GT.0.0) GO TO 50
*
*       Increase event counter and update total energy loss.
      NDISS = NDISS + 1
      ECOLL = ECOLL + (DE(1) + DE(2))
      E(10) = E(10) + (DE(1) + DE(2))
*
*       Set new energy and pericentre.
      HI = H(IPAIR)
      H(IPAIR) = H(IPAIR) + DH
      PERI1 = SEMI1*(1.0D0 - ECC1)
*
*       Print first energy change and activate indicator.
      IF (KSTAR(I).EQ.0) THEN
          P = DAYS*SEMI1*SQRT(SEMI1/BODY(I))
          if(rank.eq.0)
     &    WRITE (6,5)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2),
     &                 TIME+TOFF, ECC, ECC1, P, SEMI1, RX
    5     FORMAT (' NEW CIRC    NAM K* T E0 EF P AF R* ',
     &                          2I6,2I4,F9.2,2F8.3,F7.1,1P,2E10.2)
          KSTAR(I) = 19
      END IF
*
*       Form KS coordinate scaling factor from pericentre ratio.
      C1 = SQRT(PERI1/PERI)
*       Specify KS velocity scaling from angular momentum conservation.
      C2 = 1.0/C1
*
*       See whether circular orbit condition applies.
*     AM = SEMI1*(1.0D0 - ECC1**2)
*     IF (ECC1.LE.ECCM) C2 = SQRT(AM/AM0)/C1
*
*       Transform KS variables to yield the prescribed elements.
      R(IPAIR) = 0.0D0
      DO 10 K = 1,4
          U(K,IPAIR) = C1*U(K,IPAIR)
          UDOT(K,IPAIR) = C2*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
   10 CONTINUE
*
*       Rectify the orbital elements.
      CALL KSRECT(IPAIR)
*
*       Form new perturber list after significant energy loss.
      NP0 = LIST(1,I1)
      IF (ABS(SEMI1/SEMI).LT.0.5) THEN
          CALL KSLIST(IPAIR)
      END IF
*
*       Re-initialize KS polynomials at pericentre for perturbed case.
      T0(I1) = TIME
      IF (NP0.GT.0) THEN
          CALL RESOLV(IPAIR,1)
          IMOD = KSLOW(IPAIR)
          CALL KSPOLY(IPAIR,IMOD)
      END IF
*
*       Check for hierarchical configuration with eccentric inner binary.
      IF (ECC.GT.0.99.AND.HI.LT.0.0) THEN
          NP1 = LIST(1,I1) + 1
          DO 30 L = 2,NP1
              J = LIST(L,I1)
              RIJ2 = 0.0
              VIJ2 = 0.0
              RDOT = 0.0
              DO 25 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
                  VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
                  RDOT = (X(K,I) - X(K,J))*(XDOT(K,I) - XDOT(K,J))
   25         CONTINUE
              RIP = SQRT(RIJ2)
              A1 = 2.0/RIP - VIJ2/(BODY(I) + BODY(J))
              A1 = 1.0/A1
              IF (1.0/A1.GT.0.5/RMIN) THEN
                  ECC2 = (1.0 - RIP/A1)**2 +
     &                                  RDOT**2/(A1*(BODY(I) + BODY(J)))
                  RP = A1*(1.0 - SQRT(ECC2))
                  RA = SEMI*(1.0 + ECC)
                  SR = RP/RA
                  if(rank.eq.0)
     &            WRITE (6,28)  IPAIR, H(IPAIR), SEMI, A1, RP,
     &                          SQRT(ECC2), SR
   28             FORMAT (' HIERARCHY:    IPAIR H A0 A1 RP E1 SR ',
     &                                    I4,F7.0,1P,3E9.1,0P,F6.2,F6.1)
              END IF
   30     CONTINUE
      END IF
*
*       Specify one unperturbed period for small apocentre perturbation.
      GA = GAMMA(IPAIR)*(SEMI1*(1.0 + ECC1)/R(IPAIR))**3
      IF (GA.LT.GMIN.AND.SEMI1.GT.0.0) THEN
          STEP(I1) = TWOPI*SEMI1*SQRT(SEMI1/BODY(I))
          LIST(1,I1) = 0
      END IF
*
*       Ensure T'' = 0 for pericentre test in KSINT and disspation in UNPERT.
      IF (TDOT2(IPAIR).LT.0.0D0) THEN
          TDOT2(IPAIR) = 0.0D0
      END IF
*
*       Count any hyperbolic captures.
      IF (SEMI.LT.0.0.AND.SEMI1.GT.0.0) THEN
          NTIDE = NTIDE + 1
          QPS = QPERI/MAX(RADIUS(I1),RADIUS(I2))
          if(rank.eq.0)
     &    WRITE (6,35)  TIME+TOFF, NAME(I1), NAME(I2), ECC, ECC1, QPS
   35     FORMAT (' NEW CAPTURE    T NM E EF QP/R*  ',
     &                             F9.2,2I6,F9.4,2F6.2)
      END IF
*
*       Record diagnostics for new synchronous orbit and activate indicator.
      IF (ECC.GT.ECCM.AND.ECC1.LE.ECCM.AND.KSTAR(I).NE.20) THEN
          NSYNC = NSYNC + 1
          ESYNC = ESYNC + ZMU*H(IPAIR)
          KSTAR(I) = 20
*       Set look-up time and use KZ(34) to prevent large value in MDOT.
          TEV(I) = TIME
          KZ(34) = 1
          SEMI2 = -0.5*BODY(I)/H(IPAIR)
          P = DAYS*SEMI2*SQRT(SEMI2/BODY(I))
          if(rank.eq.0)
     &    WRITE (6,40)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2),
     &                  TIME+TOFF, ECC, ECC1, P, SEMI1, RX
   40     FORMAT (' END CIRC    NAM K* T E0 EF P AF R* ',
     &                          2I6,2I4,F9.2,2F8.3,F7.1,1P,2E10.2)
      ELSE
*       See whether a low-eccentricity synchronous state has been reached.
          RCOLL = 0.75*(RADIUS(I1) + RADIUS(I2))
          IF (ABS(SEMI1).LT.1.5*RCOLL.AND.ECC.LT.ECCM) THEN
              KSTAR(I) = 20
              if(rank.eq.0)
     &        WRITE (6,45)  ECC1, SEMI1, R(IPAIR), RCOLL
   45         FORMAT (' INACTIVE PHASE    E A R RCOLL ',F7.3,1P,3E9.1)
          END IF
      END IF
*
*       Include warning if new eccentricity exceeds old value.
      ECC2 = 1.0 - R(IPAIR)/SEMI1
      IF (ECC2.GT.MAX(ECC,ECCM)) THEN
          if(rank.eq.0)
     &    WRITE (12,48)  TIME+TOFF, IPAIR, ECC2, ECC, R(IPAIR), SEMI1
   48     FORMAT (' WARNING!    E > E0    T IP E E0 R A ',
     &                                    F9.3,I4,2F8.4,1P,2E9.1)
          CALL FLUSH(12)
      END IF
*
!VAMPINST return
50        RETURN
      RETURN
*
!VAMPINST subroutine_end
      END
