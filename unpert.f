      SUBROUTINE UNPERT(IPAIR)
*
*
*       Unperturbed two-body motion.
*       ----------------------------
*
      INCLUDE 'common6.h'
*
*
*       Set first component & c.m. index and semi-major axis.
      I1 = 2*IPAIR - 1
      I = N + IPAIR
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
*
*       Add elapsed unperturbed Kepler periods and update the time T0(I1).
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
      K = NINT((TIME - T0(I1))/TK)
      T0(I1) = TIME
*
*       Reset unperturbed counter if > 2*10**9 and update the frequency.
      IF (NKSPER.GT.2000000000.OR.NKSPER.LT.0) THEN
          NKSPER = 0
          NPRECT = NPRECT + 1
      END IF
      NKSPER = NKSPER + K
*
*       Include case of tidal dissipation or partial reflection (suppressed).
      IF (TDOT2(IPAIR).GE.0.0D0) THEN
*         IF (KZ(25).GT.0.AND.ABS(TDOT2(IPAIR).GT.1.0E-10) THEN
*             DT = 0.0
*             GO TO 10
*         END IF
*       Ensure perturbation check at least once every c.m. step.
          IF (KZ(27).GT.0) THEN
              DTR = (TIME - T0(I)) - 0.5*STEP(I)
              IF (ABS(DTR).GT.STEP(I1)) THEN
                  KPERT = 0
                  GO TO 20
              END IF
          END IF
      END IF
*
*       Evaluate interval for unperturbed motion (GAMMA < GMIN).
      KPERT = 1
      CALL TPERT(IPAIR,GMIN,DT)
*
*       Restore KS indicator and re-initialize if interval < period.
   10 IF (DT.LT.TK) THEN
*       Form perturber list and restart KS motion if required.
          CALL KSLIST(IPAIR)
          IF (LIST(1,I1).GT.0) THEN
*       Transform to apocentre variables in case of tidal dissipation.
              IF (KZ(27).GT.0.AND.R(IPAIR).LT.SEMI) THEN
                  CALL KSAPO(IPAIR)
              END IF
              KSLOW(IPAIR) = 1
              CALL RESOLV(IPAIR,1)
              CALL KSPOLY(IPAIR,1)
          END IF
          GO TO 30
      END IF
*
*       Consider case of tidal capture (both peri and apo included).
   20 IF (KZ(27).GT.0.AND.
     &   (SEMI.LT.0.01*RMIN.OR.TDOT2(IPAIR).GT.0.0D0).AND.
     &    KSTAR(I).NE.20) THEN
*
*       Compare pericentre and effective capture distance.
          ECC2 = (1.0 - R(IPAIR)/SEMI)**2 +
     &                                    TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          ECC = SQRT(ECC2)
          RP = SEMI*(1.0D0 - ECC)
          RT = 4.0*MAX(RADIUS(I1),RADIUS(I1+1))
*
*       See whether dissipation is active (A*(1 - E) < 4*MAX(R) & E > 0.015).
          IF (RP.LT.0.99*RT.AND.ECC.GT.0.015) THEN
*       Transform to pericentre if R > A (KSAPO works both ways with PI/2).
              IF (R(IPAIR).GT.SEMI) THEN
                  CALL KSAPO(IPAIR)
              END IF
*       Implement energy loss at pericentre.
              CALL KSTIDE(IPAIR,RP)
              SEMI = -0.5D0*BODY(I)/H(IPAIR)
              TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
*       Assign one Kepler period for active state and impose TDOT2 > 0.
              IF (R(IPAIR).LT.0.99*RT.AND.KSTAR(I).NE.20) THEN
                  STEP(I1) = TK
                  TDOT2(IPAIR) = 1.0E-20
              ELSE
*       Define apocentre and set 1/2 period for inactive or synchronous case.
                  CALL KSAPO(IPAIR)
                  STEP(I1) = 0.5*TK
*       Note next unperturbed check at apocentre since T'' < 0 in KSAPO.
                  WRITE(6,25)  ECC, SEMI, R(IPAIR), RP
   25             FORMAT (' INACTIVE PHASE    E A R RP ',F7.3,1P,3E9.1)
                  IF (ECC.LT.0.002.AND.SEMI.LT.0.01*RMIN) THEN
                      KSTAR(I) = 20
                  END IF
              END IF
              GO TO 30
          END IF
      END IF
*
*       Specify a conservative number of unperturbed orbits from TPERT.
      IF (KPERT.GT.0) THEN
*       Adopt c.m. step instead if integer argument exceeds 10**9.
          IF (DT.LT.2.0D+09*TK) THEN
              K = 1 + INT(0.5D0*DT/TK)
*       Restrict Kepler period to c.m. step (case of wide orbit).
              STEP(I1) = FLOAT(K)*MIN(TK,DT)
          ELSE
              STEP(I1) = STEP(I)
          END IF
      END IF
*
*       Check for optional magnetic braking or gravitational radiation.
      IF (KZ(28).GT.0.AND.SEMI*SU.LT.10.0.AND.KSTAR(I).GT.0) THEN
          CALL BRAKE(IPAIR)
          IF (IPHASE.LT.0) GO TO 30
      END IF
*
*       Check merger condition before continuing unperturbed motion.
      IF (KZ(15).GT.0.AND.STEP(I).LT.DTMIN) THEN
*       Skip possible second call during termination.
          IF (IPHASE.EQ.0) THEN
              CALL IMPACT(I)
          END IF
      END IF
*
   30 RETURN
*
      END

