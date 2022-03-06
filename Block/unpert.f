      SUBROUTINE UNPERT(IPAIR,KCASE)
*
*       Unperturbed two-body motion.
*       ----------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'kscmn6.h'
      REAL*8  UI(4),UIDOT(4)
      INTEGER KLIST(2)
      SAVE TCALL
      DATA TCALL /0.0D0/
*
*
*       Set first component & c.m. index and semi-major axis.
      I1 = 2*IPAIR - 1
      I = N + IPAIR
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
*       Restrict spurious large values of SEMI, mainly after NEW KSREG.
      SEMI = MIN(SEMI,10.0*RMIN)
      JPHASE = 0
      KSDEF = 0
*       Check for GR shrinkage of small circularized binaries.
      DT0 = TIME - T0(I1)
*     IF (KSTAR(I).EQ.10.AND.SEMI*SU.LT.6.0) GO TO 100
*     IF (KSTAR(I).EQ.0.AND.SEMI*SU.LT.6.0) GO TO 100
      IF (SEMI*SU.LT.0.0) GO TO 100      ! Reverse for future experiments.
*
*       Add elapsed unperturbed Kepler periods and update the time T0(I1).
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
      IF (TIME - T0(I1).LT.2.0E+09*TK) THEN
          K = NINT((TIME - T0(I1))/TK)
      ELSE
          K = 0
!$omp atomic
          NPRECT = NPRECT + 1
      END IF
      T0PREV = T0(I1)
      T0(I1) = TIME
*
*       Reset unperturbed counter if > 2*10**9 and update the frequency.
      IF (NKSPER.GT.2000000000.OR.NKSPER.LT.0) THEN
          NKSPER = 0
!$omp atomic
          NPRECT = NPRECT + 1
      END IF
!$omp atomic
      NKSPER = NKSPER + K
*
*       Set a nominal time-step at apocentre.
      STEP(I1) = 10.0*TK
      IF (R(IPAIR).LT.0.1*SEMI) THEN
*       Transform unperturbed orbit from pericentre to apocentre.
          CALL KSAPO(IPAIR)
*         WRITE (6,5)  NAME(I), TIME, R(IPAIR), STEP(I1)
*   5     FORMAT (' UNPERT    NM T R S1 ',I7,F10.4,1P,E10.2)
      END IF
*
*       Include case of tidal dissipation or partial reflection (suppressed).
      IF (TDOT2(IPAIR).GE.0.0D0) THEN
*         IF (KZ(25).GT.0.AND.ABS(TDOT2(IPAIR).GT.1.0E-10) THEN
*             DT = 0.0
*             GO TO 10
*         END IF
*       Ensure perturbation check at least once every c.m. step (KSTAR < 11).
          IF (KZ(27).GT.0) THEN
              IF (KSTAR(I).LT.11) THEN
                  IF (STEP(I1).LT.TK) THEN
                      IF (STEP(I1).LT.0.0001*STEP(I)) GO TO 10
                  END IF
                  IF (TIME - T0(I).GT.2.0*STEP(I1)) THEN
                      DT = MIN(3.0D0*STEP(I1),STEP(I))
                      KPERT = 0
                      JCL = 0
                      GO TO 20
                  END IF
              END IF
          END IF
      END IF
*
*       Evaluate interval for unperturbed motion (GAMMA < GMIN).
   10 KPERT = 1
*       Determine new perturber index (thread-safe procedure).
      CALL TPERT(IPAIR,GMIN,DT,JCL)
*
*       Restore KS indicator and re-initialize if interval < period.
      IF (DT.LT.TK.OR.T0PREV.EQ.T0(I)) THEN
*       Form perturber list and restart KS motion if required.
          CALL KSLIST(IPAIR)
          IF (LIST(1,I1).GT.0) THEN
*       Transform to apocentre variables in case of tidal dissipation.
              IF (R(IPAIR).LT.SEMI.AND.SEMI.LT.5.0*RMIN) THEN
                  NP = LIST(1,I1)
                  LIST(1,I1) = 0
*       Do not allow backwards integration on switch from unperturbed motion.
                  CALL KSPERI(IPAIR)
                  LIST(1,I1) = NP
                  CALL KSAPO(IPAIR)
              END IF
              KSLOW(IPAIR) = 1
              CALL RESOLV(IPAIR,1)
              CALL KSPOLY(IPAIR,1)
*       Include differential correction (experimental 10/02).
              KLIST(1) = I
              NNB = LIST(1,I1)
              CALL NBPOT2(1,NNB,POT1,KLIST,LIST(2,I1))
              KLIST(1) = I1
              KLIST(2) = I1 + 1
              CALL NBPOT2(2,NNB,POT2,KLIST,LIST(2,I1))
*             WRITE (6,13)  NAME(I1), POT1-POT2, GAMMA(IPAIR)
*  13         FORMAT ('  CORRECT    NAM DPHI G ',I6,1P,2E10.2)
!$omp atomic
              EMERGE = EMERGE + (POT2 - POT1)
!$omp atomic
              BE(3) = BE(3) + (POT2 - POT1)
              GO TO 40
          END IF
          IR = 1
          GO TO 28
      ELSE
*       Check rectifying large eccentricity (> 1.01) for unperturbed orbit.
          ECC2 = (1.0 - R(IPAIR)/SEMI)**2 +
     &                                    TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          IF (ECC2.GT.1.02) THEN
              ECC = SQRT(ECC2)
      IF (R(IPAIR).GT.2.0*SEMI) THEN
          WRITE (6,13)  ECC, R(IPAIR), SEMI
   13     FORMAT (' MAX-WARNING    E R A ',F9.5,1P,2E10.2)
      END IF
      IF (R(IPAIR).LT.SEMI) CALL KSAPO(IPAIR)
              CALL DEFORM(IPAIR,ECC,0.99D0)
              WRITE (6,15)  IPAIR, ECC, R(IPAIR)/SEMI, H(IPAIR)
   15         FORMAT (' UNPERT DEFORM    KS E R/A H ',I6,F7.3,F9.5,F8.2)
              KSDEF = 1
              NEWKS = NEWKS + 1
          END IF
      END IF
*
*       Check for tidal dissipation (peri & apocentre included; skip merger).
   20 IF (KZ(27).GT.0.AND.
     &   (SEMI.LT.RMIN.OR.TDOT2(IPAIR).GT.0.0D0).AND.
     &   (KSTAR(I).EQ.0.OR.KSTAR(I).EQ.-1).AND.
     &   (NAME(I).GT.0)) THEN
*
*       Compare pericentre and effective capture distance.
          ECC2 = (1.0 - R(IPAIR)/SEMI)**2 +
     &                                    TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          ECC = SQRT(ECC2)
          RP = SEMI*(1.0D0 - ECC)
          RT = 4.0*MAX(RADIUS(I1),RADIUS(I1+1))
*
*       Specify circularization index (skip SPIRAL but include CHAOS).
          ICIRC = 0
          IF (KZ(27).EQ.1.AND.TTOT.GT.TCALL) THEN
              IF (RP.LT.RT) ICIRC = 1
              TCALL = TTOT + 0.01
          ELSE IF (RP.LT.2.5*RT.AND.KSTAR(I).EQ.0.AND.
     &        TTOT.GT.TCALL) THEN
*       Perform occasional checks of TC (NB! small periods).
              CALL TCIRC(RP,ECC,I1,I1+1,ICIRC,TC)
              TCALL = TTOT + 0.01
          ELSE IF (KSTAR(I).EQ.-1) THEN
              ICIRC = 1
          END IF
*
*       See whether dissipation is active (A*(1 - E) < 4*MAX(R) & E > 0.002).
*         IF (RP.LT.0.99*RT.AND.ECC.GT.0.002) THEN
          IF (ICIRC.GT.0.AND.ECC.GT.0.002) THEN
*       Transform to pericentre if R > A (KSAPO works both ways with PI/2).
              IF (R(IPAIR).GT.SEMI) THEN
                  CALL KSAPO(IPAIR)
              END IF
*       Implement energy loss at pericentre (exit on collision).
              CALL KSTIDE(IPAIR,KCASE,RP)
              IF (ISTAT(KCASE).LT.0) THEN
                  JPHASE = -1
                  GO TO 30
              END IF
*             IF (IPHASEX(IPAIR).LT.0) GO TO 30
*             SEMI = -0.5D0*BODY(I)/H(IPAIR)
*             TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
*       Assign one Kepler period for active state and impose TDOT2 > 0.
              IF (R(IPAIR).LT.0.99*RT.AND.KSTAR(I).NE.10) THEN
                  STEP(I1) = TK
                  TDOT2(IPAIR) = 1.0E-20
              ELSE
*       Define apocentre and set 1/2 period for inactive or synchronous case.
                  CALL KSAPO(IPAIR)
                  STEP(I1) = 0.5*TK
                  ECC2 = (1.0 - R(IPAIR)/SEMI)**2 +
     &                                    TDOT2(IPAIR)**2/(BODY(I)*SEMI)
                  ECC = SQRT(ECC2)
*       Note next unperturbed check at apocentre since T'' < 0 in KSAPO.
                  IF (ABS(RT - RP)/RT.GT.0.1.AND.KZ(27).LE.1) THEN
                      WRITE(6,25)  ECC, SEMI, R(IPAIR), RP, RADIUS(I1)
   25                 FORMAT (' INACTIVE PHASE    E A R RP R* ',
     &                                            F7.3,1P,4E10.2)
*       Include safety condition for extremely small SEMI (i.e. tiny period).
                      IF (SEMI.LT.0.001*RMIN) STEP(I1) = STEP(I)
                  END IF
                  IF (ECC.LT.0.002.AND.SEMI.LT.0.01*RMIN) THEN
                      KSTAR(I) = 10
                      TEV(I) = TIME
                  END IF
              END IF
              GO TO 30
          END IF
      END IF
*
*       Check for Roche overflow from synchronous orbit (but KSTAR < 13).
      IR = 0
      IF (KSTAR(I).EQ.13.AND.MIN(TEV(I1),TEV(I1+1)).LT.TIME + DT) THEN
          IR = -1
      END IF
      IF (KSTAR(I).GT.13) IR = 1
      IF (KZ(34).GT.0.AND.(KSTAR(I).GE.10.AND.KSTAR(I).LE.12)) THEN
          TM = MIN(TEV(I1),TEV(I1+1),TEV(I))
          IF (TM.LT.TIME) THEN
              CALL TRFLOW(IPAIR,DTR)
*       Exit if ROCHE is indicated (change 16/08/2006).
*             IF (DTR.LT.STEP(I1)) THEN
*                 CALL ROCHE(IPAIR)
*       Exit on coalescence (KPERT > 0 would cause trouble).
*                 IF (IPHASEX(IPAIR).LT.0) GO TO 30
*             ELSE
*                 IR = 1
*             END IF
              IR = 1
              IF (DTR.LT.STEP(I1)) GO TO 30
*
              IF (DTR.GT.TK) THEN
                  DT = MIN(DTR,DT)
                  KPERT = 2
              END IF
          ELSE
              IR = 1
          END IF
      END IF
*
*       Perform general two-body collision test at pericentre (< 3/4 R1+R2).
      IF (KZ(19).GE.3.AND.NAME(I).GT.0) THEN
          RI = 0.0
          DO 26 K = 1,4
              UI(K) = U0(K,IPAIR)
              UIDOT(K) = UDOT(K,IPAIR)
              RI = RI + UI(K)**2
   26     CONTINUE
          CALL PERI(UI,UIDOT,RI,BODY(I1),BODY(I1+1),QPERI)
          IF (QPERI.LT.0.75*(RADIUS(I1) + RADIUS(I1+1))) THEN
*       Obtain KS variables at pericentre before coalescence to one body.
*             CALL KSPERI(IPAIR)
*             KSPAIR = IPAIR
*       Set indicator for skipping ECOLL updating in COAL.
*             IQCOLL = -2
              JPHASE = -1
*             CALL CMBODY(QPERI,2)
              IF (JPHASE.LT.0) GO TO 30
          END IF
      END IF
*
*       Specify a conservative number of unperturbed orbits from TPERT.
      IF (KPERT.GT.0.AND.DT.GT.0.0) THEN
*       Ensure that next look-up time is not exceeded (option #27).
          IF (KZ(27).GT.0.AND.KPERT.GT.1) THEN
              TM = MIN(TEV(I1),TEV(I1+1))
              IF (TM - T0(I1).GT.0.0.AND.TM - T0(I1).LT.0.5*DT) THEN
!$omp atomic
                  NWARN = NWARN + 1
                  IF (NWARN.LT.1000) THEN
                      WRITE (25,27)  IPAIR, KSTAR(I1), KSTAR(I1+1),
     &                               KSTAR(I), TM-T0(I1), STEP(I1)
   27                 FORMAT (' UNPERT WARNING!    KS K* TEV-T0 STEP1 ',
     &                                             4I4,1P,2E10.2)
                      CALL FLUSH(25)
                  END IF
                  DT = MIN(2.0*(TM - T0(I1)),DT)
                  DT = MAX(DT,TK,STEPX)
              END IF
          END IF
*       Adopt c.m. step instead if integer argument exceeds 10**9.
          IF (DT.LT.2.0E+09*TK) THEN
              K = 1 + INT(0.5D0*DT/TK)
*       Restrict Kepler period to c.m. step (case of very wide orbit).
              STEP(I1) = MIN(FLOAT(K)*TK,STEP(I))
          ELSE
              STEP(I1) = STEP(I)
          END IF
*       Include optional treatment for spiralling of chaotic binary orbit.
          IF (KZ(27).GT.1.AND.KSTAR(I).EQ.-2) THEN
*       Ensure pericentre position after possible perturbed motion.
              IF (R(IPAIR).GT.SEMI) THEN
                  CALL KSRECT(IPAIR)
*       Reduce eccentric anomaly by pi for inward motion.
                  IF (TDOT2(IPAIR).LT.0.0D0) THEN
                      CALL KSAPO(IPAIR)
                  END IF
*       Transform from outward motion (anomaly < pi) to exact pericentre.
                  CALL KSPERI(IPAIR)
              END IF
              CALL SPIRAL(IPAIR)
              IF (IPHASEX(IPAIR).LT.0) JPHASE = -1
              IF (IPHASEX(IPAIR).LT.0) GO TO 30
          END IF
      END IF
*
*       Check merger condition before continuing (skip small Roche steps).
   28 IF (KZ(15).GT.0.AND.IR.GE.0.AND.TIME+STEP(I1).GT.TBLOCK) THEN
          IF (STEP(I).LT.DTMIN) THEN
              CALL IMPACT(I,JPHASE)
          ELSE IF (JCL.GT.0.AND.STEP(I).LT.10.0*DTMIN) THEN
              CALL HISTAB(IPAIR,JCL,PMIN,RSTAB)
              IF (RSTAB.LT.PMIN) THEN
                  CALL IMPACT(I,JPHASE)
              END IF
*       Include a more generous test for massive quadruples.
          ELSE IF (JCL.GT.N) THEN
              FAC = 2.0*(BODY(I) + BODY(JCL))/BODYM
              IF (STEP(I).LT.FAC*DTMIN) THEN
                  CALL HISTAB(IPAIR,JCL,PMIN,RSTAB)
                  IF (RSTAB.LT.PMIN) THEN
                      CALL IMPACT(I,JPHASE)
                  END IF
              END IF
          END IF
      END IF
*
*       Check collision criterion for special case (unusual option).
      IF (KZ(27).EQ.-1.AND.KZ(13).LT.0) THEN
          ECC2 = (1.0 - R(IPAIR)/SEMI)**2 +
     &                                TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          ECC = SQRT(ECC2)
          QPERI = SEMI*(1.0 - ECC)
          RFAC = 2.0
          I2 = I1 + 1
          IF (QPERI.LT.RFAC*MAX(RADIUS(I1),RADIUS(I2))) THEN
              J1 = I1
              IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
              FAC = 0.5*BODY(I)/BODY(J1)
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
              RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
              IF (QPERI.LT.RCOLL) THEN
                  CALL TOUCH(IPAIR,I1,I2,RCOLL)
              END IF
          END IF
      END IF
*
*       Set quantized time-step from DT > TK on continued unperturbed motion.
   30 DT = MAX(DT,TK)
      DT = MIN(DT, 8.0*STEP(I), STEPR(I))
      CALL STEPK(DT,DTN)
      DO WHILE (DMOD(TIME,DTN).NE.0.0D0)
          DTN = 0.5 * DTN
      END DO
      STEP(I1) = DTN
*
   40 ISTAT(KCASE) = JPHASE
*
      IF (KSDEF.GT.0) THEN
          ECC2 = (1.0 - R(IPAIR)/SEMI)**2 +
     &                                    TDOT2(IPAIR)**2/(BODY(I)*SEMI)
              WRITE (6,55)  IPAIR, ECC, R(IPAIR)/SEMI, H(IPAIR)
   55         FORMAT (' UNPERT MODIFY    KS E R/A H ',I6,F7.3,F9.5,F8.2)
*     END IF
          IF (R(IPAIR).GT.2.0*SEMI) THEN
              WRITE (6,60)  IPAIR, ECC, R(IPAIR), SEMI
   60         FORMAT (' UNP-WARNING    KS E R A ',I5,F9.5,1P,2E10.2)
              CALL FLUSH(6)
          END IF
      END IF
      GO TO 200
* 
*       Include GR shrinkage to eliminate small circularized orbits.
  100 CONTINUE
*     IF ((KSTAR(I).EQ.10.OR.KSTAR(I).EQ.0).AND.SEMI*SU.LT.6.0) THEN
      IF (SEMI*SU.LT.6.0) THEN
          KSPAIR = IPAIR
          CVEL = 3.0D+05/VSTAR
          ZFAC = BODY(I1)*BODY(I1+1)*BODY(I)*SMU**3
*         ZFAC = BODY(I1)*BODY(I1+1)*BODY(I)
          ADOT =  64.0/5.0*ZFAC/(CVEL**5*SEMI**3)
*       Shrink circular orbit according to Peters 1964.
          SEMI0 = SEMI - ADOT*DT0
      DA = ADOT*DT0/SEMI0
      WRITE (6,104) KSTAR(I), SEMI, SEMI0, ADOT*DT0, DA, ADOT
  104 FORMAT (' GR SHRINK    K* A A0 AD*DT DA/A ',I4,1P,5E10.2)
          IF (SEMI0.LT.0.9*SEMI) SEMI0 = 0.9*SEMI
*         WRITE (6,105)  NAME(I1), TIME, SEMI, SEMI0, ADOT
* 105     FORMAT (' EXPAND   NM1 T A A0 AD ',I6,F11.7,1P,5E10.2)
*       Employ contraction for safety scaling of U & UDOT.
          CALL EXPAND2(IPAIR,SEMI0)
*
*       Assume circular orbit for energy conservation via ECOLL.
          HI = H(IPAIR)
          H(IPAIR) = -0.5*BODY(I)/SEMI0
          ZMU = BODY(I1)*BODY(I1+1)/BODY(I)
          ECOLL = ECOLL + ZMU*(HI - H(IPAIR))
*
*       Adopt enforced coalescence condition.
          R12 = RADIUS(I1) + RADIUS(I1+1)
          R12 = 0.9*R12
          R(IPAIR) = SEMI0
*         WRITE (6,110)  TIME, NAME(I1), SEMI0, R(IPAIR)
* 110     FORMAT (' CONTRACT    T NM A0 R ',F10.6,I6,1P,2E10.2)
          IF (SEMI0.LT.R12) THEN
  120         FORMAT (' GR COAL    NAM R12 A ',I6,1P,2E10.2)
              CALL CMBODY(R12,2)
              IPHASE = -1
          END IF
      END IF
*
  200 RETURN
*
      END
