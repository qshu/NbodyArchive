      SUBROUTINE IMPACT(I)
*
*
*       Multiple collision or merger search.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX)
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      REAL*8  XX(3,3),VV(3,3)
      CHARACTER*8  WHICH1
      SAVE  NMARG,NWARN
      DATA  NMARG,NWARN /0,0/
*
*
*       Set index of KS pair & both components of c.m. body #I.
      IPAIR = I - N
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      NTTRY = NTTRY + 1
      PERT1 = 0.0
      PERT2 = 0.0
      JCOMP = I
      NP = 0
      KS2 = 0
      RMAX2 = 1.0
*
*       Search c.m. neighbours if binary has at most two perturbers.
      IF (LIST(1,I1).LE.2) I1 = I
      NNB2 = LIST(1,I1) + 1
*
*       Find the dominant body (JCOMP) and nearest perturber (JMAX).
      DO 10 L = 2,NNB2
          J = LIST(L,I1)
          NP = NP + 1
          JLIST(NP) = J
          RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                  (X(3,I) - X(3,J))**2
          PERT = BODY(J)/(RIJ2*SQRT(RIJ2))
          IF (PERT.GT.PERT2) THEN 
              IF (PERT.GT.PERT1) THEN
                  RJMIN2 = RIJ2
                  JMAX = JCOMP
                  JCOMP = J
                  PERT2 = PERT1
                  PERT1 = PERT
              ELSE
                  JMAX = J
                  PERT2 = PERT
                  RMAX2 = RIJ2
              END IF
          END IF
   10 CONTINUE
*
      RDOT = (X(1,I) - X(1,JCOMP))*(XDOT(1,I) - XDOT(1,JCOMP)) +
     &       (X(2,I) - X(2,JCOMP))*(XDOT(2,I) - XDOT(2,JCOMP)) +
     &       (X(3,I) - X(3,JCOMP))*(XDOT(3,I) - XDOT(3,JCOMP))
*
*       Specify larger perturbation for optional chain regularization.
      IF (KZ(30).GT.0.AND.NCH.EQ.0) THEN
          GSTAR = 100.0*GMIN
          KCHAIN = 1
      ELSE
          GSTAR = GMIN
          KCHAIN = 0
      END IF
*
*       Only accept inward motion or small secondary perturbation.
      PERT3 = 2.0*R(IPAIR)**3*PERT2/BODY(I)
      IF (RDOT.GT.0.0.OR.PERT3.GT.100.0*GSTAR) GO TO 100
*
*       Include impact parameter test to distinguish different cases.
      A2 = (XDOT(1,I) - XDOT(1,JCOMP))**2 + 
     &     (XDOT(2,I) - XDOT(2,JCOMP))**2 +
     &     (XDOT(3,I) - XDOT(3,JCOMP))**2
      RIJ = SQRT(RJMIN2)
      A3 = 2.0/RIJ - A2/(BODY(I) + BODY(JCOMP))
      SEMI1 = 1.0/A3
      A4 = RDOT**2/(SEMI1*(BODY(I) + BODY(JCOMP)))
      ECC1 = SQRT((1.0D0 - RIJ/SEMI1)**2 + A4)
      PMIN = SEMI1*(1.0D0 - ECC1)
*
*       Set semi-major axis, eccentricity & apocentre of inner binary.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      A0 = SEMI
      ECC2 = (1.0D0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
      APO = ABS(SEMI)*(1.0 + ECC)
*
*       Form binding energy of inner & outer binary.
      I1 = 2*IPAIR - 1
      EB = BODY(I1)*BODY(I2)*H(IPAIR)/BODY(I)
      EB1 = -0.5*BODY(JCOMP)*BODY(I)/SEMI1
*
*       Obtain the total perturbing force acting on body #I & JCOMP.
      CALL FPERT(I,JCOMP,NP,PERT)
*
*       Choose maximum of dominant scalar & total vectorial perturbation.
      PERT = PERT*RJMIN2/(BODY(I) + BODY(JCOMP))
      PERT4 = 2.0*RJMIN2*RIJ*PERT2/(BODY(I) + BODY(JCOMP))
      PERTM = MAX(PERT4,PERT)
*
*       Use combined semi-major axis for binary-binary collision.
      IF (JCOMP.GT.N) THEN
          JPAIR = JCOMP - N
          SEMI2 = -0.5D0*BODY(JCOMP)/H(JPAIR)
          J1 = 2*JPAIR - 1
          EB2 = -0.5*BODY(J1)*BODY(J1+1)/SEMI2
*       Define SEMI0 as smallest binary in case IPAIR denotes widest pair.
          SEMI0 = MIN(ABS(SEMI),ABS(SEMI2))
          SEMI = SEMI + SEMI2
          APO = APO + MAX(ABS(SEMI2),R(JPAIR))
*       Do not allow negative or soft cross section.
          IF (1.0/SEMI.LT.1.0/RMIN) GO TO 100
*       Retain KS treatment for PMIN > SEMI and large semi-major axis ratio.
          IF (PMIN.GT.SEMI.AND.SEMI2.GT.20.0*SEMI0) GO TO 30
      END IF
*
*       Check separation in case of chain regularization.
      IF (KCHAIN.GT.0) THEN
*       Form effective gravitational radius (combine triple & quad).
          EBT = EB + EB1
          ZMM = BODY(I1)*BODY(I2) + BODY(I)*BODY(JCOMP)
*       Set length of chain for decision-making (also used at termination).
          RSUM = R(IPAIR) + RIJ
          RI = R(IPAIR)
          IF (JCOMP.GT.N) THEN
              EBT = EBT + EB2
              ZMM = ZMM + BODY(J1)*BODY(J1+1)
              RSUM = RSUM + R(JPAIR)
              RI = MAX(R(JPAIR),RI)
          END IF
          RGRAV = ZMM/ABS(EBT)
          RGRAV = MIN(RGRAV,RMIN)
*       Save initial energy in binaries for routine SETSYS.
          EBCH0 = EBT - EB1
*       Use RIJ instead of RSUM in 3*RGRAV test (increases initial RIJ).
          IF (RIJ.GT.MAX(3.0*RGRAV,RMIN).OR.RGRAV.GT.RMIN) GO TO 30
*       Allow perturber distance 40*RGRAV because of the short duration.
          IF (40.0*RGRAV.GT.RS(I)) GO TO 30
          GI = 2.0*BODY(JCOMP)*(RI/RIJ)**3/BODY(I)
*       Enforce KS orbit using MERGE for high eccentricity if PMIN > 10*RI.
          IF (ECC1.GT.0.99.AND.PMIN.GT.10.0*RI.AND.
     &        PERTM.LT.GMAX) GO TO 40
*         IF (GI.LT.0.001) GO TO 30
          IF (KZ(27).GT.0.AND.JCOMP.GT.N) THEN
              IF (SEMI0.LT.SEMI2) THEN
                  J1 = 2*(I - N) - 1
              END IF
              RT = 4.0*MAX(RADIUS(J1),RADIUS(J1+1))
*       Do not allow large distance ratio for nearly synchronous binary.
              IF (SEMI0.GT.RT.AND.RI.GT.25.0*SEMI0) GO TO 30
          END IF
      END IF
*
*       Include special case of strong interaction and large ECC1.
      IF (ECC1.GT.0.9.AND.GAMMA(IPAIR).GT.0.01) THEN
          IF (APO.LT.0.01*RMIN.AND.PMIN.LT.2.5*APO) GO TO 17
      END IF
*
*       Adopt triple, quad or chain regularization for strong interactions.
      IF ((APO.GT.0.01*RMIN.OR.JCOMP.GT.N).AND.PMIN.GT.1.5*APO) GO TO 30
      IF ((RIJ.GT.RMIN.AND.SEMI1.GT.0.0).OR.RIJ.GT.2.0*RMIN) GO TO 100
      IF (PERTM.GT.100.0*GSTAR) GO TO 30
   17 IF (JCOMP.GT.N.AND.PMIN.GT.0.1*RMIN) THEN
          IF (PMIN.GT.A0 + SEMI2) GO TO 30
      END IF
      IF (R(IPAIR).GT.RMIN) GO TO 30
*
*       Check almost stable triples (factor 1.2 is experimental).
      IF (JCOMP.LE.N.AND.PMIN.GT.2.5*SEMI) THEN
          CALL HISTAB(IPAIR,JCOMP,PMIN,RSTAB)
          RA = SEMI1*(1.0 + ECC1)
          IF (SEMI1.LT.0.0) RA = RIJ
          GI = PERT*(RA/RIJ)**3
*       Use estimated apocentre perturbation for decision-making.
          IF (PMIN.GT.1.2*RSTAB) THEN
              IF (GI.LT.0.05) GO TO 30
*       Choose chain for critical case of highly eccentric outer orbit.
              IF (ECC1.LT.0.95) GO TO 100
          ELSE IF (PMIN.GT.RSTAB) THEN
*       Treat marginally stable triple according to external perturbation.
              IF (GI.LT.0.05) GO TO 30
              IF (GI.LT.1.0.OR.ECC1.LT.0.9) GO TO 100
          END IF
      END IF
*
*       Specify maximum size of unperturbed motion.
      IF (PERT2.GT.0.0) THEN
         RPERT = (100.0*GSTAR*(BODY(I) + BODY(JCOMP))/(2.0*PERT2))**0.33
      ELSE
         RPERT = 10.0*SEMI
      END IF
*
*       Skip chain if merged binary or chain c.m. (denoted by NAME <= 0).
      IF (NAME(I).LE.0.OR.NAME(JCOMP).LE.0) GO TO 100
*
*       Compare with existing subsystem of same type (if any).
      IF (NSUB.GT.0.AND.KCHAIN.EQ.0) THEN
          IGO = 0
          CALL PERMIT(RPERT,IGO)
          IF (IGO.GT.0) THEN
              NWARN = NWARN + 1
              IF (NWARN.LT.50) WRITE (6,19)  RPERT
   19         FORMAT (' IMPACT    TERMINATION REQUEST    RPERT',1P,E9.1)
              GO TO 100
          END IF
      END IF
*
      WHICH1 = ' TRIPLE '
      IF (JCOMP.GT.N) WHICH1 = ' QUAD   '
      IF (KCHAIN.GT.0) WHICH1 = ' CHAIN '
      IF (H(IPAIR).GT.0.0) THEN
          WRITE (6,18)  I, JCOMP, ECC, ECC1, SEMI1, RIJ, GAMMA(IPAIR)
   18     FORMAT (' HYP CHAIN    I J E E1 A1 RIJ G  ',
     &                           2I6,2F7.3,1P,3E9.1)
      END IF
*
      IF (KZ(15).GT.1.OR.KZ(30).GT.1) THEN
          WRITE (6,20)  WHICH1, IPAIR, TIME+TOFF, H(IPAIR), R(IPAIR),
     &                  BODY(I), BODY(JCOMP), PERT4, RIJ, PMIN,
     &                  EB1/EB, LIST(1,I1)
   20     FORMAT (/,' NEW',A8,I5,'  T =',F8.2,'  H =',F7.0,
     &              '  R =',1P,E8.1,'  M =',0P,2F7.4,'  G4 =',1P,E8.1,
     &              '  R1 =',E8.1,'  P =',E8.1,'  E1 =',0P,F6.3,
     &              '  NP =',I2)
          CALL FLUSH(6)
      END IF
*
*       Include any close single perturber (cf. routine SETSYS).
      IF (JMAX.NE.JCOMP.AND.SQRT(RMAX2).LT.2.0*RIJ.AND.
     &    JMAX.LE.N.AND.NAME(JMAX).GT.0) THEN
          WRITE (6,22)  NAME(JCOMP), NAME(JMAX), RSUM, SQRT(RMAX2)
   22     FORMAT (' B+2 CHAIN    NAM RSUM RMX ',2I6,1P,2E10.2)
          CALL XVPRED(JMAX,-1)
          JCMAX = JMAX
      ELSE
          JCMAX = 0
      END IF
*
*       Predict coordinates and velocities of #JCOMP & c.m. to F3DOT.
      CALL XVPRED(JCOMP,-1)
      CALL XVPRED(I,0)
*
*       Save global index of intruder for TRIPLE or CHAIN.
      JCLOSE = JCOMP
*
*       Replace unperturbed near-synchronous binary by inert body in CHAIN.
      IF (KCHAIN.GT.0.AND.JCOMP.GT.N) THEN
          K1 = 2*JPAIR - 1
          WRITE (6,24)  NAME(I1), NAME(I2), NAME(K1), NAME(K1+1),
     &                  KSTAR(I), KSTAR(JCOMP), ECC, ECC1, A0, SEMI2,
     &                  RIJ, SEMI1, PMIN
   24     FORMAT (' CHAIN B-B    NAME K* E0 E1 A0 A2 RIJ A1 PM ',
     &                           4I6,2I4,2F7.3,1P,5E10.2)
          IF (SEMI0.LT.4.0*RT.AND.LIST(1,J1).EQ.0) THEN
              IF (SEMI0.LT.SEMI2) THEN
                  KPAIR = JPAIR
                  JPAIR = IPAIR
                  IPAIR = KPAIR
                  JCLOSE = N + JPAIR
              END IF
*       Check reduction of c.m. index (JPAIR becomes JPAIR - 1 if > IPAIR).
              IF (JPAIR.GT.IPAIR) JCLOSE = JCLOSE - 1
              IF (KZ(26).LT.2) THEN
*       Replace unperturbed near-synchronous binary by inert body in CHAIN.
                  JCOMP = 0
                  WRITE (6,25)  SEMI0, RIJ, R(JPAIR), GAMMA(JPAIR)
   25             FORMAT (' INERT BINARY    A0 RIJ R G ',1P,4E10.2)
              END IF
          ELSE
              JCLOSE = 0
          END IF
      END IF
*
*       Set phase indicator for calling TRIPLE or QUAD from MAIN.
      IPHASE = 4
      KSPAIR = IPAIR
*
*       Include the case of two interacting KS pairs.
      IF (JCOMP.GT.N) THEN
          IPHASE = 5
*       Switch pair indices and rename JCOMP if JPAIR has smaller step.
          IF (STEP(J1).LT.STEP(I1)) THEN
              KSPAIR = JPAIR
              JCOMP = I
              KS2 = IPAIR
          ELSE
              KS2 = JPAIR
          END IF
*       Terminate first KS pair (replaced by delay procedure).
*         CALL KSTERM
*       Reduce second index if higher than first.
          IF (KS2.GT.KSPAIR) KS2 = KS2 - 1
      END IF
*
*       See whether chain regularization indicator should be activated.
      IF (KCHAIN.GT.0) THEN
          IPHASE = 8
      END IF
*
*       Save KS indices and delay initialization until end of block step.
      CALL DELAY(KCHAIN,KS2)
*
*       Terminate binary in triple or widest B-B binary (replaced by delay).
*     CALL KSTERM
*
      IF (N.GT.10) GO TO 999
*       Prepare procedure for chain between hierarchy and single body (9/99).
      IF (NAME(I).LT.0.AND.NAME(I).GE.-NZERO.AND.JCOMP.LE.N) THEN
*       Indentify merged ghost particle JG.
          CALL FINDJ(I1,JG,IM)
          WRITE (6,28)  NAME(I), NAME(JCOMP), NAME(JG), ECC1, PMIN, RIJ
   28     FORMAT (' HI CHAIN    NAM E1 PM RIJ ',I7,2I6,F7.3,1P,2E10.2)
          JJ = JCOMP
*       Terminate the merger in the usual way.
          KSPAIR = IPAIR
          IPHASE = 7
          CALL RESET
          ZMU = BODY(2*NPAIRS-1)*BODY(2*NPAIRS)/BODY(NTOT)
          EBCH0 = EBCH0 + ZMU*H(NPAIRS)
*       Specify chain indicator and define the two single particles.
          IPHASE = 8
          JCMAX = JG
          JCLOSE = JJ
          KSPAIR = NPAIRS
*       Set relevant variables in DELAY before terminating inner binary.
          CALL DELAY(KCHAIN,KS2)
          CALL DELAY(IPHASE,-1)
*       Initialize new chain of the 4 members JMAX, JCLOSE & KS components.
          ISUB = 0
          CALL CHAIN(ISUB)
*       Note that IPHASE = -1 now and INTGRT goes back to the beginning.
      ELSE IF (NAME(I).LT.-NZERO.OR.NAME(JCOMP).LT.0.OR.
     &        (NAME(I).LT.0.AND.JCOMP.GT.N)) THEN
*       Continue until KS termination on MERGE2 or merger with JCOMP > N.
          IPHASE = 0
      END IF
  999 CONTINUE
*
      GO TO 100
*
*       Begin check for merger of stable hierarchical configuration.
   30 RA = SEMI1*(1.0 + ECC1)
      IF (SEMI1.LT.0.0) RA = RIJ
*
*       Do not allow merger in the inner region of perturbed eccentric orbit.
      IF (RIJ.LT.SEMI1.AND.LIST(1,I1).GT.0) THEN
          IF (ECC1.GT.0.95.AND.RIJ.LT.2.0*PMIN) THEN
              GO TO 100
          END IF
      END IF
*
*       Allow temporary merger of inner part of extremely eccentric orbit.
      RFAC = 10.0*RMIN
      IF (ECC1.GT.0.99.AND.RA.GT.RFAC) THEN
          IF (RIJ.LT.0.1*SEMI1) RFAC = RA
      END IF
*
*       Increase apocentre tolerance to local scale factor for EB1 < EBS.
      RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                               (X(3,I) - RDENS(3))**2
      EBS = 0.25*EBH/SQRT(1.0 + SQRT(RI2)/RSCALE)
      IF (EB1.LT.EBS) THEN
          H2 = (RC**2 + RI2)/FLOAT(NC+10)**0.66667
          RH = 6.0*SQRT(H2/CMSEP2)
          RFAC = MAX(RFAC,RH)
      END IF
 
*       Skip merger for hyperbolic & soft binding energy or large apocentre.
*     IF (EB.GT.EBH.OR.EB1.GT.EBS.OR.RA.GT.RFAC) THEN
      IF (EB.GT.EBH.OR.EB1.GT.EBS) THEN
          GO TO 100
      ELSE
*       Set close encounter indices for hard hierarchy (redundant ?).
          ICLOSE = I
          JCLOSE = JCOMP
      END IF
      IF (SEMI1.LT.0.0.AND.RA.GT.RFAC) GO TO 100
*
*       Estimate the relative apocentre perturbations on body #I & JCOMP.
      PERT = PERT*(RA/RIJ)**3
      PERTA = PERT4*(RA/RIJ)**3
*
*       Check for circularization or collision.
      IF (KZ(27).GT.0) THEN
*       Skip merger for active phase of tidal dissipation.
          RT = 4.0*MAX(RADIUS(I1),RADIUS(I2))
          PERI = SEMI*(1.0 - ECC)
*       Include circularization test (binary may be unperturbed).
          IF (ECC.LT.0.0021.AND.KSTAR(I).EQ.0) KSTAR(I) = 20
          IF (PERI.LT.RADIUS(I1) + RADIUS(I2)) THEN
              KSPAIR = IPAIR
              CALL CMBODY(PERI,2)
              IF (IPHASE.LT.0) GO TO 100
          END IF
          IF (ECC.GT.0.1.AND.KSTAR(I).EQ.20) THEN
              WRITE (6,35)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2),
     &                      KSTAR(I), LIST(1,I1), ECC, SEMI, RT, PERI,
     &                      GAMMA(IPAIR)
   35         FORMAT (' NON-CIRCULAR    NM K* NP E A RT PERI G ',
     &                                  2I5,4I4,F7.3,1P,4E9.1)
              KSTAR(I) = 0
          END IF
          IF (PERI.LT.0.99*RT.AND.KSTAR(I).EQ.0) GO TO 100
      END IF
*
*       Ensure consistency of estimated perturbations with termination.
      PERT = PERT + PERTA
*       Allow highly eccentric outer orbit by reducing estimated PERT.
*     IF (ECC1.GT.0.98.AND.RIJ.LT.0.1*SEMI1) PERT = PERTM
      IF (PERT4.GT.GMAX.OR.PERT.GT.0.05) GO TO 100
*
*       Check for sufficient perturbers at apocentre (RP < 2*RS).
      RP = SQRT(CMSEP2)*RA
      NNB = LIST(1,I)
*       Allow for increased neighbour radius (updated in routine MERGE).
      RSI = RS(I)*MAX((FLOAT(NNBMAX)/FLOAT(NNB))**0.3333,1.0)
      IF (RP.GT.2.0*RSI) GO TO 100
*
*       Skip merger if an outer binary is fairly perturbed or not hard.
      IF (JCOMP.GT.N) THEN
          IF (GAMMA(JPAIR).GT.1.0E-03.OR.EB2.GT.EBH) GO TO 100
      END IF
*
*       Ensure the inner semi-major axis is used for subsequent tests.
   40 SEMI = -0.5*BODY(I)/H(IPAIR)
*
*       Form coefficients for stability test (Valtonen, Vistas Ast 32, 1988).
*     AM = (2.65 + ECC)*(1.0 + BODY(JCOMP)/BODY(I))**0.3333
*     FM = (2.0*BODY(JCOMP) - BODY(I))/(3.0*BODY(I))
*
*       Expand natural logarithm for small arguments.
*     IF (ABS(FM).LT.0.67) THEN
*         BM = FM*(1.0 - (0.5 - ONE3*FM)*FM)
*     ELSE
*         BM = LOG(1.0D0 + FM)
*     END IF
*
*       Adopt mass dependent criterion of Harrington (A.J. 82, 753) & Bailyn.
*     PCRIT = AM*(1.0 + 0.7*BM)*SEMI
*
*       Employ the new stability criterion (MA 1997).
      Q = BODY(JCOMP)/BODY(I)
      IF (ECC1.LT.1.0) THEN
          XFAC = (1.0 + Q)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
      ELSE
          XFAC = 40.0*(1.0 + Q)
      END IF
*       Include correction at small eccentricity (f(E) = 1.0 for now).
      FE  = 1.0
      PCRIT = 2.8*FE*XFAC**0.4*SEMI
*
*       Choose the most dominant triple in case of two binaries.
      YFAC = 1.0
      IF (JCOMP.GT.N) THEN
          SFAC = (1.0 + Q)**0.4*SEMI
          SFAC2 = (1.0 + BODY(I)/BODY(JCOMP))**0.4*SEMI2
*       Adopt 10% fudge factor with linear dependence on smallest ratio.
          YFAC = 1.0 + 0.1*MIN(SEMI2/SEMI,SEMI/SEMI2)
          IF (SFAC2.GT.SFAC) THEN
              PCRIT = PCRIT*(SFAC2/SFAC)
          END IF
      END IF
*
*       Prepare inclination evaluation for triple or widest inner binary.
      IF (JCOMP.GT.N) THEN
*       Ensure widest inner binary (swap is OK for termination or ZARE).
          IF (SEMI.LT.SEMI2) THEN
              ECC2 = (1.0 - R(JPAIR)/SEMI2)**2 +
     &                             TDOT2(JPAIR)**2/(BODY(JCOMP)*SEMI2)
              ECC = SQRT(ECC2)
              IPAIR = JPAIR
              I1 = 2*IPAIR - 1
              I2 = I1 + 1
              JJ = I
              I = JCOMP
              JCOMP = JJ
          END IF
      END IF
*
*       Resolve weakly perturbed binary (prevent X(K,I1) = X(K,I2)).
      IF (GAMMA(IPAIR).LT.GMIN.OR.X(1,I1).EQ.X(1,I2)) THEN
          CALL RESOLV(IPAIR,1)
      END IF
*
*       Copy coordinates and velocities to local variables.
      DO 42 K = 1,3
          XX(K,1) = X(K,I1)
          XX(K,2) = X(K,I2)
          XX(K,3) = X(K,JCOMP)
          VV(K,1) = XDOT(K,I1)
          VV(K,2) = XDOT(K,I2)
          VV(K,3) = XDOT(K,JCOMP)
  42  CONTINUE
*
*       Determine the inclination (OK to use c.m. instead of inner c.m.).
      CALL INCLIN(XX,VV,X(1,I),XDOT(1,I),ANGLE)
*
*       Adopt an emperical fudge factor for the inclination.
      IF (ECC1.LT.1.0) THEN
          YFAC = YFAC - 0.3*ANGLE/180.0
*       Employ additional gradual reduction above ECC1 = 0.98.
*         IF (ECC1.GT.0.98) THEN
          IF (ECC1.GT.0.96) THEN
              YFAC = YFAC - 10.0*(ECC1 - 0.96)
          END IF
      ELSE IF (RIJ.GT.RMIN) THEN
          GO TO 100
      END IF
*
*       Check whether the main perturber dominates the outer component.
      RIJ2 = (X(1,JMAX) - X(1,JCOMP))**2 + (X(2,JMAX) - X(2,JCOMP))**2 +
     &                                     (X(3,JMAX) - X(3,JCOMP))**2
      FMAX = (BODY(JMAX) + BODY(JCOMP))/RIJ2
      IF (FMAX.GT.(BODY(I) + BODY(JCOMP))/RJMIN2) GO TO 100
*
*       Include inclination angle procedure for marginal stability.
      IF (PMIN.LT.YFAC*PCRIT.AND.PMIN.GT.0.6*PCRIT) THEN
*    &    JCOMP.GT.N.AND.SEMI2.GT.10.0*SEMI) THEN
          NMARG = NMARG + 1
          YF = 0.99*PMIN*(1.0 - PERT)/PCRIT
          TK = TWOPI*SEMI*SQRT(SEMI/(BODY(I) + BODY(JCOMP)))
*         WRITE (6,44)  NMARG, ANGLE, YF, PMIN, YF*PCRIT, TK
*  44     FORMAT (' MARGINAL    # ANGLE YF PM YF*PCR TK ',
*    &                          I7,2F7.2,1P,3E10.2)
*         CALL FLUSH(6)
          IF (PMIN*(1.0 - PERT).GT.YF*PCRIT.AND.
     &        (NMARG.GT.10000.OR.NMARG*TK.GT.0.5*DTADJ)) THEN
*    &        (NMARG.GT.100000.OR.NMARG*TK.GT.0.5*DTADJ)) THEN
              YFAC = 0.99*PMIN*(1.0 - PERT)/PCRIT
              WRITE (6,45)  NMARG, NAME(I), NAME(JCOMP), ANGLE, Q,
     &                      YFAC, PMIN, PCRIT, YFAC*PCRIT
   45         FORMAT (' NEW HIERARCHY    # NAM ANGLE Q YF PM PC1 PC2 ',
     &                                   I8,2I6,F7.1,2F6.2,1P,3E10.2)
              NMARG = 0
              TK1 = TWOPI*SEMI1*SQRT(SEMI1/(BODY(I) + BODY(JCOMP)))
              WRITE (6,46)  ECC1, PERT, RP, RSI, TK1
   46         FORMAT (' WATCH!    E1 G RP RSI TK1 ',F8.4,F8.4,1P,5E10.2)
          END IF
      END IF
*
*       Check perturbed stability condition (factor 1.01 avoids switching).
      IF (PMIN*(1.0 - PERT).LT.1.01*YFAC*PCRIT) GO TO 100
*
*       Check Zare exchange stability criterion and create diagnostics.
      IF (SEMI1.GT.0.0) THEN
          CALL ZARE(I1,I2,SP)
          IF (SP.LT.1.0.AND.ANGLE.LT.10.0) THEN
              WRITE (6,48)  TIME+TOFF, Q, ECC, ECC1, SEMI, PMIN, PCRIT,
     &                      YFAC, SP
   48         FORMAT (' ZARE TEST    T Q E E1 A PM PCR YF SP ',
     &                               F8.2,F5.1,2F7.3,1P,3E9.1,0P,2F6.2)
              GO TO 100
          END IF
          WRITE (73,49)  TIME+TOFF, Q, ECC, ECC1, SEMI, PMIN, PCRIT,
     &                   SP, ANGLE, KSTAR(I)
   49     FORMAT (' STAB    T Q E E1 A PM PCR SP IN K* ',
     &                      F8.2,F5.1,2F7.3,1P,3E9.1,0P,F6.2,F7.1,I4)
          CALL FLUSH(73)
      END IF
*
*       Specify the final critical pericentre using the fudge factor.
      PCRIT = YFAC*PCRIT
*
      IF (NMERGE.EQ.MMAX) THEN
          WRITE (6,50)  NMERGE
   50     FORMAT (5X,'WARNING!   MERGER LIMIT REACHED   NMERGE =',I4)
          GO TO 100
      END IF
*
*       Skip if #JCOMP is a chain c.m. but allow bound double hierarchy.
      IF (NAME(JCOMP).EQ.0) GO TO 100
      IF (ECC1.GT.1.0.AND.MIN(NAME(I),NAME(JCOMP)).LT.0) GO TO 100
*
      DO 55 ISUB = 1,NSUB
          IF (NAME(JCOMP).EQ.NAMES(1,ISUB)) GO TO 100
   55 CONTINUE
*       Do not allow the formation of a SEPTUPLET.
      IF ((NAME(I).LT.-2*NZERO.AND.JCOMP.GT.N).OR.
     &     NAME(JCOMP).LT.-2*NZERO) GO TO 100
*
      IF (NAME(I).LT.0) THEN
*       Avoid frequent mergers due to short look-up times.
          IF (KZ(19).GT.0) THEN
              TM = MIN(TEV(I1),TEV(I2),TEV(JCOMP))
              IF (TM - TIME.LT.0.05) GO TO 100
          END IF
          CALL HIMAX2(I1,ECC,SEMI,ECC1,SEMI1,EMAX,EMIN,ZI,TG,EDAV)
          Q = BODY(JCOMP)/BODY(I)
          XFAC = (1.0 + Q)*(1.0 + EMAX)/SQRT(1.0 - EMAX)
          FE  = 1.0
          PCRIT2 = 2.8*FE*XFAC**0.4*SEMI
          IF (SEMI1*(1.0 - EMAX).LT.YFAC*PCRIT2) THEN
              PM = SEMI1*(1.0 - EMAX)
              WRITE (6,58)  ECC1, EMAX, YFAC, PM, PCRIT, YFAC*PCRIT2, TG
   58         FORMAT (' WATCH HIMAX2    E1 EX YF PM PC PC2 TG ',
     &                                  2F7.3,1P,5E10.2)
          ELSE
              PM = SEMI1*(1.0 - EMAX)
              WRITE (6,59)  ECC1, EMAX, YFAC, PM, PCRIT, YFAC*PCRIT2, TG
   59         FORMAT (' CHECK HIMAX2    E1 EX YF PMX PC PC2 TG ',
     &                                  2F7.3,1P,5E10.2)
          END IF
      END IF
*
*       Include diagnostics for double hierarchy or optional standard case.
      IF (NAME(I).LT.0.OR.NAME(JCOMP).LT.0) THEN
          IF (KZ(15).GT.1) THEN
              WHICH1 = ' MERGE2 '
              WRITE (6,20)  WHICH1, IPAIR, TIME+TOFF, H(IPAIR),R(IPAIR),
     &                      BODY(I), BODY(JCOMP), PERT4, RIJ, PMIN,
     &                      EB1/EB, LIST(1,I1)
          END IF
*       Note rare case of two hierarchies merging and identify ghost names.
          IF (NAME(I).LT.0.AND.NAME(JCOMP).LT.0) THEN
              CALL FINDJ(I1,JI,IM)
              J1 = 2*JPAIR - 1
              CALL FINDJ(J1,JJ,JM)
              WRITE (6,60)  NAME(JI), NAME(JJ), ECC, ECC1, R(JPAIR)
   60         FORMAT (' HI MERGE    NM E E1 RJ ',2I6,2F7.3,1P,E10.2)
          END IF
      ELSE IF (KZ(15).GT.1) THEN
          WHICH1 = ' MERGER '
          WRITE (6,20)  WHICH1, IPAIR, TIME+TOFF, H(IPAIR), R(IPAIR),
     &                  BODY(I), BODY(JCOMP), PERT4, RIJ, PMIN,
     &                  EB1/EB, LIST(1,I1)
      END IF
*
*       Copy RA & pair index and set indicator for calling MERGE from MAIN.
      RMAX = RA
      KSPAIR = IPAIR
      IPHASE = 6
*
*       Save KS indices and delay merger until end of block step.
      CALL DELAY(KS2,KS2)
*
  100 IF (IPHASE.NE.8) JCMAX = 0
*
      RETURN
*
      END
