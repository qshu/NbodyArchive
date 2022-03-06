      SUBROUTINE IMPACT(I)
*
*
*       Multiple collision or merger search.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      CHARACTER*8  WHICH1
      EXTERNAL SIGNAL
*
*
*       Set index of KS pair & first component of c.m. body #I.
      IPAIR = I - N
      I1 = 2*IPAIR - 1
      NCOUNT(21) = NCOUNT(21) + 1
      NTTRY = NTTRY + 1
      PERT1 = 0.0
      PERT2 = 0.0
      JCOMP = I
      NP = 0
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
              END IF
          END IF
   10 CONTINUE
*
      RDOT = (X(1,I) - X(1,JCOMP))*(XDOT(1,I) - XDOT(1,JCOMP)) +
     &       (X(2,I) - X(2,JCOMP))*(XDOT(2,I) - XDOT(2,JCOMP)) +
     &       (X(3,I) - X(3,JCOMP))*(XDOT(3,I) - XDOT(3,JCOMP))
*
*       Specify larger perturbation for optional chain regularization.
      GSTAR = GMIN
      KCHAIN = 0
      IF (KZ(30).GT.0.AND.NCH.EQ.0) THEN
          GSTAR = 100.0*GMIN
          KCHAIN = 1
      END IF
*
*       Only accept inward motion or small secondary perturbation.
      PERT3 = 2.0*R(IPAIR)**3*PERT2/BODY(I)
      IF (RDOT.GT.0.0.OR.PERT3.GT.100.0*GSTAR) GO TO 100
*
*       Skip rare case of merged binary or chain c.m. (denoted by NAME <= 0).
      IF (NAME(I).LE.0.OR.NAME(JCOMP).LE.0) GO TO 100
*
*       Include impact parameter test to distinguish different cases.
      A2 = (XDOT(1,I) - XDOT(1,JCOMP))**2 +
     &     (XDOT(2,I) - XDOT(2,JCOMP))**2 +
     &     (XDOT(3,I) - XDOT(3,JCOMP))**2
      RIJMIN = SQRT(RJMIN2)
      A3 = 2.0/RIJMIN - A2/(BODY(I) + BODY(JCOMP))
      SEMI1 = 1.0/A3
      A4 = RDOT**2/(SEMI1*(BODY(I) + BODY(JCOMP)))
      ECC1 = SQRT((1.0D0 - RIJMIN/SEMI1)**2 + A4)
      PMIN = SEMI1*(1.0D0 - ECC1)
*
*       Set semi-major axis & eccentricity of inner binary.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      SEMI0 = SEMI
      ECC2 = (1.0D0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
*
*       Form binding energy of inner & outer binary.
      EB = BODY(2*IPAIR-1)*BODY(2*IPAIR)*H(IPAIR)/BODY(I)
      EB1 = -0.5*BODY(JCOMP)*BODY(I)/SEMI1
*
*       Obtain the total perturbing force acting on body #I & JCOMP.
      CALL FPERT(I,JCOMP,NP,PERT)
*
*       Choose maximum of dominant scalar & total vectorial perturbation.
      PERT = PERT*RJMIN2/(BODY(I) + BODY(JCOMP))
      PERT4 = 2.0*RJMIN2*RIJMIN*PERT2/(BODY(I) + BODY(JCOMP))
      PERTM = MAX(PERT4,PERT)
*
*       Use combined semi-major axis for binary-binary collision.
      IF (JCOMP.GT.N) THEN
          JPAIR = JCOMP - N
          SEMI2 = -0.5D0*BODY(JCOMP)/H(JPAIR)
          J1 = 2*JPAIR - 1
          EB2 = -0.5*BODY(J1)*BODY(J1+1)/SEMI2
*       Ensure SEMI0 is smallest binary in case IPAIR denotes widest pair.
          SEMI0 = MIN(ABS(SEMI),ABS(SEMI2))
          SEMI = SEMI + SEMI2
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
          I1 = 2*IPAIR - 1
          ZMM = BODY(I1)*BODY(I1+1) + BODY(I)*BODY(JCOMP)
          IF (JCOMP.GT.N) THEN
              EBT = EBT + EB2
              ZMM = ZMM + BODY(J1)*BODY(J1+1)
          END IF
          RGRAV = ZMM/ABS(EBT)
*       Estimate length of chain (RSUM used for termination).
          RSUM = RIJMIN + SEMI
*         WRITE (6,15)  I, JCOMP, SEMI0, SEMI1, RIJMIN, RSUM, RGRAV, EBT
*  15     FORMAT (' IMPACT:   I J A0 A1 RIJ RSUM RG EBT ',2I4,6F8.4)
          IF (RSUM.GT.MIN(3.0*RGRAV,RMIN)) GO TO 30
          IF (RGRAV.GT.RMIN.OR.CMSEP2*RGRAV**2.GT.RS(I)**2) GO TO 30
      END IF
*
*       Adopt triple or four-body regularization for strong interactions.
      IF (PMIN.GT.2.0*SEMI.OR.PERTM.GT.100.0*GSTAR) GO TO 30
      IF (RIJMIN.GT.RMIN) GO TO 100
*
*       Switch to merger test for tidal dissipation configuration.
      IF (KZ(27).GT.0) THEN
          IF (SEMI0.LT.RSYNC) GO TO 30
      END IF
*
*       Specify maximum size of unperturbed motion.
      IF (PERT2.GT.0.0) THEN
          RPERT = (100.0*GSTAR*(BODY(I) + BODY(JCOMP))/(2.0*PERT2))**0.33
      ELSE
          RPERT = 10.0*SEMI
      END IF
*
*       Compare with existing subsystem of same type (if any).
      IGO = 0
      CALL SIGNAL(RPERT,IGO)
      IF (IGO.GT.0) GO TO 100
*
      WHICH1 = ' TRIPLE '
      IF (JCOMP.GT.N) WHICH1 = ' QUAD   '
      IF (KCHAIN.GT.0) WHICH1 = ' CHAIN '
*
      IF (KZ(15).GT.1) THEN
          WRITE (6,20)  WHICH1, IPAIR, TIME, H(IPAIR), R(IPAIR),
     &                  BODY(I), BODY(JCOMP), PERT4, RIJMIN, PMIN,
     &                  EB1/EB, LIST(1,2*IPAIR-1)
   20     FORMAT (/,' NEW',A8,I2,'  T =',F8.2,'  H =',F6.0,
     &              '  R =',1PE8.1,'  M =',0P2F7.4,'  G4 =',1PE8.1,
     &              '  R1 =',E8.1,'  P =',E8.1,'  E1 =',0PF6.3,
     &              '  NP =',I2)
      CALL FLUSH(6)
      END IF
*
*       Improve coordinates & velocities of the outer component.
      CALL XVPRED(JCOMP,0)
      IF (GAMMA(IPAIR).LT.0.0001) CALL XVPRED(I,0)
*
*       Save location of intruder for routine START3 in case of triple.
      JCLOSE = JCOMP
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
*       Terminate smallest pair first and reduce second index if higher.
          CALL KSTERM
          IF (KS2.GT.KSPAIR) KS2 = KS2 - 1
          KSPAIR = KS2
      END IF
*
*       Terminate binary in triple or widest binary-binary collision pair.
      CALL KSTERM
*
*       See whether chain regularization indicator should be switched on.
      IF (KCHAIN.GT.0) THEN
          IPHASE = 8
      END IF
      GO TO 100
*
*       Begin check for merger of stable hierarchical configuration.
   30 RA = SEMI1*(1.0 + ECC1)
*       Do not allow merger in the inner region of eccentric orbit.
      IF (RIJMIN.LT.SEMI1) GO TO 100
*
*       Estimate the relative apocentre perturbations on body #I & JCOMP.
      PERT = PERT*(RA/RIJMIN)**3
      PERTA = PERT4*(RA/RIJMIN)**3
*
*       Set close encounter indices for hard hierarchy (reset on merger).
      IF (EB1.LT.-0.5*BODYM*ECLOSE) THEN
          ICLOSE = I
          JCLOSE = JCOMP
      ELSE
*       Skip merger for soft binding energy or hyperbolic orbit.
          GO TO 100
      END IF
*
*       Check tidal capture option (synchronous or evolving binary orbit).
      IF (KZ(27).GT.0) THEN
*       Use basic perturbation test for synchronous case (set nominal PCRIT).
          IF (SEMI0.LT.RSYNC.AND.PERT.LT.0.25) THEN
              PCRIT = 3.0*SEMI0
              GO TO 40
          ELSE
*       Do not permit merger if predicted pericentre < 5 tidal capture radii.
              RT = 3.0*MAX(RADIUS(2*IPAIR-1),RADIUS(2*IPAIR))
              IF (SEMI0*(1.0D0 - ECC).LT.5.0*RT) GO TO 100
          END IF
      END IF
*
*       Ensure consistency of estimated perturbations with termination.
      PERT = PERT + PERTA
      IF (PERT4.GT.GMAX.OR.PERT.GT.0.25) GO TO 100
*
*       Check for sufficient perturbers at apocentre (RP < 2*RS).
      RP = SQRT(CMSEP2)*RA
      NNB = LIST(1,I)
*       Allow for increased neighbour radius (updated in routine MERGE).
      RSI = RS(I)*MAX((ZNBMAX/FLOAT(NNB))**0.3333,1.D0)
      IF (RP.GT.2.0*RSI) GO TO 100
*
*       Skip merger if an outer binary is fairly perturbed or not hard.
      IF (JCOMP.GT.N) THEN
          IF (GAMMA(JPAIR).GT.1.0E-04.OR.H(JPAIR).GT.-ECLOSE) GO TO 100
      END IF
*
*       Form coefficients for stability test (Valtonen, Vistas Ast 32, 1988).
      AM = (2.65 + ECC)*(1.0 + BODY(JCOMP)/BODY(I))**0.3333
      FM = (2.0*BODY(JCOMP) - BODY(I))/(3.0*BODY(I))
*
*       Expand natural logarithm for small arguments.
      IF (ABS(FM).LT.0.67) THEN
          BM = FM*(1.0 - (0.5 - ONE3*FM)*FM)
      ELSE
          BM = LOG(1.0D0 + FM)
      END IF
*
*       Adopt mass dependent criterion of Harrington (A.J. 82, 753) & Bailyn.
      PCRIT = AM*(1.0 + 0.7*BM)*SEMI
*       Check perturbed stability condition (PCRIT used by routine MERGE).
      IF (PMIN*(1.0 - PERT).LT.1.1*PCRIT) GO TO 100
*
*       Also check whether the main perturber dominates the outer component.
      RIJ2 = (X(1,JMAX) - X(1,JCOMP))**2 + (X(2,JMAX) - X(2,JCOMP))**2 +
     &                                     (X(3,JMAX) - X(3,JCOMP))**2
      FMAX = (BODY(JMAX) + BODY(JCOMP))/RIJ2
      IF (FMAX.GT.(BODY(I) + BODY(JCOMP))/RJMIN2) GO TO 100
*
   40 IF (NMERGE.EQ.MMAX) THEN
          WRITE (6,50)  NMERGE
   50     FORMAT (5X,'WARNING!   MERGER LIMIT REACHED   NMERGE =',I4)
          GO TO 100
      END IF
*
      WHICH1 = ' MERGER '
      IF (KZ(15).GT.1) THEN
          WRITE (6,20)  WHICH1, IPAIR, TIME, H(IPAIR), R(IPAIR),
     &                  BODY(I), BODY(JCOMP), PERT4, RIJMIN, PMIN,
     &                  EB1/EB, LIST(1,2*IPAIR-1)
      END IF
*
*       Copy RA & pair index and set indicator for calling MERGE from MAIN.
      RMAX = RA
      KSPAIR = IPAIR
      IPHASE = 6
*
  100 RETURN
*
      END



