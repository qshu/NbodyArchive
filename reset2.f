      SUBROUTINE RESET2
*
*
*      Termination of double hierarchy. 
*      --------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX)
*
*
*       Set index of disrupted pair and save output parameters.
      IPAIR = KSPAIR
      I = N + IPAIR
      E1 = BODY(2*IPAIR-1)*BODY(2*IPAIR)*H(IPAIR)/BODY(I)
      G1 = GAMMA(IPAIR)
      R1 = R(IPAIR)
      SEMI1 = -0.5*BODY(I)/H(IPAIR)
*
*       Locate current position in the merger table.
      IMERGE = 0
      DO 1 K = 1,NMERGE
          IF (NAMEM(K).EQ.NAME(I)) IMERGE = K
    1 CONTINUE
*
*       Check optional diagnostics for hierarchy.
      IF (KZ(11).GT.0.AND.KSTAR(IMERGE).LE.20) THEN
          CALL HIARCH(IPAIR)
      END IF
*
*       Save neighbours for correction procedure.
      NNB = LIST(1,I) + 1
      DO 4 L = 2,NNB
          J = LIST(L,I)
          JPERT(L) = J
    4 CONTINUE
*
*       Ensure that c.m. coordinates are known to highest order.
      CALL XVPRED(I,0)
*
*       Predict neighbour coordinates & velocities (XDOT used by FPOLY1).
      DO 5 L = 2,NNB
          J = JPERT(L)
          CALL XVPRED(J,0)
    5 CONTINUE
*
*       Obtain current coordinates & velocities and specify KS components.
      CALL RESOLV(IPAIR,2)
      ICOMP = 2*IPAIR - 1
      JCOMP = ICOMP + 1
*
*       Add outer component to neighbour list.
      JPERT(1) = JCOMP
*       Sum first part of potential energy correction due to tidal effect.
      JLIST(1) = ICOMP
      CALL NBPOT(1,NNB,POT1)
*
*       Find correct location of ghost particle using identification name.
      ICM = I
      JCOMP1 = JCOMP
      DO 10 I = 1,NTOT
          IF (BODY(I).EQ.0.0D0.AND.NAME(I).EQ.NAMEG(IMERGE)) JCOMP1 = I
   10 CONTINUE
*
*       Regularize two-body configuration if JCOMP1 cannot be identified.
      IF (JCOMP.EQ.JCOMP1) THEN
          WRITE (6,12)  IMERGE, NAMEG(IMERGE), JCOMP
   12     FORMAT (/,5X,'WARNING!    RESET2    JCOMP NOT IDENTIFIED ',
     &                       '   IM =',I3,'  NAMEG =',I6,'  JCOMP =',I6)
          GO TO 100
      END IF
*
*       Initialize basic variables for ghost and new c.m (JCOMP -> JCOMP1).
      J1 = JCOMP1
      J = JCOMP
   13 T0(J1) = TIME
      BODY(J1) = BODY(J)
      DO 14 K = 1,3
          X(K,J1) = X(K,J)
          X0(K,J1) = X(K,J)
          XDOT(K,J1) = XDOT(K,J)
          X0DOT(K,J1) = XDOT(K,J)
   14 CONTINUE
      IF (J.EQ.JCOMP) THEN
          J1 = ICM
          J = ICOMP
          GO TO 13
      END IF
*
*       Restore masses, coordinates & velocities of inner binary.
      BODY(ICOMP) = CM(1,IMERGE)
      BODY(JCOMP) = CM(2,IMERGE)
      ZM = -BODY(ICOMP)/(BODY(ICOMP) + BODY(JCOMP))
*
*       Begin with second component since ICOMP holds new c.m. variables.
      I = JCOMP
      DO 20 KCOMP = 1,2
          DO 15 K = 1,3
              X(K,I) = X(K,ICOMP) + ZM*XREL(K,IMERGE)
              X0DOT(K,I) = X0DOT(K,ICOMP) + ZM*VREL(K,IMERGE)
              XDOT(K,I) = X0DOT(K,I)
*       Note that XDOT is only needed for improved polynomials of JCOMP.
   15     CONTINUE
          I = ICOMP
          ZM = BODY(JCOMP)/(BODY(ICOMP) + BODY(JCOMP))
   20 CONTINUE
*
*       Copy KS variables for inner binary (small TDOT2 near apo/peri).
      I1 = 2*IPAIR - 1
      T0(I1) = TIME
      LIST(1,I1) = 1
      H(IPAIR) = HM(IMERGE)
      R(IPAIR) = 0.0D0
      TDOT2(IPAIR) = 0.0D0
      VI2 = 0.0
      DO 30 K = 1,4
          U(K,IPAIR) = UM(K,IMERGE)
          U0(K,IPAIR) = U(K,IPAIR)
          UDOT(K,IPAIR) = UMDOT(K,IMERGE)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
          VI2 = VI2 + UDOT(K,IPAIR)**2
   30 CONTINUE
*
*       Initialize force polynomial for outer component using resolved c.m.
      CALL FPOLY1(JCOMP1,JCOMP1,0)
      CALL FPOLY2(JCOMP1,JCOMP1,0)
*
*       See whether body #JCOMP1 should be included in NLIST.
      IF (T0(JCOMP1) + STEP(JCOMP1).LT.TLIST) THEN
          CALL NLMOD(JCOMP1,1)
      END IF
*
*       Initialize c.m. polynomials and activate inner binary.
      CALL KSIN2(2)
*
*       Evaluate stability parameter from current elements (no fudge factor).
      CALL HISTAB(IPAIR,JCOMP1,PMIN,RSTAB)
      R0(IPAIR) = RSTAB
*
*       Restore original name of inner hierarchy (c.m. NAME set in MERGE2).
      NAME(ICM) = NAME(ICM) + 2*NZERO
*
*       Rename perturber list for routine NBPOT.
      JPERT(1) = JCOMP1
*
*       Restore Roche stage indicator (any ghost c.m. is OK).
      KSTAR(ICM) = KSTARM(IMERGE)
*
*       See whether the outer component is a single or composite particle.
      POT3 = 0.0D0
      POT4 = 0.0D0
      IF (JCOMP1.LE.N) GO TO 50
*
*       Restore component masses for outer binary.
      JPAIR = JCOMP1 - N
      BODY(2*JPAIR-1) = CM(3,IMERGE)
      BODY(2*JPAIR) = CM(4,IMERGE)
*
*       Obtain coordinates & velocities of unperturbed binary components.
      CALL RESOLV(JPAIR,1)
*
*       Select new perturbers and initialize polynomials for KS motion.
      CALL KSLIST(JPAIR)
      CALL KSPOLY(JPAIR,1)
*
*       Apply tidal correction for outer binary perturbers.
      JLIST(1) = 2*JPAIR - 1
      JLIST(2) = 2*JPAIR
      CALL NBPOT(2,NNB,POT3)
      JLIST(1) = JCOMP1
      CALL NBPOT(1,NNB,POT4)
*
*       Update the merger energy.
      EB2 = BODY(2*JPAIR-1)*BODY(2*JPAIR)*H(JPAIR)/BODY(JCOMP1)
      EMERGE = EMERGE - EB2
*
      E2 = E1/EB2
      EB2 = EB2/BE(3)
      DP = POT4 - POT3
      IF (KZ(15).GT.1) THEN
          WRITE (6,45)  JPAIR, H(JPAIR), BODY(2*JPAIR-1),
     &                  BODY(2*JPAIR), E2, EB2, R(JPAIR), GAMMA(JPAIR),
     &                  DP
   45     FORMAT (' END OUTER MERGER',I4,'  H =',F7.1,'  M =',2F7.4,
     &            '  E1 =',F6.3,'  EB2 =',F6.3,'  RB2 =',1PE8.1,
     &            '  G2 =',E8.1,'  DP =',E8.1)
      END IF
*
*       Include interaction of body #ICOMP & JCOMP with perturbers.
   50 JLIST(1) = ICOMP
      JLIST(2) = JCOMP
      CALL NBPOT(2,NNB,POT2)
*
*       Form square of c.m. velocity correction due to tidal effects.
*     VI2 = X0DOT(1,ICM)**2 + X0DOT(2,ICM)**2 + X0DOT(3,ICM)**2
      DPHI = (POT2 - POT1) + (POT4 - POT3)
*     CORR = 1.0 + 2.0*DPHI/(BODY(ICM)*VI2)
*     IF (CORR.LE.0.0D0) CORR = 0.0
*
*       Adjust c.m. velocity by net tidal energy correction.
*     DO 60 K = 1,3
*         X0DOT(K,ICM) = SQRT(CORR)*X0DOT(K,ICM)
*  60 CONTINUE
*
*       Modify the merger energy to maintain conservation.
      EB = BODY(2*IPAIR-1)*BODY(2*IPAIR)*H(IPAIR)/BODY(ICM)
      EMERGE = EMERGE - EB + DPHI
*
      E1 = E1/EB
      EB = EB/BE(3)
      IF (KZ(15).GT.1) THEN
          WRITE (6,65)  IMERGE, TIME, BODY(2*IPAIR-1),
     &                  BODY(2*IPAIR), R1, SEMI1, EB, E1,
     &                  GAMMA(IPAIR), G1, NNB-1
   65     FORMAT (' END MERGE2',I3,'  T =',F8.2,'  M =',2F7.4,
     &            '  R1 =',1PE8.1,'  A1 =',E8.1,'  EB =',0PF6.3,
     &            '  E1 =',F6.3,'  GB =',1PE8.1,'  G =',0PF6.3,
     &            '  NB =',I3)
      END IF
*
*       Reduce merger counter and update tables (unless last or only pair).
   70 NMERGE = NMERGE - 1
      DO 80 L = IMERGE,NMERGE
          L1 = L + 1
          HM(L) = HM(L1)
          NAMEG(L) = NAMEG(L1)
          NAMEM(L) = NAMEM(L1)
          KSTARM(L) = KSTARM(L1)
          DO 74 K = 1,3
              XREL(K,L) = XREL(K,L1)
              VREL(K,L) = VREL(K,L1)
   74     CONTINUE
          DO 75 K = 1,4
              CM(K,L) = CM(K,L1)
              UM(K,L) = UM(K,L1)
              UMDOT(K,L) = UMDOT(K,L1)
   75     CONTINUE
   80 CONTINUE
*
*       Examine merger list for possible escapers (retain double mergers).
      DO 90 L = 1,NMERGE
          DO 85 J = 1,NPAIRS
              IF (NAMEM(L).EQ.NAME(N+J).OR.
     &            NAMEM(L).EQ.NAME(N+J) + 2*NZERO) GO TO 90
   85     CONTINUE
*       Remove tables for any merger not identified.
          IMERGE = L
          GO TO 70
   90 CONTINUE
*
*       Set IPHASE < 0 for new NLIST.
      IPHASE = -1
*
  100 RETURN
*
      END
