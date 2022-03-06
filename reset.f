      SUBROUTINE RESET
*
*
*       Restore hierarchical configuration.
*       -----------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX)
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
      DO 1 K = 1,NMERGE
          IF (NAMEM(K).EQ.NAME(I)) IMERGE = K
    1 CONTINUE
*
*       Identify the merger pair index corrected for any escapers.
      JMERGE = 0
      DO 2 JPAIR = 1,IPAIR
          IF (NAME(N+JPAIR).LT.0) JMERGE = JMERGE + 1
    2 CONTINUE
*
*       Save neighbours for correction procedure and rename if moved up.
      NNB = LIST(1,I) + 1
      DO 4 L = 2,NNB
          J = LIST(L,I)
          IF (J.GT.I) J = J - 1
          IF (J.LE.2*NPAIRS.AND.J.GT.2*IPAIR) J = J - 2
          JPERT(L) = J
    4 CONTINUE
*
*       Ensure that c.m. coordinates are known to highest order.
      CALL XVPRED(I,0)
*
*       Restore original c.m. name and terminate outer KS pair.
      NAME(I) = NZERO - NAME(I)
      CALL KSTERM
*
*       Predict neighbour coordinates & velocities (XDOT used by FPOLY1).
      DO 5 L = 2,NNB
          J = JPERT(L)
          CALL XVPRED(J,0)
    5 CONTINUE
*
*       Add outer component to perturber list and set old dipole range.
      JPERT(1) = JCOMP
      RCRIT2 = CMSEP2*(XREL(1,IMERGE)**2 + XREL(2,IMERGE)**2 +
     &                                                XREL(3,IMERGE)**2)
*
*       Sum first part of potential energy correction due to tidal effect.
      JLIST(1) = ICOMP
      CALL NBPOT(1,NNB,POT1)
*
*       Find the nearest neighbour and reduce steps of active perturbers.
      RJMIN2 = 100.0
      JMIN = N
      DO 8 L = 1,NNB
          J = JPERT(L)
          RIJ2 = (X(1,J) - X(1,ICOMP))**2 + (X(2,J) - X(2,ICOMP))**2 +
     &                                      (X(3,J) - X(3,ICOMP))**2
*       Identify the closest perturber for dominant force calculation.
          IF (RIJ2.LT.RJMIN2.AND.J.NE.JCOMP) THEN
              RJMIN2 = RIJ2
              JMIN = J
          END IF
*         IF (RIJ2.LT.RCRIT2) THEN
*       Reduce step of inner binary perturbers (c.m. approximation used).
*             STEP(J) = MAX(0.5D0*STEP(J),TIME - T0(J))
*       Include particle #J in time-step list unless already present.
*             IF (T0(J) + STEP(J).LT.TLIST) THEN
*                 CALL NLMOD(J,1)
*             END IF
*         END IF
    8 CONTINUE
*
*       Find correct location of ghost particle using identification name.
      JCOMP1 = JCOMP
*       Note that ghost may be saved in an old binary c.m. (search NTOT).
      DO 10 I = 1,NTOT
          IF (BODY(I).EQ.0.0D0.AND.NAME(I).EQ.NAMEG(IMERGE)) JCOMP = I
   10 CONTINUE
*
*       Regularize two-body configuration if JCOMP cannot be identified.
      IF (JCOMP.EQ.JCOMP1) THEN
          WRITE (6,12)  NAMEG(IMERGE)
   12     FORMAT (/,5X,'WARNING!   JCOMP NOT IDENTIFIED IN RESET',
     &                                                  '   NAMEG =',I5)
*       Set IPHASE = -1 to follow standard proceure in routine KSINIT.
          IPHASE = -1
          CALL KSREG
          GO TO 70
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
*       Add #JCOMP to neighbour lists containing ICOMP (KSREG sets c.m.).
      JLIST(1) = JCOMP
      CALL NBREST(ICOMP,1,NNB)
*
*       Ensure that any rare case of missing second component is included.
      DO 30 I = 1,NTOT
          NNB1 = LIST(1,I) + 1
          DO 25 L = 2,NNB1
              IF (LIST(L,I).GT.ICOMP) GO TO 30
*
*       Now see whether #JCOMP has already been added.
              DO 22 K = L,NNB1
                  IF (LIST(K,I).EQ.JCOMP) GO TO 30
   22         CONTINUE
*
*       Specify index #I and add body #JCOMP as above.
              JPERT(1) = I
              CALL NBREST(ICOMP,1,1)
   25     CONTINUE
   30 CONTINUE
*
*       Set dominant F & FDOT on inner components including main intruder.
      JLIST(1) = ICOMP
      JLIST(2) = JCOMP
      JLIST(3) = JCOMP1
      JLIST(4) = JMIN
*
      CALL FCLOSE(ICOMP,4)
      CALL FCLOSE(JCOMP,4)
*
*       Initialize force polynomials for outer component using resolved c.m.
      CALL FPOLY1(JCOMP1,JCOMP1,0)
      CALL FPOLY2(JCOMP1,JCOMP1,0)
*
*       See whether body #JCOMP1 should be include in NLIST.
*     IF (T0(JCOMP1) + STEP(JCOMP1).LT.TLIST) THEN
*         CALL NLMOD(JCOMP1,1)
*     END IF
*
*       Rename perturber list for routine NBPOT.
      JPERT(1) = JCOMP
*
*       Copy basic KS variables for inner binary (small TDOT2 near apo/peri).
      JP1 = NPAIRS + 1
      H(JP1) = HM(IMERGE)
      R(JP1) = 0.0D0
      TDOT2(JP1) = 0.0D0
      DO 40 K = 1,4
          U(K,JP1) = UM(K,IMERGE)
          U0(K,JP1) = U(K,JP1)
          UDOT(K,JP1) = UMDOT(K,IMERGE)
          R(JP1) = R(JP1) + U(K,JP1)**2
   40 CONTINUE
*
*       Save ghost index and reactivate inner binary (JCOMP <-> JCOMP1).
      JCOMP1 = JCOMP
      CALL KSREG
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
      VI2 = X0DOT(1,NTOT)**2 + X0DOT(2,NTOT)**2 + X0DOT(3,NTOT)**2
      DPHI = (POT2 - POT1) + (POT4 - POT3)
      CORR = 1.0 + 2.0*DPHI/(BODY(NTOT)*VI2)
      IF (CORR.LE.0.0D0) CORR = 0.0
*
*       Adjust c.m. velocity by net tidal energy correction.
      DO 60 K = 1,3
          X0DOT(K,NTOT) = SQRT(CORR)*X0DOT(K,NTOT)
   60 CONTINUE
*
*       Modify the merger energy to maintain conservation.
      EB = BODY(2*NPAIRS-1)*BODY(2*NPAIRS)*H(NPAIRS)/BODY(NTOT)
      EMERGE = EMERGE - EB
*       Note that EMERGE may contain escaped mergers.
*
      E1 = E1/EB
      EB = EB/BE(3)
      IF (KZ(15).GT.1) THEN
          WRITE (6,65)  IMERGE, TIME, BODY(2*NPAIRS-1),
     &                  BODY(2*NPAIRS), R1, SEMI1, EB, E1,
     &                  GAMMA(NPAIRS), G1, NNB
   65     FORMAT (' END MERGER',I3,'  T =',F7.2,'  M =',2F7.4,
     &            '  R1 =',1PE8.1,'  A1 =',E8.1,'  EB =',0PF6.3,
     &            '  E1 =',F6.3,'  GB =',1PE8.1,'  G =',0PF6.3,
     &            '  NB =',I3)
      END IF
*
*       Reduce merger counter and update tables (unless last or only pair).
   70 NMERGE = NMERGE - 1
      IF (JMERGE.LE.NMERGE) THEN
*       Old escapers above IMERGE are removed first if IMERGE > JMERGE.
          DO 80 L = JMERGE,NMERGE
              HM(L) = HM(L+1)
              NAMEG(L) = NAMEG(L+1)
              NAMEM(L) = NAMEM(L+1)
              DO 74 K = 1,3
                  XREL(K,L) = XREL(K,L+1)
                  VREL(K,L) = VREL(K,L+1)
   74         CONTINUE
              DO 76 K = 1,4
                  CM(K,L) = CM(K,L+1)
                  UM(K,L) = UM(K,L+1)
                  UMDOT(K,L) = UMDOT(K,L+1)
   76         CONTINUE
   80     CONTINUE
*
*       See whether another update is needed (current c.m. removed last).
          IMERGE = IMERGE - 1
          IF (IMERGE.GE.JMERGE) GO TO 70
      END IF
*
*       Set phase indicator = -1 for new NLIST in routine INTGRT.
      IPHASE = -1
*
      RETURN
*
      END
