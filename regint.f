      SUBROUTINE REGINT(I,KLIST)
*
*       Regular integration.
*       --------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX,1)
*       Calculate potential with little extra cost.
      COMMON/POTENT/PHII(NMAX),PHIR(NMAX),PHIR1(NMAX)
      REAL*8  W0(4),W1(4),W2(4),W3(4)
      REAL*8  XI(3),XIDOT(3),FIRR(3),FREG(3),DV(3),FD(3),FDR(3)
*
      INTEGER IKSRL,JKSRL,KLIST(LMAX)
      COMMON/LOCKSR/IKSRL(NMAX),JKSRL(NMAX)
*
      DO 6 K = 1,3
          XI(K) = X(K,I)
          XIDOT(K) = XDOT(K,I)
   6  CONTINUE
*
*       Copy uncorrected X and set time-step & central distance.
      NNB0 = KLIST(1)
      DTR = TIME - T0R(I)
      IRSKIP = 0
      RI2 = (XI(1) - RDENS(1))**2 + (XI(2) - RDENS(2))**2 +
     &                              (XI(3) - RDENS(3))**2
*
*       Obtain irregular & regular force and determine current neighbours.
      RSPH = RS(I)
      RS2 = RSPH**2
*       Start count at 2 and subtract 1 at the end to avoid IKSRL(NNB+1).
    1 NNB = 1
*
*       Initialize scalars for force & first derivative.
      DO 5 K = 1,3
          FIRR(K) = 0.0D0
          FREG(K) = 0.0D0
          FD(K) = 0.0
          FDR(K) = 0.0
    5 CONTINUE
      PHII(I) = 0.D0
      PHIR(I) = 0.D0
      PHIR1(I) = 0.D0
*
*       Choose appropriate force loop for single particle or c.m. body.
      IF (I.GT.N) THEN
*       See whether perturbation allows single particle approximation.
          IF (GAMMA(I-N).GE.GMIN) THEN
*       Obtain total force on c.m. particle.
          CALL CMFREG(I,XI,XIDOT,RS2,NNB,FIRR,FREG,FD,FDR)
              GO TO 20
          END IF
      END IF
*
*       Perform fast force loop over single particles.
          call cputim(tt1)
      DO 10 J = IFIRST,N
*RSP
          IF (J.EQ.I) GO TO 10
*RSP
          A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
*       Predicted coordinates avoids spurious force differences.
          DV(1) = XDOT(1,J) - XIDOT(1)
          DV(2) = XDOT(2,J) - XIDOT(2)
          DV(3) = XDOT(3,J) - XIDOT(3)
*
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDP = A1*DV(1) + A2*DV(2) + A3*DV(3)
          DRDV = 3.0*DRDP*DR2I
*
*       See whether the distance exceeds the outer shell radius.
          IF (RIJ2.GT.RS2) GO TO 8
*RSP          IF (J.EQ.I) GO TO 10
*
*       Increase neighbour counter and obtain current irregular force.
          NNB = NNB + 1
          IKSRL(NNB) = J
          FIRR(1) = FIRR(1) + A1*DR3I
          FIRR(2) = FIRR(2) + A2*DR3I
          FIRR(3) = FIRR(3) + A3*DR3I
          FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
          FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
          FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
*       Obtain potential.
          PHII(I) = PHII(I) - DR3I*RIJ2
          GO TO 10
*
*       Obtain the regular force.
    8     FREG(1) = FREG(1) + A1*DR3I
          FREG(2) = FREG(2) + A2*DR3I
          FREG(3) = FREG(3) + A3*DR3I
          FDR(1) = FDR(1) + (DV(1) - A1*DRDV)*DR3I
          FDR(2) = FDR(2) + (DV(2) - A2*DRDV)*DR3I
          FDR(3) = FDR(3) + (DV(3) - A3*DRDV)*DR3I
*       Obtain potential and derivative.
          PHIR(I) = PHIR(I) - DR3I*RIJ2
          PHIR1(I) = PHIR1(I) + DRDP*DR3I
   10 CONTINUE
          call cputim(tt2)
          ttfrc = ttfrc + (tt2-tt1)*60.
*
*       Add any contributions from regularized c.m. particles.
      IF (NPAIRS.GT.0) THEN
          CALL CMFREG(I,XI,XIDOT,RS2,NNB,FIRR,FREG,FD,FDR)
      END IF
*
*       Include treatment for regularized clump.
      IF (NCH.GT.0) THEN
*       Distinguish between chain c.m. and any other particle.
          IF (NAME(I).EQ.0) THEN
              CALL CHFIRR(I,1,XI,XIDOT,FIRR,FD)
          ELSE
*       Search the chain perturber list for #I.
              DO 15 L = 2,NNB+1
                  J = IKSRL(L)
                  IF (J.GT.ICH) GO TO 20
                  IF (J.EQ.ICH) CALL FCHAIN(I,1,XI,XIDOT,FIRR,FD)
   15         CONTINUE
          END IF
      END IF
*
*       See whether an external force should be added.
   20 CONTINUE
*
      IF (KZ(14).NE.0) THEN
          CALL XTRNLF(I,XI,XIDOT,FIRR,FREG,FD,FDR,1)
      END IF
*
*       Check whether cloud forces should be included.
*     IF (KZ(13).NE.0) THEN
*         CALL FCLOUD(I,FREG,FDR,1)
*     END IF
*
      NNB = NNB - 1
      IF (NNB.EQ.0) THEN
*       Double the neighbour sphere and try again unless RS > 50*RSCALE.
          IF (RSPH.GT.50.0*RSCALE) THEN
              IRSKIP = 1
*       Assume small mass at centre for rare case of no neighbours.
              FIJ = 0.01*BODYM/(RI2*SQRT(RI2))
              DO 25 K = 1,3
                  FIRR(K) = FIRR(K) - FIJ*XI(K)
                  FD(K) = FD(K) - FIJ*XIDOT(K)
   25         CONTINUE
              KLIST(1) = 0
              GO TO 50
          ELSE
              RS2 = 1.59*RS2
          END IF
          RSPH = SQRT(RS2)
          NCOUNT(6) = NCOUNT(6) + 1
          NBVOID = NBVOID + 1
          IF (RSPH.GT.10.0*RSCALE) IRSKIP = 1
          GO TO 1
      END IF
*
*       Restrict neighbour number < NNBMAX to permit one normal addition.
      IF (NNB.LT.NNBMAX) GO TO 40
*
*       Reduce search radius by cube root of conservative volume factor.
   30 NNB2 = 0.8*NNBMAX
      A1 = FLOAT(NNB2)/FLOAT(NNB)
      IF (RSPH.GT.5.0*RSCALE) THEN
          A1 = MIN(5.0*A1,0.9D0)
          IRSKIP = 1
      END IF
      RS2 = RS2*A1**0.66667
      RSPH = SQRT(RS2)
      NNB1 = 0
*
      DO 35 L = 1,NNB
          J = IKSRL(L+1)
          IF (L + NNB2.GT.NNB + NNB1) GO TO 32
*       Sum of neighbours (NNB1) & those left (NNB+1-L) set to NNB2.
          A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
          DV(1) = XDOT(1,J) - XIDOT(1)
          DV(2) = XDOT(2,J) - XIDOT(2)
          DV(3) = XDOT(3,J) - XIDOT(3)
*
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDP = A1*DV(1) + A2*DV(2) + A3*DV(3)
          DRDV = 3.0*DRDP*DR2I
          IF (RIJ2.GT.RS2) GO TO 34
*
   32     NNB1 = NNB1 + 1
          JKSRL(NNB1+1) = J
          GO TO 35
*
*       Subtract neighbour force included above and add to regular force.
   34     FIRR(1) = FIRR(1) - A1*DR3I
          FIRR(2) = FIRR(2) - A2*DR3I
          FIRR(3) = FIRR(3) - A3*DR3I
          FD(1) = FD(1) - (DV(1) - A1*DRDV)*DR3I
          FD(2) = FD(2) - (DV(2) - A2*DRDV)*DR3I
          FD(3) = FD(3) - (DV(3) - A3*DRDV)*DR3I
          FREG(1) = FREG(1) + A1*DR3I
          FREG(2) = FREG(2) + A2*DR3I
          FREG(3) = FREG(3) + A3*DR3I
          FDR(1) = FDR(1) + (DV(1) - A1*DRDV)*DR3I
          FDR(2) = FDR(2) + (DV(2) - A2*DRDV)*DR3I
          FDR(3) = FDR(3) + (DV(3) - A3*DRDV)*DR3I
*       Obtain potential and derivative.
          PHII(I) = PHII(I) - DR3I*RIJ2
          PHIR(I) = PHIR(I) - DR3I*RIJ2
          PHIR1(I) = PHIR1(I) + DRDP*DR3I
   35 CONTINUE
*
      DO 38 L = 2,NNB1+1
          IKSRL(L) = JKSRL(L)
   38 CONTINUE
      NNB = NNB1
      NCOUNT(5) = NCOUNT(5) + 1
      NBFULL = NBFULL + 1
*       See whether to reduce NNB further.
      IF (NNB.GE.NNBMAX) GO TO 30
*
*       Stabilize NNB between ZNBMIN & ZNBMAX by square root of contrast.
*       Include optional stabilization to increase neighbour number.
*       Take input parameter NNBOPT as optimal neighbour number (R.Sp.)
*       Note that it substitutes input parameter NNBMAX, which
*       is now a parameter NNBMAX=LMAX-3
*
   40 CONTINUE
*
      FAC = 1.D0
      IF (KZ(40).GT.1) THEN
          FAC = 1.0 + (FLOAT(NNBOPT) - FLOAT(NNB))/FLOAT(NNBOPT)
      END IF
*
*     A3 = ALPHA*FAC*SQRT(FLOAT(NNB)*RSPH)/RS2
      A3 = FLOAT(NNB)*FAC
      A3 = MIN(A3,0.9*ZNBMAX)
      A4 = MAX(A3,ZNBMIN)/FLOAT(NNB)
*
*       Include inertial factor to prevent resonance oscillations in RS.
      IF ((A3 - FLOAT(NNB0))*(A3 - FLOAT(NNB)).LT.0.0) A4 = SQRT(A4)
*
*       Modify volume ratio by radial velocity factor outside the core.
      IF (RI2.GT.RC2) THEN
          RIDOT = (XI(1) - RDENS(1))*XIDOT(1) +
     &            (XI(2) - RDENS(2))*XIDOT(2) +
     &            (XI(3) - RDENS(3))*XIDOT(3)
          A4 = A4*(1.0 + RIDOT*DTR/RI2)
      END IF
*
*       See whether neighbour radius of c.m. body should be increased.
      IF (I.GT.N) THEN
*       Set perturber range (soft binaries & H > 0 have short duration).
          A2 = 100.0*BODY(I)/ABS(H(I-N))
*       Stabilize NNB on ZNBMAX if too few perturbers.
          IF (A2.GT.RSPH) THEN
              A3 = MAX(1.0 - FLOAT(NNB)/ZNBMAX,0.0D0)
*       Modify volume ratio by approximate square root factor.
              A4 = 1.0 + 0.5*A3
          END IF
      END IF
*
*       Limit change of RS for small steps (RSFAC = MAX(25/TCR,1.5*VC/RSMIN).
      A5 = MIN(RSFAC*DTR,0.25D0)
*       Restrict volume ratio by inertial factor of 25 per cent either way.
      A4 = A4 - 1.0
      IF (A4.GT.A5) THEN
          A4 = A5
      ELSE
          A4 = MAX(A4,-A5)
      END IF
*
*       Modify neighbour sphere radius by volume factor.
      IF (IRSKIP.EQ.0) THEN
          IF (RSPH.GT.50.0*RSCALE) GO TO 50
          A3 = ONE3*A4
          A1 = 1.0 + A3 - A3*A3
*       Second-order cube root expansion (maximum volume error < 0.3 %).
          IF (RSPH.GT.5.0*RSCALE) A1 = SQRT(A1)
*       Inertial factor for distant particle avoids oscillations in RS.
          RSPH = A1*RSPH
      END IF
*
*       Calculate the radial velocity with respect to at most 3 neighbours.
      IF (NNB.LE.3) THEN
          A1 = 2.0*RSPH
*
          DO 45 L = 1,NNB
              J = IKSRL(L+1)
              RIJ = SQRT((XI(1) - X(1,J))**2 + (XI(2) - X(2,J))**2 +
     &                                         (XI(3) - X(3,J))**2)
              RSDOT = ((XI(1) - X(1,J))*(XIDOT(1) - XDOT(1,J)) +
     &                 (XI(2) - X(2,J))*(XIDOT(2) - XDOT(2,J)) +
     &                 (XI(3) - X(3,J))*(XIDOT(3) - XDOT(3,J)))/RIJ
*       Find smallest neighbour distance assuming constant regular step.
              A1 = MIN(A1,RIJ + RSDOT*DTR)
   45     CONTINUE
*
*       Increase neighbour sphere if all members are leaving inner region.
          RSPH = MAX(A1,1.1*RSPH)
      END IF
*
*       Check optional procedures for adding neighbours.
      IF (KZ(18).EQ.0) GO TO 50
      IF (KZ(18).EQ.1.AND.LISTV(1).EQ.0) GO TO 50
      CALL CHECKL(I,NNB,XI,XIDOT,RS2,FIRR,FREG,FD,FDR)
*
*       Find loss or gain of neighbours at the same time.
   50 NBLOSS = 0
      NBGAIN = 0
*
*       Check case of zero old or new membership (skip if both are zero).
      IF (NNB0.EQ.0) THEN
          IF (NNB.EQ.0) GO TO 70
          KLIST(2) = 0
      END IF
*
      JMIN = 0
      L = 2
      LG = 2
*       Set termination value in IKSRL(NNB+2) and save last list member.
      IKSRL(NNB+2) = NTOT + 1
      IKSRL(1) = KLIST(NNB0+1)
*
*       Compare old and new list members in locations L & LG.
   56 IF (KLIST(L).EQ.IKSRL(LG)) GO TO 58
*
*       Now check whether inequality means gain or loss.
      IF (KLIST(L).GE.IKSRL(LG)) THEN
          NBGAIN = NBGAIN + 1
          JKSRL(NNB0+NBGAIN) = IKSRL(LG)
*       Number of neighbour losses can at most be NNB0.
          L = L - 1
*       The same location will be searched again after increasing L below.
      ELSE
          NBLOSS = NBLOSS + 1
          J = KLIST(L)
          JKSRL(NBLOSS) = J
*       Check SMIN step indicator (rare case permits fast skip below).
          IF (STEP(J).LT.SMIN) JMIN = J
          LG = LG - 1
      END IF
*
*       See whether the last location has been checked.
   58 IF (L.LE.NNB0) THEN
          L = L + 1
          LG = LG + 1
*       Last value of second search index is NNB + 2 which holds NTOT + 1.
          GO TO 56
      ELSE IF (LG.LE.NNB) THEN
          LG = LG + 1
          KLIST(L) = NTOT + 1
*       Last location of list holds termination value (saved in IKSRL(1)).
          GO TO 56
      END IF
*
*       See whether any old neighbour with small step should be retained.
      IF (JMIN.EQ.0) GO TO 70
*
      K = 1
   60 IF (NNB.GT.NNBMAX.OR.I.GT.N) GO TO 70
      J = JKSRL(K)
*       A single regularized component will be replaced by the c.m.
      IF (STEP(J).GT.SMIN.OR.J.LT.IFIRST.OR.J.GT.N) GO TO 68
*       Retain old neighbour inside 2*RS to avoid large correction terms.
      RIJ2 = (XI(1) - X(1,J))**2 + (XI(2) - X(2,J))**2 +
     &                             (XI(3) - X(3,J))**2
      IF (RIJ2.GT.4.0*RS2) GO TO 68
*
      L = NNB + 1
   62 IF (IKSRL(L).LT.J) GO TO 64
      IKSRL(L+1) = IKSRL(L)
      L = L - 1
      IF (L.GT.1) GO TO 62
*
*       Save index of body #J and update NNB & NBLOSS.
   64 IKSRL(L+1) = J
      NNB = NNB + 1
      NBLOSS = NBLOSS - 1
*       Restore last old neighbour in case NBLOSS = 0 at end of search.
      KLIST(NNB0+1) = IKSRL(1)
      NBSMIN = NBSMIN + 1
      NCOUNT(11) = NCOUNT(11) + 1
*
*       Perform correction to irregular and regular force components.
      A1 = X(1,J) - XI(1)
      A2 = X(2,J) - XI(2)
      A3 = X(3,J) - XI(3)
      DV(1) = XDOT(1,J) - XIDOT(1)
      DV(2) = XDOT(2,J) - XIDOT(2)
      DV(3) = XDOT(3,J) - XIDOT(3)
*
      RIJ2 = A1*A1 + A2*A2 + A3*A3
      DR2I = 1.0/RIJ2
      DR3I = BODY(J)*DR2I*SQRT(DR2I)
      DRDP = A1*DV(1) + A2*DV(2) + A3*DV(3)
      DRDV = 3.0*DRDP*DR2I
*
      FIRR(1) = FIRR(1) + A1*DR3I
      FIRR(2) = FIRR(2) + A2*DR3I
      FIRR(3) = FIRR(3) + A3*DR3I
      FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
      FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
      FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
      FREG(1) = FREG(1) - A1*DR3I
      FREG(2) = FREG(2) - A2*DR3I
      FREG(3) = FREG(3) - A3*DR3I
      FDR(1) = FDR(1) - (DV(1) - A1*DRDV)*DR3I
      FDR(2) = FDR(2) - (DV(2) - A2*DRDV)*DR3I
      FDR(3) = FDR(3) - (DV(3) - A3*DRDV)*DR3I
*       Obtain potential and derivative.
      PHII(I) = PHII(I) - DR3I*RIJ2
      PHIR(I) = PHIR(I) - DR3I*RIJ2
      PHIR1(I) = PHIR1(I) + DRDP*DR3I
*
*       Remove body #J from JLIST unless it is the last or only member.
      IF (K.GT.NBLOSS) GO TO 70
      DO 66 L = K,NBLOSS
          JKSRL(L) = JKSRL(L+1)
   66 CONTINUE
*       Index of last body to be moved up is L = NBLOSS + 1.
      K = K - 1
*       Check the same location again since a new body has appeared.
   68 K = K + 1
*       Last member to be considered is in location K = NBLOSS.
      IF (K.LE.NBLOSS) GO TO 60
*
*       Form time-step factors and update regular force time.
   70 DTSQ = DTR**2
      DT6 = 6.0/(DTR*DTSQ)
      DT2 = 2.0/DTSQ
      DTSQ12 = ONE12*DTSQ
      DTR13 = ONE3*DTR
      T0R(I) = TIME
*
      TIRR = STEP(I)
*       Suppress the corrector for large time-step ratios (experimental).
      IF (DTR.GT.50.0*TIRR) THEN
          DTR13 = 0.0
          DTSQ12 = 0.0
      END IF
*
*       Include the corrector and save F, FI, FR & polynomial derivatives.
      DO 75 K = 1,3
*       Subtract change of neighbour force to form actual first derivative.
          DFR = FR(K,I) - (FIRR(K) - FI(K,I)) - FREG(K)
          FDR0 = FDR(K) - (FIDOT(K,I) - FD(K))
*
          FRD = FRDOT(K,I)
	  SUM = FRD + FDR0
	  AT3 = 2.0*DFR + DTR*SUM
	  BT2 = -3.0*DFR - DTR*(SUM + FRD)
*
          X0(K,I) = X0(K,I) + (0.6*AT3 + BT2)*DTSQ12
          X0DOT(K,I) = X0DOT(K,I) + (0.75*AT3 + BT2)*DTR13
*
          FI(K,I) = FIRR(K)
          FR(K,I) = FREG(K)
          FIDOT(K,I) = FD(K)
          FRDOT(K,I) = FDR(K)
*
          D0(K,I) = FIRR(K)
          D0R(K,I) = FREG(K)
          D0R(K,I) = FREG(K) - (FI(K,I) - FIRR(K))
          D1R(K,I) = FDR0
          D3R(K,I) = AT3*DT6
          D2R(K,I) = (3.0*AT3 + BT2)*DT2
   75 CONTINUE
*
*       Correct force polynomials due to neighbour changes (KZ(38) or I > N).
      IF (KZ(38).GT.0.OR.I.GT.N) THEN
      CALL FPCORR(I,NBLOSS,NBGAIN,XI,XIDOT,FIRR,FREG,FD,FDR,KLIST)
      END IF
*
*       Copy new neighbour list if different from old list.
      IF (NBLOSS + NBGAIN.GT.0) THEN
          KLIST(1) = NNB
          DO 80 L = 2,NNB+1
              KLIST(L) = IKSRL(L)
   80     CONTINUE
      END IF
*
*       Check for boundary reflection (RI2 < 0 denotes new polynomials set).
*     IF (KZ(29).GT.0.AND.RI2.GT.RSPH2) THEN
*         CALL REFLCT(I,RI2)
*         IF (RI2.LT.0.0) GO TO 120
*     END IF
*
*       Obtain new regular integration step using composite expression.
*       STEPR = (ETAR*(F*F2DOT + FDOT**2)/(FDOT*F3DOT + F2DOT**2))**0.5.
      DO 100 K = 1,3
          W1(K) = FDR(K)
          W2(K) = D2R(K,I)
          W3(K) = D3R(K,I)
  100 CONTINUE
*
      W0(4) = FREG(1)**2 + FREG(2)**2 + FREG(3)**2
      W1(4) = W1(1)**2 + W1(2)**2 + W1(3)**2
      W2(4) = W2(1)**2 + W2(2)**2 + W2(3)**2
      W3(4) = W3(1)**2 + W3(2)**2 + W3(3)**2
*
*       Form new step by relative criterion (extra SQRT for large F3DOT).
      IF (W3(4).LT.1.0E+20) THEN
          W0(1) = (SQRT(W0(4)*W2(4)) + W1(4))/
     &                                       (SQRT(W1(4)*W3(4)) + W2(4))
      ELSE
          W0(1) = (SQRT(W0(4)*W2(4)) + W1(4))/
     &                                 (SQRT(W1(4))*SQRT(W3(4)) + W2(4))
      END IF
      W0(1) = ETAR*W0(1)
      TTMP = SQRT(W0(1))
*       Winston Sweatman's suggestion
*     DVV = (XDOT(1,I)-X0DOT(1,I))**2 + (XDOT(2,I)-X0DOT(2,I))**2 +
*    &     (XDOT(3,I)-X0DOT(3,I))**2
*     FFD = FREG(1)**2 + FREG(2)**2 + FREG(3)**2
*     ETARW = ETAR
*     TTMPW = ETARW*DVV*BODY(I)/FFD
*
*     PRINT*,' Reg I=',I,' TTMP,TTMPW,RATIO=',
*    &  TTMP,TTMPW,TTMP/TTMPW
*
*     IF(TTMP.GT.TTMPW)THEN
*     IGT = IGT + 1
*     ELSE
*     ILE = ILE + 1
*     END IF
*     IF(MOD(IGT+ILE,100).EQ.0)PRINT*,' irr IGT,ILE=',IGT,ILE
*
*     TTMP = MAX(TTMPW,TTMP)
*     TTR = TSTEP(FREG,FDR,D2R(1,I),D3R(1,I),ETAR)
*
*       Adopt FAC*MIN(FREG,FIRR) (or tidal force) for convergence test.
      FAC = MIN(ETAR,0.04D0)
      IF (TIDAL(1).EQ.0.0D0) THEN
          FI2 = FIRR(1)**2 + FIRR(2)**2 + FIRR(3)**2
          W0(4) = FAC**2*MIN(DBLE(W0(4)),FI2)
      ELSE
          W0(1) = (TIDAL(1)*XI(1))**2
          W0(4) = FAC**2*MAX(W0(4),W0(1))
      END IF
*
*       Obtain regular force change using twice the predicted step.
      DTC = 2.0*TTMP
      S2 = 0.5*DTC
      S3 = ONE3*DTC
      W0(1) = 0.0
      DO 105 K = 1,3
          W0(2) = ((W3(K)*S3 + W2(K))*S2 + W1(K))*DTC
          W0(1) = W0(1) + W0(2)**2
  105 CONTINUE
*
*       See whether regular step can be increased by factor 2.
      IF (W0(1).LT.W0(4)) THEN
          TTMP = DTC
          NRCONV = NRCONV + 1
          NCOUNT(18) = NCOUNT(18) + 1
      END IF
*
*       Impose a smooth step reduction inside compact core.
      IF (NC.LT.50.AND.RI2.LT.RC2) THEN
          TTMP = TTMP*MIN(1.0D0,0.5D0*(1.0D0 + RI2*RC2IN))
      END IF
*
      TREG = STEPR(I)
*       Select discrete value (increased by 2, decreased by 2 or unchanged).
		IF (TTMP .GT. 2.0*TREG) THEN
			IF (DMOD(TIME,2.0*TREG) .EQ. 0.0D0) THEN
              TTMP = MIN(2.0*TREG,1.D0)
          ELSE
              TTMP = TREG
          END IF
      ELSE IF (TTMP .LT. TREG) THEN
          TTMP = 0.5*TREG
      ELSE
          TTMP = TREG
      END IF
*
*       Set new regular step and reduce irregular step if STEPR < STEP.
*     PRINT*,' New Step = ',TTMP,' Old ',TREG,' Quot ',TTMP/TREG
      TREG = TTMP
*
  110 IF (TTMP.LT.TIRR) THEN
          TIRR = 0.5D0*TIRR
          TIMENW(I) = T0(I) + TIRR
          NICONV = NICONV + 1
          GO TO 110
      END IF
*
      STEPR(I) = TREG
      STEP(I) = TIRR
      RS(I) = RSPH
*       Include any genuine high velocity particle in special list.
      IF (KZ(18).EQ.0.OR.TREG.GT.20.0*DTMIN) GO TO 120
*
      A1 = X0DOT(1,I)**2 + X0DOT(2,I)**2 + X0DOT(3,I)**2
      IF (A1.LT.16.0*ECLOSE.OR.TIRR.LT.DTMIN) GO TO 120
*
*       Check limit and see whether body #I has already been detected.
      NNB1 = LISTV(1) + 1
      IF (NNB1.GE.MLV) GO TO 120
      DO 115 L = 2,NNB1
          IF (I.EQ.LISTV(L)) GO TO 120
  115 CONTINUE
*
*       Add body #I to LISTV and increase membership.
      LISTV(NNB1+1) = I
      LISTV(1) = LISTV(1) + 1
      NFAST = NFAST + 1
      NCOUNT(30) = NCOUNT(30) + 1
*
  120 RETURN
*
      END
