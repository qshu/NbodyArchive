      SUBROUTINE CHTERM(ISUB)
*
*
*       Termination of chain system.
*       ----------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,R2(NMX,NMX)
      INTEGER  IJ(NMX)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX,1)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
      COMMON/ECHAIN/  ECH
*
*
*       Decide between standard termination or collision (ISUB > 0 or < 0).
      IF (ISUB.NE.0) THEN
          ITERM = ISUB
          ISUB = IABS(ISUB)
      END IF
*
*       Prepare KS regularization and direct integration of two bodies.
      CALL R2SORT(IJ,R2)
      I1 = IJ(1)
      I2 = IJ(2)
      I3 = IJ(3)
      IF (NCH.EQ.2) I3 = I2
      IF (NCH.LE.3) THEN
          I4 = I1
          R2(I2,I4) = R2(I2,I3)
      ELSE
          I4 = IJ(4)
      END IF
*
      IF (KZ(30).GT.2) THEN
          WRITE (6,1)  SQRT(R2(I1,I2)), SQRT(R2(I1,I3)),SQRT(R2(I2,I3)),
     &                 SQRT(R2(I2,I4)), SQRT(R2(I3,I4))
    1     FORMAT (' CHTERM:   RIJ (1-2 1-3 2-3 2-4 3-4)  ',1P,5E9.1)
      END IF
*
      JLIST(5) = NAMEC(I1)
      JLIST(6) = NAMEC(I2)
      JLIST(7) = NAMEC(I3)
      JLIST(8) = NAMEC(I4)
*
*       Specify chain phase indicator and restore original name to c.m. body.
      IPHASE = 8
      NAME(ICH) = NAME0
*
*       Indentify current global indices by searching all single particles.
      DO 102 J = IFIRST,N
          DO 101 L = 1,NCH
              IF (NAME(J).EQ.JLIST(L+4)) THEN
                  JLIST(L) = J
                  IF (BODY(J).GT.0.0D0) ICM = J
              END IF
  101     CONTINUE
  102 CONTINUE
*
*       Modify identification list for special cases NCH = 2 & NCH = 3.
      IF (NCH.EQ.2) THEN
          I3 = I1
          I4 = I2
          JLIST(3) = JLIST(1)
          JLIST(4) = JLIST(2)
      ELSE IF (NCH.EQ.3) THEN
          JLIST(4) = JLIST(3)
      END IF
*
*       Ensure ICOMP < JCOMP for KS regularization.
      ICOMP = MIN(JLIST(1),JLIST(2))
      JCOMP = MAX(JLIST(1),JLIST(2))
*
*       Copy final coordinates & velocities to standard variables.
      LK = 0
      DO 104 L = 1,NCH
          DO 103 K = 1,3
              LK = LK + 1
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
  103     CONTINUE
  104 CONTINUE
*
*       Define discrete time for new polynomials (subject to TIME <= TBLOCK).
      TIME2 = T0S(ISUB) + TIMEC - TPREV
      DT = 0.1*STEP(ICM)
      CALL STEPK(DT,DTN)
      TIME = TPREV + INT((TIME2 + DT)/DTN)*DTN
      TIME = MIN(TBLOCK,TIME)
*
*       Predict current coordinates & velocities to F3DOT before termination.
      CALL XVPRED(ICM,-1)
*
*       Copy c.m. coordinates & velocities.
      DO 105 K = 1,3
          CM(K) = X(K,ICM)
          CM(K+3) = XDOT(K,ICM)
  105 CONTINUE
*
*       Set configuration pointers for KS candidates & distant bodies.
      JLIST(5) = I1
      JLIST(6) = I2
      JLIST(7) = I3
      JLIST(8) = I4
*
*       Place new coordinates in the original locations.
      DO 120 L = 1,NCH
          J = JLIST(L)
*       Compare global name & subsystem name to restore the mass & T0.
          DO 112 K = 1,NCH
              IF (NAME(J).EQ.NAMEC(K)) THEN
                  BODY(J) = BODYC(K)
                  T0(J) = TIME
              END IF
  112     CONTINUE
*       Transform to global coordinates & velocities using c.m. values.
          LL = JLIST(L+4)
          DO 115 K = 1,3
              X(K,J) = X4(K,LL) + CM(K)
              XDOT(K,J) = XDOT4(K,LL) + CM(K+3)
              X0(K,J) = X(K,J)
              X0DOT(K,J) = XDOT(K,J)
  115     CONTINUE
  120 CONTINUE
*
*       Set one ghost index for skipping neighbour list search.
      JG = ICOMP
      IF (JG.EQ.ICM) JG = JCOMP
*
*       Search all neighbour lists for splitting chain c.m. into components.
      NNB1 = 0
      DO 130 J = IFIRST,NTOT
          NNB2 = LIST(1,J) + 1
          DO 125 L = 2,NNB2
              IF (LIST(L,J).GT.ICM) GO TO 130
*       Form list of all ICM without one ghost and predict X & XDOT to FDOT.
              IF (LIST(L,J).EQ.ICM) THEN
                  DO 122 K = 2,NNB2
                      IF (LIST(K,J).EQ.JG) GO TO 130
  122             CONTINUE
*       Skip subsystem members.
                  DO 124 K = 1,NCH
                      IF (J.EQ.JLIST(K)) GO TO 130
  124             CONTINUE
                  IF (NNB1.LT.LMAX) THEN
                      NNB1 = NNB1 + 1
                      JPERT(NNB1) = J
                      CALL XVPRED(J,0)
                  END IF
              END IF
  125     CONTINUE
  130 CONTINUE
*
*       Update subsystem COMMON variables unless last or only case.
      IF (ISUB.LT.NSUB) THEN
          DO 150 L = ISUB,NSUB
              DO 145 K = 1,6
                  BODYS(K,L) = BODYS(K,L+1)
                  NAMES(K,L) = NAMES(K,L+1)
  145         CONTINUE
              T0S(L) = T0S(L+1)
              TS(L) = TS(L+1)
              STEPS(L) = STEPS(L+1)
              RMAXS(L) = RMAXS(L+1)
              ISYS(L) = ISYS(L+1)
  150     CONTINUE
      END IF
*
*       Reduce subsystem counter and initialize membership & internal energy.
      NSUB = NSUB - 1
      NCH0 = NCH
      NCH = 0
      ECH = 0.0
*     ESUB = ESUB - ENERGY - ECOLL3
*
*       Assign new neighbours for dominant KS, I3 & I4 (Note: JLIST is *2).
      RS0 = RS(ICM)
      J3 = JLIST(3)
      J4 = JLIST(4)
      CALL NBLIST(ICOMP,RS0)
      IF (NCH0.GT.2) CALL NBLIST(J3,RS0)
      IF (NCH0.GT.3) CALL NBLIST(J4,RS0)
*
*       Check for stellar collision (only needs coordinates & velocities).
      IF (ITERM.LT.0) THEN
          JLIST(1) = ICOMP
          JLIST(2) = JCOMP
*
*       See whether relabelling is required (indices I1 - I4 still local).
          IF (R2(I1,I4).LT.R2(I1,I3).OR.R2(I3,I4).LT.R2(I1,I3)) THEN
              IF (R2(I1,I4).LT.R2(I3,I4)) THEN
*       Switch body #I3 & I4 to give new dominant pair I1 & I3.
                  I = JLIST(4)
                  JLIST(4) = JLIST(3)
                  JLIST(3) = I
              ELSE
*       Set JLIST(5) < 0 to denote that body #I3 & I4 will be new KS pair.
                  JLIST(5) = -1
              END IF
          END IF
          GO TO 170
      END IF
*
*       Replace ICM in perturber neighbour lists by all subsystem members.
      CALL NBREST(ICM,NCH0,NNB1)
*
*       Exclude the dominant interaction for c.m. approximation (large FDOT).
      IF (MIN(R2(I1,I3),R2(I2,I4)).GT.CMSEP2*R2(I1,I2)) THEN
          JLIST(1) = JLIST(3)
          JLIST(2) = JLIST(4)
          NNB = 2
          IF (NCH0.EQ.3) JLIST(2) = JCLOSE
      ELSE
          NNB = NCH0
      END IF
*
      IF (NCH0.EQ.2) THEN
          JLIST(1) = JCLOSE
          NNB = 1
      END IF
*
*       Set dominant F & FDOT on body #ICOMP & JCOMP for #I3 & I4 in FPOLY2.
      CALL FCLOSE(ICOMP,NNB)
      CALL FCLOSE(JCOMP,NNB)
*
*       See whether a second binary is present (NCH0 = 4 & RB < RMIN).
      KS2 = 0
      IF (NCH0.EQ.4.AND.R2(I3,I4).LT.RMIN2) THEN
*       Accept second KS pair if well separated (> 2) from smallest binary.
          IF (R2(I3,I4).LT.0.25*MIN(R2(I1,I3),R2(I2,I4))) THEN
              KS2 = 1
          END IF
      END IF
*
*       Specify global indices of least dominant bodies.
      I3 = JLIST(3)
      I4 = JLIST(4)
*
*       Set global names and dominant F & FDOT of #I3 & I4 for new KS pair.
      IF (KS2.GT.0) THEN
          NAME3 = NAME(I3)
          NAME4 = NAME(I4)
          JLIST(1) = ICOMP
          JLIST(2) = JCOMP
          CALL FCLOSE(I3,4)
          CALL FCLOSE(I4,4)
      END IF
*
*       Initialize force polynomials & time-steps for body #I3 & #I4.
      IF (NCH0.GT.2.AND.KS2.EQ.0) THEN
          CALL FPOLY1(I3,I3,0)
          IF (NCH0.GE.4) THEN
              CALL FPOLY1(I4,I4,0)
              CALL FPOLY2(I4,I4,0)
          END IF
          CALL FPOLY2(I3,I3,0)
      END IF
*
*       Perform KS regularization of dominant components (ICOMP < JCOMP).
      STEP3 = STEP(I3)
      STEP4 = STEP(I4)
      CALL KSREG
*
      IF (KZ(30).GT.1) THEN
          WRITE (6,155)  TIME, I3, I4, NTOT, STEP3, STEP4, STEP(NTOT)
  155     FORMAT (' CHTERM:   T I3 I4 NT DT ',F9.4,3I5,1P,3E9.1)
      END IF
*
*       Initialize second KS pair if separation is small.
      IF (KS2.GT.0) THEN
          I3 = 0
*       Re-determine the global indices after updating in KSREG.
          DO 164 I = IFIRST,N
              IF (NAME(I).EQ.NAME3) I3 = I
              IF (NAME(I).EQ.NAME4) THEN
                  I4 = I
                  IF (I3.GT.0) GO TO 165
              END IF
  164     CONTINUE
*
*       Define components and perform new regularization.
  165     ICOMP = MIN(I3,I4)
          JCOMP = MAX(I3,I4)
          CALL KSREG
*
          IF (KZ(30).GT.1) THEN
              WRITE (6,166)  I3, I4, STEP(NTOT), STEPR(NTOT),
     &                       R(NPAIRS), H(NPAIRS), GAMMA(NPAIRS)
  166         FORMAT (' CHTERM:   SECOND BINARY   I3 I4 DT DTR R H G ',
     &                                            2I5,1P,5E9.1)
          END IF
      END IF
*
*       Check minimum two-body distance (NB! RCOLL not determined yet).
*     DMINC = MIN(DMINC,RCOLL)
*
*       Update net binary energy change.
      CHCOLL = CHCOLL + CM(9)
*
*       Update number of DIFSY calls, tidal dissipations & collision energy.
  170 NSTEPC = NSTEPC + NSTEP1
*     NDISS = NDISS + NDISS4
*     ECOLL = ECOLL + ECOLL3
*
*       Check for subsystem at last COMMON dump (no restart with NSUB > 0).
      IF (NSUB.EQ.0.AND.KZ(2).GE.1) THEN
          IF (TIME - TDUMP.LT.TIMEC) THEN
              TDUMP = TIME
              CALL MYDUMP(1,2)
          END IF
      END IF
*
*       Set phase indicator < 0 to ensure new NLIST in routine INTGRT.
      IPHASE = -1
*
      RETURN
*
      END
