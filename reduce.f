      SUBROUTINE REDUCE(IESC,JESC,ISUB)
*
*
*       Reduction of chain.
*       -------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,XCM(3),VCM(3),DXC(3),DVC(3),CG(6)
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
      REAL*8  FIRR(3),FD(3)
*
*
*       Define discrete time for new polynomial.
      TIME0 = TIME
      TIME2 = T0S(ISUB) + TIMEC - TPREV
      DT = 0.01*STEP(ICH)
      CALL STEPK(DT,DTN)
      TIME = TPREV + INT((TIME2 + DT)/DTN)*DTN
      TIME = MIN(TBLOCK,TIME)
*       Re-define initial epoch for consistency (ignore phase error).
      T0S(ISUB) = TIME - TIMEC
      IF (KZ(30).GT.2) THEN
          WRITE (6,1)  TIME0, TIME2, TIME, DTN
    1     FORMAT (' REDUCE:   TIME0 TIME2 TIME DTN ',3F10.6,1P,E9.1)
      END IF
*
*       Make new chain from remaining members.
    2 LK = 0
      DO 10 L = 1,NCH
          IF (L.EQ.IESC) GO TO 10
          DO 5 K = 1,3
              LK = LK + 1
              XCH(LK) = X4(K,L)
              VCH(LK) = XDOT4(K,L)
    5     CONTINUE
   10 CONTINUE
*
*       Reduce chain membership and mass (global & local COMMON).
      NCH = NCH - 1
      NN = NCH
      MASS = MASS - M(IESC)
*
*       Improve coordinates & velocities of c.m. body to order F3DOT.
      CALL XVPRED(ICH,-1)
*
*       Set new c.m. for reduced system and save old c.m. variables.
      DO 20 K = 1,3
          DXC(K) = -BODYC(IESC)*X4(K,IESC)/(BODY(ICH) - BODYC(IESC))
          DVC(K) = -BODYC(IESC)*XDOT4(K,IESC)/(BODY(ICH) - BODYC(IESC))
          XCM(K) = X(K,ICH) + DXC(K)
          VCM(K) = XDOT(K,ICH) + DVC(K)
          CM(K) = X(K,ICH)
          CM(K+3) = XDOT(K,ICH)
   20 CONTINUE
*
*       Re-define new chain variables w.r. to modified c.m. (NB! retain X4).
      LK = 0
      DO 30 L = 1,NCH
          DO 25 K = 1,3
              LK = LK + 1
              XCH(LK) = XCH(LK) - DXC(K)
              VCH(LK) = VCH(LK) - DVC(K)
   25     CONTINUE
   30 CONTINUE
*
*       Save original mass & neighbour radius of c.m. body.
      BODYCH = BODY(ICH)
      RS0 = RS(ICH)
*
*       Search for global index of escaper (assume single body for now).
      DO 40 J = IFIRST,N
          IF (NAME(J).EQ.NAMEC(IESC)) THEN
              I = J
              IF (BODY(J).GT.0.0D0) WRITE (6,35)  I, IESC, NAMEC(IESC)
   35         FORMAT (' WARNING!   NON-ZERO GHOST    I IESC NAMEC ',3I5)
              GO TO 55
          END IF
   40 CONTINUE
*
*       Switch to another reference body since body #ICH is escaping.
      I = ICH
      IF (IESC.GT.1) THEN
          NEW = 1
      ELSE
          NEW = 2
      END IF
*
*       Identify global index of new reference body.
      DO 45 J = IFIRST,N
          IF (NAME(J).EQ.NAMEC(NEW)) THEN
              ICH = J
              GO TO 50
          END IF
   45 CONTINUE
*
*       Include warning if no reference body (this should not occur).
      WRITE (6,48)  IESC, NAMEC(NEW)
   48 FORMAT (' REDUCE:   DANGER!   NO REFERENCE BODY    IESC NAME',2I5)
      NCH = NCH + 1
      NN = NCH
      MASS = MASS + M(IESC)
      GO TO 100
*
*       Copy neighbour list of old c.m. body for routine NBREST.
   50 NNB = LIST(1,I)
      DO 52 L = 2,NNB+1
          JPERT(L-1) = LIST(L,I)
   52 CONTINUE
*
*       Restore ghost to lists of neighbours (body #I will be skipped).
      JLIST(1) = ICH
      CALL NBREST(I,1,NNB)
*
      WRITE (6,53)  NAME(I), NAME0, NAME(ICH), ICH
   53 FORMAT (' SWITCH REFERENCE    NAMEI NAME0 NAMECH ICH ',4I5)
*
*       Exchange name of reference body and initialize new c.m. name.
      NAME(I) = NAME0
      NAME0 = NAME(ICH)
      NAME(ICH) = 0
*
*       Update total mass and initialize new c.m. body variables.
   55 BODY(ICH) = BODYCH - BODYC(IESC)
      CM(7) = BODY(ICH)
      T0(ICH) = TIME
      DO 60 K = 1,3
          X(K,ICH) = XCM(K)
          X0(K,ICH) = XCM(K)
          XDOT(K,ICH) = VCM(K)
          X0DOT(K,ICH) = VCM(K)
   60 CONTINUE
*
*       Restore the mass and transform to global coordinates & velocities.
      BODY(I) = BODYC(IESC)
      T0(I) = TIME
      DO 65 K = 1,3
          X(K,I) = X4(K,IESC) + CM(K)
          XDOT(K,I) = XDOT4(K,IESC) + CM(K+3)
          X0(K,I) = X(K,I)
          X0DOT(K,I) = XDOT(K,I)
   65 CONTINUE
*
*       Remove chain (and clump) mass & reference name of escaper.
      DO 70 L = IESC,NCH
          M(L) = M(L+1)
          BODYC(L) = BODYC(L+1)
          NAMEC(L) = NAMEC(L+1)
          BODYS(L,ISUB) = BODYS(L+1,ISUB)
          NAMES(L,ISUB) = NAMES(L+1,ISUB)
   70 CONTINUE
*
*       Initialize neighbour list using current radius of (old) #ICH.
      CALL NBLIST(ICH,RS0)
*
*       Set dominant F & FDOT on body #I for #ICH in FPOLY2.
*     JLIST(1) = ICH
*     JLIST(2) = JCLOSE
*     CALL FCLOSE(I,2)
*
*       Perform re-initialization of c.m. polynomials & perturber list.
      CALL REINIT(ISUB)
*
*       Copy new clump coordinates & velocities to standard variables.
      LK = 0
      DO 80 L = 1,NCH
          DO 75 K = 1,3
          LK = LK + 1
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
   75     CONTINUE
   80 CONTINUE
*
*       Copy neighbour list of cm. body for routine NBREST.
      NNB = LIST(1,ICH)
      DO 85 L = 2,NNB+1
          JPERT(L-1) = LIST(L,ICH)
   85 CONTINUE
*
*       Restore ghost to lists of neighbours (body #ICH will be skipped).
      JLIST(1) = I
      CALL NBREST(ICH,1,NNB)
*
*       Make new neighbour list.
      CALL NBLIST(I,RS0)
*
*       Distinguish between single particle and binary (JESC = 0 & > 0).
      IF (JESC.EQ.0) THEN
*       Initialize force polynomials & time-steps (add differential F & FD).
          IPHASE = 8
          CALL FPOLY1(I,I,0)
          DO 82 K = 1,3
              FIRR(K) = 0.0
              FD(K) = 0.0
   82     CONTINUE
          CALL FCHAIN(I,0,X(1,I),XDOT(1,I),FIRR,FD)
          DO 84 K = 1,3
              F(K,I) = F(K,I) + 0.5*FIRR(K)
              FDOT(K,I) = FDOT(K,I) + ONE6*FD(K)
   84     CONTINUE
          CALL FPOLY2(I,I,0)
          IPHASE = -1
*       Include case of escaping binary (set JESC < 0 for new KS).
      ELSE IF (JESC.GT.0) THEN
          ICLOSE = I
          IF (JESC.GT.IESC) JESC = JESC - 1
          IESC = JESC
          JESC = -1
          GO TO 2
      ELSE
*       Initialize KS regularization after second reduction.
          ICOMP = MIN(ICLOSE,I)
          JCOMP = MAX(ICLOSE,I)
          CALL KSREG
      END IF
*
*       Prepare c.m. check.
      DO 88 K = 1,6
          CG(K) = 0.0
   88 CONTINUE
*
      LK = 0
      DO 95 L = 1,NCH
          DO 90 K = 1,3
          LK = LK + 1
              CG(K) = CG(K) + BODYC(L)*XCH(LK)
              CG(K+3) = CG(K+3) + BODYC(L)*VCH(LK)
   90     CONTINUE
   95 CONTINUE
*
      DO 96 K = 1,6
          CG(K) = CG(K)/BODY(ICH)
   96 CONTINUE
*
*     WRITE (6,97)  TIME, (CG(K),K=1,6)
*  97 FORMAT (' REDUCE:   T CG ',F10.5,1P,6E9.1)
*
*       Ensure current coordinates & velocities for chain components.
      CALL XCPRED(1)
*
  100 RETURN
*
      END
