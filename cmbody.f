      SUBROUTINE CMBODY(ENERGY,NSYS)
*
*
*       Formation of c.m. body by collision.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3
      COMMON/EBSAVE/  EBS
      REAL*8  CM(6)
      REAL  LUMS(6),TSCLS(6),ZM0,TM,TN
      CHARACTER*8  WHICH1
      DATA  AU  /2.0627E+05/
*
*
*       Set phase indicator to denote the case of collision.
      IPHASE = 9
*
*       Specify global indices of subsystem (membership: NSYS = 2, 3, 4).
      IF (NSYS.EQ.2) THEN
*       Save binding energy & separation and terminate KS pair.
          EB = BODY(2*KSPAIR-1)*BODY(2*KSPAIR)*H(KSPAIR)/BODY(N+KSPAIR)
          RB = R(KSPAIR)
          I = N + KSPAIR
*
*       Check for hierarchical configuration.
          NP1 = LIST(1,2*KSPAIR-1) + 1
          DO 5 L = 2,NP1
              J = LIST(L,2*KSPAIR-1)
              RIJ2 = 0.0
              VIJ2 = 0.0
              RDOT = 0.0
              DO 2 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
                  VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
                  RDOT = (X(K,I) - X(K,J))*(XDOT(K,I) - XDOT(K,J))
    2         CONTINUE
              RIP = SQRT(RIJ2)
              A1 = 2.0/RIP - VIJ2/(BODY(I) + BODY(J))
              A1 = 1.0/A1
              IF (1.0/A1.GT.0.5/RMIN) THEN
                  ECC2 = (1.0 - RIP/A1)**2 +
     &                                  RDOT**2/(A1*(BODY(I) + BODY(J)))
                  RP = A1*(1.0 - SQRT(ECC2))
                  A0 = -0.5*BODY(I)/H(KSPAIR)
                  ECC = 1.0 - R(KSPAIR)/A0
                  RA = A0*(1.0 + ECC)
                  SR = RP/RA
                  WRITE (6,4)  KSPAIR, H(KSPAIR), A0, A1, RP,
     &                         SQRT(ECC2), SR
    4             FORMAT (' HIERARCHY:   IPAIR H A0 A1 RP E1 SR ',
     &                                      I4,F7.0,1P3E9.1,0PF6.2,F6.1)
              END IF
    5     CONTINUE
*
          CALL KSTERM
          I1 = 2*NPAIRS + 1
          I2 = I1 + 1
          I3 = 0
          ICOMP = I1
          RCOLL = RB
          WHICH1 = ' BINARY '
      ELSE
          I1 = JLIST(1)
          I2 = JLIST(2)
          I3 = JLIST(3)
          I4 = JLIST(4)
*       Ignore case of three-body system here (JLIST(4) = 0).
          KSPAIR = NPAIRS
      END IF
*
*       Form global c.m. coordinates & velocities from body #I1 & I2.
      ZM = BODY(I1) + BODY(I2)
      DO 10 K = 1,3
          CM(K) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/ZM
          CM(K+3) = (BODY(I1)*XDOT(K,I1) + BODY(I2)*XDOT(K,I2))/ZM
   10 CONTINUE
*
*       Change stellar evolution parameters if mass loss option is active.
      IF (KZ(19).GE.3) THEN
          LI = NAME(I1)
          LJ = NAME(I2)
*
*       Update the original mass just in case.
          IF (LI.LE.500) THEN
              BODY0(LI) = ZM
          END IF
*
*       Define Thorne-Zytkow object or blue straggler from evolution type.
          IF (MAX(KSTAR(LI),KSTAR(LJ)).GE.5) THEN
              NTZ = NTZ + 1
              KSTAR(LI) = 9
*       Adopt mass loss rate of 1 Sun/Myr for 10 % loss and radius of 10 Sun.
              TEV(LI) = TIME + 0.1*ZM*ZMBAR/TSTAR
              RADIUS(I1) = 0.05*(RBAR/AU)
          ELSE
              NBS = NBS + 1
              ZM0 = ZM*ZMBAR
              CALL STAR(ZM0,TM,TN,TSCLS,LUMS)
*        Specify main sequence evolution time and modify radius & type.
              TEV(LI) = TIME + TM/TSTAR
              RADIUS(I1) = 1.5*MAX(RADIUS(I1),RADIUS(I2))
              KSTAR(LI) = MAX(KSTAR(LI),KSTAR(LJ))
          END IF
*
          WRITE (6,12)  NAME(I1), KSTAR(LI), ZM*ZMBAR, TEV(LI)*TSTAR,
     &                  (RADIUS(I1)/0.005)*(AU/RBAR)
   12     FORMAT (' NEW STAR    NAME KW M TEV R ',I5,I3,F6.1,F7.1,F6.1)
          TEV(LJ) = 1.0E+10
      END IF
*
*       Create new body from c.m. and initialize zero mass ghost in #I2.
      BODY(I1) = ZM
      BODY(I2) = 0.0D0
      LIST(1,I2) = 0
      T0(I2) = TADJ + DTADJ
      STEP(I2) = 1.0D+06
      RI = SQRT(X(1,I2)**2 + X(2,I2)**2 + X(3,I2)**2)
      VI = SQRT(XDOT(1,I2)**2 + XDOT(2,I2)**2 + XDOT(3,I2)**2)
      NAME1 = NAME(I1)
      NAME2 = NAME(I2)
*
      DO 20 K = 1,3
          X(K,I1) = CM(K)
          XDOT(K,I1) = CM(K+3)
          X0DOT(K,I1) = CM(K+3)
*       Ensure that ghost will escape next output (far from fast escapers).
          X0(K,I2) = 1000.0*RSCALE*X(K,I2)/RI
          X(K,I2) = X0(K,I2)
          X0DOT(K,I2) = SQRT(0.004*ZMASS/RSCALE)*XDOT(K,I2)/VI
          XDOT(K,I2) = X0DOT(K,I2)
          F(K,I2) = 0.0D0
          FDOT(K,I2) = 0.0D0
          D2(K,I2) = 0.0D0
          D3(K,I2) = 0.0D0
          D2R(K,I2) = 0.0D0
          D3R(K,I2) = 0.0D0
   20 CONTINUE
*
*       Copy members of neighbour or perturber list.
      NNB = LIST(1,I1)
      DO 30 L = 1,NNB
          JPERT(L) = LIST(L+1,I1)
   30 CONTINUE
*
*       Remove the ghost particle from neighbour lists containing #I1.
      JLIST(1) = I2
      CALL NBREM(I1,1,NNB)
*
*       Also remove ghost from list of I1 (use NTOT as dummy here).
      JPERT(1) = I1
      CALL NBREM(NTOT,1,1)
*
*       Replace body #I3 by spurious member (> 0) if it is an only neighbour.
      IF (LIST(1,I1).EQ.1.AND.LIST(2,I1).EQ.I3) THEN
          LIST(2,I1) = MIN(I1+I2+I3,N)
      END IF
*
*       Decide appropriate path for each case.
      IF (NSYS.EQ.2) GO TO 40
      IF (NSYS.EQ.3) GO TO 45
*
*       Switch KS components if body #I3 & I4 is closer than #I1 & I3.
      IF (JLIST(5).LT.0) THEN
          I4 = I1
          I1 = JLIST(4)
      END IF
*
*       Obtain dominant F & FDOT on body #I1 & I3 for #I4 in FPOLY2.
      JLIST(1) = I1
      JLIST(2) = I3
      JLIST(3) = I4
      ICOMP = I4
*
      CALL FCLOSE(I1,3)
      CALL FCLOSE(I3,3)
*
*       Predict neighbour coordinates & velocities in case of KS collision.
   40 IF (NSYS.EQ.2) THEN
          NNB = LIST(1,ICOMP)
          CALL XVPRED(ICOMP,NNB)
      END IF
*
*       Initialize force polynomials for new single or third body (ICOMP).
      CALL FPOLY1(ICOMP,ICOMP,0)
      CALL FPOLY2(ICOMP,ICOMP,0)
*
*       See whether body #ICOMP should be added to NLIST.
      IF (T0(ICOMP) + STEP(ICOMP).LT.TLIST) THEN
          CALL NLMOD(ICOMP,1)
      END IF
      IF (NSYS.EQ.2) GO TO 80
*
*       Obtain binding energy of the subsystem (2 or 3 members).
   45 ZKE = 0.0D0
      POTS = 0.0D0
      I = I1
      J = I3
   50 ZKE = ZKE + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
*       Note that variable ENERGY also contains c.m. kinetic energy.
      IF (I.EQ.I3.AND.NSYS.EQ.3) GO TO 60
      IF (I.EQ.I4) GO TO 60
   55 RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                              (X(3,I) - X(3,J))**2
      POTS = POTS + BODY(I)*BODY(J)/SQRT(RIJ2)
*
*       Include all interactions.
      IF (I.EQ.I1.AND.NSYS.EQ.4) THEN
          IF (J.EQ.I3) THEN
              J = I4
              GO TO 55
          ELSE
              I = I3
              J = I4
              GO TO 50
          END IF
      END IF
*
*       Add kinetic energy from last body and check DMIN (only #I3 or I4).
      IF (NSYS.EQ.3) THEN
          I = I3
          WHICH1 = ' TRIPLE '
          DMIN3 = MIN(DMIN3,RCOLL)
      ELSE
          I = I4
          WHICH1 = '   QUAD '
          DMIN4 = MIN(DMIN4,RCOLL)
      END IF
*
      GO TO 50
*
*       Form net energy correction for triple or quad case.
   60 EB = ENERGY - (0.5D0*ZKE - POTS)
      RB = SQRT(RIJ2)
*
*       Set global components for new KS regularization (ICOMP < JCOMP).
      ICOMP = MIN(I1,I3)
      JCOMP = MAX(I1,I3)
*
*       Initialize new KS pair.
      CALL KSREG
*
*       Update energy loss & collision counters.
   80 ECOLL = ECOLL + EB
      E(10) = E(10) + EB
      NPOP(8) = NPOP(8) + 1
      NCOUNT(28) = NCOUNT(28) + 1
      NCOLL = NCOLL + 1
*
      WRITE (6,90)  WHICH1, KSPAIR, NAME1, NAME2, ZM, RCOLL, RB, EB,
     &              ECOLL
   90 FORMAT (/,A8,'COLLISION    KSPAIR =',I3,'  NAME =',2I6,
     &             '  M =',F7.4,'  RCOLL =',1PE8.1,'  R =',E8.1,
     &             '  EB =',0PF9.5,'  ECOLL =',F10.6)
*
      IF (NSYS.EQ.4) WRITE (6,95)  EB,EBS,(EB-EBS)/EB
   95 FORMAT (' CMBODY:  EB EBS DE/E  ',1P,3E10.2)
      RETURN
*
      END
