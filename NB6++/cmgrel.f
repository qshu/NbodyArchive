      SUBROUTINE CMGREL(ENERGY,NSYS)
*
*
*       Formation of c.m. body by collision.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX4=4*NMX)
      COMMON/CLOSE/  RIJ4(4,4),RCOLL4,QPERI4,SIZE4(4),ECOLL4,IP(4)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),ASYNC,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/EBSAVE/  EBS
      REAL*8  CM(6),LUMS(10),TSCLS(20),GB(10)
      REAL*8  M01,M02,M03,M1,M2,M3,LUM,MC1,MC2,MC3
C      SAVE KTYPE
C      INTEGER  KTYPE(0:14,0:14)
      DATA KTYPE /1,1,1,23,24,25,25,24,24,24,313,315,315,5,5,
     &            1,1,1,23,24,25,25,24,24,24,313,315,315,5,5,
     &            1,1,1,23,24,25,25,24,24,24,313,315,315,5,5,
     &            23,23,23,23,24,25,25,23,23,23,313,315,315,5,5,
     &            24,24,24,24,24,25,25,24,24,24,313,315,315,5,5,
     &            25,25,25,25,25,25,25,25,25,25,313,315,315,5,5,
     &            25,25,25,25,25,25,25,25,25,25,313,315,315,5,5,
     &            24,24,24,23,24,25,25,4,4,4,4,4,4,5,5,
     &            24,24,24,23,24,25,25,4,4,4,4,4,4,5,5,
     &            24,24,24,23,24,25,25,4,4,4,4,4,4,5,5,
     &            323,323,323,323,323,323,323,4,4,4,5,5,5,5,5,
     &            325,325,325,325,325,325,325,4,4,4,5,5,5,5,5,
     &            325,325,325,325,325,325,325,4,4,4,5,5,5,5,5,
     &            5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
     &            5,5,5,5,5,5,5,5,5,5,5,5,5,5,5/
      CHARACTER*8  WHICH1
*
*
*       Distinguish between chain and triple or quad case (ICH > 0 or = 0).
      IF (IPHASE.EQ.9) THEN
          ICH = 1
      ELSE
*       Activate collision indicator (otherwise done in CHTERM).
          ICH = 0
          IPHASE = 9
      END IF
*
*       Specify global indices of subsystem (membership: NSYS = 2 - 5).
      IF (NSYS.EQ.2) THEN
*
*       Define discrete time for prediction & new polynomials (T <= TBLOCK).
          DT2 = TIME - TPREV
          DT8 = (TBLOCK - TPREV)/8.0D0
          TIME = TPREV + NINT(DT2/DT8)*DT8
          TIME = MIN(TBLOCK,TIME)
          TIME0 = TIME
*
*       Save binding energy and terminate KS pair.
          EB = BODY(2*KSPAIR-1)*BODY(2*KSPAIR)*H(KSPAIR)/BODY(N+KSPAIR)
          I = N + KSPAIR
          JCLOSE = 0
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
                  if(rank.eq.0)
     &            WRITE (6,4)  KSPAIR, NAME(J), H(KSPAIR), ECC, A0, A1,
     &                         RP, ECC, SR
    4             FORMAT (' HIERARCHY:   KS NMJ H E A0 A1 RP E1 SR',
     &                              2I6,F7.0,F8.4,1P,3E9.1,0PF6.2,F6.1)
              END IF
*       Select closest single body inside 0.5*RMIN as KS component.
              IF (RIP.LT.0.5*RMIN.AND.J.LE.N) THEN
                  IF (JCLOSE.GT.0) THEN
                      IF (RIP.GT.RIP0) GO TO 5
                      JCLOSE = J
                      RIP0 = RIP
                  ELSE
                      JCLOSE = J
                      RIP0 = RIP
                  END IF
              END IF
    5     CONTINUE
*
*       Search for evidence of recent regularization.
          NAM1 = NAME(2*KSPAIR-1)
          NAM2 = NAME(2*KSPAIR)
          NNB = LISTD(1)
          DO 7 K = 2,NNB+1
              IF (LISTD(K).EQ.NAM1.OR.LISTD(K).EQ.NAM2) THEN
                  if(rank.eq.0)
     &            WRITE (6,6)  NAM1, NAM2, LISTD(K), K
    6             FORMAT (' KS REMNANT:    NAM LISTD K  ',3I6,I4)
              END IF
    7     CONTINUE
*
*       Update body #JCLOSE to current time for new KS with combined c.m.
          IF (JCLOSE.GT.0) THEN
              CALL XVPRED(JCLOSE,-1)
              T0(JCLOSE) = TIME
              DO 8 K = 1,3
                  X0DOT(K,JCLOSE) = XDOT(K,JCLOSE)
                  X0(K,JCLOSE) = X(K,JCLOSE)
    8         CONTINUE
          END IF
*
*       Ensure orbit is at pericentre (perturbed hyperbolic case is OK).
          SEMI = -0.5*BODY(N+KSPAIR)/H(KSPAIR)
          IF (R(KSPAIR).GT.SEMI.AND.SEMI.GT.0.0) THEN
              CALL KSAPO(KSPAIR)
              CALL KSPERI(KSPAIR)
*       Restore quantized time to avoid small STEP (KSTERM needs T0 = TIME).
              TIME = TIME0
              T0(2*KSPAIR-1) = TIME
          ELSE IF (KSTAR(N+KSPAIR).GT.0) THEN
              CALL KSPERI(KSPAIR)
              TIME = TIME0
              T0(2*KSPAIR-1) = TIME
          END IF
*
*       Save collision distance and VINF (km/sec).
          RCOLL = R(KSPAIR)
          VINF = 0.0
          IF (H(KSPAIR).GT.0.0) VINF = SQRT(2.0*H(KSPAIR))*VSTAR
          ECC = 1.0 - R(KSPAIR)/SEMI
*
*       Terminate KS pair and set relevant indices for collision treatment.
          ICOLLISION = 1
          CALL KSTERM
          I1 = 2*NPAIRS + 1
          I2 = I1 + 1
          I3 = 0
          ICOMP = I1
          DMIN2 = MIN(DMIN2,RCOLL)
          WHICH1 = 'BH_BINARY'
      ELSE
          I1 = JLIST(1)
          I2 = JLIST(2)
          I3 = JLIST(3)
          I4 = JLIST(4)
          I5 = JLIST(5)
*       Ignore case of three-body system here (JLIST(4) = 0).
          VINF = 0.0
          ECC = 1.0 + 2.0*EBS*DMINC/(BODY(I1)*BODY(I2))
          ECC = MAX(ECC,0.0D0)
          IF (EBS.GT.0.0) THEN
              HI = EBS*(BODY(I1) + BODY(I2))/(BODY(I1)*BODY(I2))
              VINF = SQRT(2.0*HI)*VSTAR
          END IF
      END IF
*
*       Set new quantized time (note: restore if problems in CHTERM).
*     TIME = TBLOCK
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
          M1 = BODY(I1)*ZMBAR
          M2 = BODY(I2)*ZMBAR
          RS1 = RADIUS(I1)*SU
          RS2 = RADIUS(I2)*SU
          TPHYS = (TIME + TOFF)*TSTAR
*
*       Determine evolution time scales for first star.
          KW1 = KSTAR(I1)
          M01 = BODY0(I1)*ZMBAR
          AGE1 = TEV(I1)*TSTAR - EPOCH(I1)
          CALL STAR(KW1,M01,M1,TMS1,TN,TSCLS,LUMS,GB,ZPARS)
          CALL HRDIAG(M01,AGE1,M1,TMS1,TN,TSCLS,LUMS,GB,ZPARS,
     &                RM,LUM,KW1,MC1,RCC)
*       Obtain time scales for second star.
          KW2 = KSTAR(I2)
          M02 = BODY0(I2)*ZMBAR
          AGE2 = TEV(I2)*TSTAR - EPOCH(I2)
          CALL STAR(KW2,M02,M2,TMS2,TN,TSCLS,LUMS,GB,ZPARS)
          CALL HRDIAG(M02,AGE2,M2,TMS2,TN,TSCLS,LUMS,GB,ZPARS,
     &                RM,LUM,KW2,MC2,RCC)
*
*       Treat different types using collision matrix (Chris Tout 4/4/02).
          M3 = M1 + M2
          IF (KTYPE(KW1,KW2).EQ.1) THEN
*       Specify new age based on complete mixing (KSTAR <= 2).
              M03 = M01 + M02
              KW = 1
              CALL STAR(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
              AGE3 = 0.1*TMS3*(AGE1*M1/TMS1 + AGE2*M2/TMS2)/M3
*       Determine proper final type for giants.
          ELSE IF (KTYPE(KW1,KW2)/10.EQ.2) THEN
              KW = KTYPE(KW1,KW2) - 20
              MC3 = MAX(MC1,MC2)
*       Distinguish between first or second component being a white dwarf.
          ELSE IF (KTYPE(KW1,KW2)/10.EQ.31) THEN
              MC3 = M1
              KW = KTYPE(KW1,KW2) - 310
          ELSE IF (KTYPE(KW1,KW2)/10.EQ.32) THEN
              MC3 = M2
              KW = KTYPE(KW1,KW2) - 320
          ELSE
*       Define neutron star or black hole in other cases (depends on M3).
              M03 = M3
              KW = 13
              AGE3 = 0.0
          END IF
*
*       Obtain initial mass and age from the general evolution package.
          IF (KTYPE(KW1,KW2).GT.20) THEN
              AGE3 = 0.0
              CALL GNTAGE(MC3,M3,KW,ZPARS,M03,AGE3)
*       Note absence of mass loss (routine COMENV not included yet).
          END IF
*
*       Update type, initial mass, evolution times and epoch.
          KSTAR(I1) = KW
          BODY0(I1) = M03/ZMBAR
          TEV0(I1) = TIME
          TEV(I1) = TIME + 0.001/TSTAR
          EPOCH(I1) = TEV0(I1)*TSTAR - AGE3
*
*       Check for blue straggler.
          IF (TMS3.LT.TPHYS.AND.KW.LE.1) THEN
              NBS = NBS + 1
          END IF
*
          if(rank.eq.0)
     &    WRITE (6,15)  NAME(I1), KSTAR(I1), TEV(I1)*TSTAR, M1, M2,
     &                  RS1, RS2, TPHYS, AGE3
   15     FORMAT (' NEW STAR:    NAME KW TEV M1 M2 R1 R2 TP AGE ',
     &                           I6,I4,F7.1,2F6.1,F7.1,F6.1,2F8.1)
          TEV(I2) = 1.0E+10
      END IF
*
*       Copy all members of neighbour or perturber list.
      NNB = LIST(1,I1)
      DO 20 L = 1,NNB
          JPERT(L) = LIST(L+1,I1)
   20 CONTINUE
*
*       Evaluate potential energy with respect to colliding bodies.
      IF (NSYS.EQ.2) THEN
          JLIST(1) = I1
          JLIST(2) = I2
*       Replace second old KS component temporarily by arbitrary body.
          JPERT(1) = N
          IF (I2.EQ.N) JPERT(1) = N + 1
          CALL NBPOT(2,NNB,POT1)
      ELSE
*       Obtain differential effect on #I1 & #I2 due to closest member #I3.
          JPERT(1) = I3
          CALL NBPOT(2,1,POT1)
*       Include possible fourth member.
          IF (NSYS.GT.3) THEN
              JPERT(1) = I4
              CALL NBPOT(2,1,POT3)
              POT1 = POT1 + POT3
          END IF
      END IF
*
*       Create new body from c.m. and initialize zero mass ghost in #I2.
      BODY(I1) = ZM
      BODY(I2) = 0.0
      LIST(1,I2) = 0
      T0(I2) = TADJ + DTADJ 
      STEP(I2) = 1.0D+06
      STEPR(I2) = 1.0D+06
      RI = SQRT(X(1,I2)**2 + X(2,I2)**2 + X(3,I2)**2)
      VI = SQRT(XDOT(1,I2)**2 + XDOT(2,I2)**2 + XDOT(3,I2)**2)
      NAME1 = NAME(I1)
      NAME2 = NAME(I2)
*
      DO 30 K = 1,3
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
          D0(K,I2) = 0.0
          D1(K,I2) = 0.0
          D2(K,I2) = 0.0D0
          D3(K,I2) = 0.0D0
          D0R(K,I2) = 0.0
          D1R(K,I2) = 0.0
          D2R(K,I2) = 0.0D0
          D3R(K,I2) = 0.0D0
   30 CONTINUE
*
*       Obtain potential energy w.r.t. new c.m. and apply tidal correction.
      IF (NSYS.EQ.2) THEN
          CALL NBPOT(1,NNB,POT2)
      ELSE
          JPERT(1) = I3
          CALL NBPOT(1,1,POT2)
          IF (NSYS.GT.3) THEN
              JPERT(1) = I4
              CALL NBPOT(1,1,POT4)
              POT2 = POT2 + POT4
          END IF
      END IF
      DP = POT2 - POT1
      ECOLL = ECOLL + DP
*
*       ---------------------------
*       Additions for PN evolution.
*
*       Identify most massive member, merger product will inherit
*       this name.
      IF (BODY(I1).GT.BODY(I2)) THEN
         IMERGER = I1
         IGHOST = I2
      ELSE
         IMERGER = I2
         IGHOST = I1
      END IF
*
*       Now stustitute I1 with IMERGER and I2 with IGHOST in the rest
*       of the code.
*
*       Remove the ghost particle from neighbour lists containing #IMERGER.
      JPERT(1) = IGHOST
      JLIST(1) = IGHOST
      CALL NBREM(IMERGER,1,NNB)
*
*       Remove ghost from list of IMERGER (use NTOT as dummy here).
      JPERT(1) = IMERGER
      CALL NBREM(NTOT,1,1)
*
*       Replace body #I3 by spurious member (> 0) if it is an only neighbour.
      IF (LIST(1,IMERGER).EQ.1.AND.LIST(2,IMERGER).EQ.I3) THEN
          LIST(2,IMERGER) = MIN(IMERGER+IGHOST+I3,N)
      END IF
*
*       Decide appropriate path for each case.
      IF (NSYS.EQ.2) GO TO 40
      IF (NSYS.EQ.3) GO TO 45
*
*       Switch KS components if body #I3 & I4 is closer than #IMERGER & I3.
      IF (JLIST(6).LT.0) THEN
          I4 = IMERGER
          IMERGERS = IMERGER
          IMERGER = JLIST(4)
          JLIST(4) = IMERGERS
      END IF
*
*       Obtain dominant F & FDOT on body #IMERGER & I3 for #I4 in FPOLY2.
      ICOMP = I4
      CALL FCLOSE(IMERGER,3)
      CALL FCLOSE(I3,3)
*
*       Predict neighbour coordinates & velocities in case of KS collision.
   40 IF (NSYS.EQ.2) THEN
          NNB = LIST(1,ICOMP)
          CALL XVPRED(ICOMP,NNB)
*       Perform immediate KS regularization with close hierarchical body.
          IF (JCLOSE.GT.0) THEN
              ICOMP = IMERGER
              JCOMP = JCLOSE
              CALL KSREG
              GO TO 80
          END IF
      END IF
*
*       Initialize force polynomial for new single or third body (ICOMP).
      CALL FPOLY1(ICOMP,ICOMP,0)
      CALL FPOLY2(ICOMP,ICOMP,0)
*
*       See whether body #ICOMP should be added to NLIST.
      IF (T0(ICOMP) + STEP(ICOMP).LT.TLIST) THEN
          CALL NLMOD(ICOMP,1)
      END IF
      IF (NSYS.EQ.2) GO TO 80
*
*       Check chain parameters or DMIN in TRIPLE or QUAD.
   45 IF (ICH.GT.0) THEN
          WHICH1 = '  CHAIN '
          RCOLL = DMINC
*       Specify new membership (< 0 for SETSYS) and remove ghost from chain.
          NCH = -(NSYS - 1)
          JLIST(1) = IMERGER
          DO 50 L = 2,NSYS-1
              JLIST(L) = JLIST(L+1)
   50     CONTINUE
*       Copy well defined binding energy and skip explicit evaluation.
          EB = EBS
          CHCOLL = CHCOLL + EB
          GO TO 80
      ELSE IF (NSYS.EQ.3) THEN
          WHICH1 = ' TRIPLE '
          DMIN3 = MIN(DMIN3,RCOLL4)
          RCOLL = RCOLL4
      ELSE
          WHICH1 = '   QUAD '
          DMIN4 = MIN(DMIN4,RCOLL4)
          RCOLL = RCOLL4
      END IF
*
*       Obtain binding energy of the subsystem (ignore ghost).
      ZKE = 0.0D0
      POTS = 0.0D0
      DO 60 L = 1,NSYS
          I = JLIST(L)
          ZKE = ZKE + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                        XDOT(3,I)**2)
   60 CONTINUE
*
      DO 70 L = 1,NSYS-1
          DO 65 LL = L+1,NSYS
              I = JLIST(L)
              J = JLIST(LL)
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
              POTS = POTS + BODY(I)*BODY(J)/SQRT(RIJ2)
   65     CONTINUE
   70 CONTINUE
*
*       Form net energy correction for triple or quad case.
      EB = ENERGY - (0.5D0*ZKE - POTS)
*
*       Set global components for new KS regularization (ICOMP < JCOMP).
      ICOMP = MIN(IMERGER,I3)
      JCOMP = MAX(IMERGER,I3)
*
*       Initialize new KS pair.
      CALL KSREG
*
*       Update energy loss & collision counters.
   80 ECOLL = ECOLL + EB
      E(10) = E(10) + EB
      NPOP(8) = NPOP(8) + 1
      NCOLL = NCOLL + 1
*
*       Set IPHASE < 0 to ensure updating of time-step sequence.
      IPHASE = -1
*
      if(rank.eq.0)
     &WRITE (6,90)  WHICH1, NSYS, NAME1, NAME2, ZM, RCOLL, EB, VINF,
     &              ECC, DP
   90 FORMAT (/,A10,'  COLLISION  NSYS =',I3,'  NAME =',2I6,
     &             '  M =',F7.4,'  RCOLL =',1P,E8.1,'  EB =',E9.1,
     &             '  VINF =',0P,F5.1,'  ECC =',F9.5,'  DP =',1P,E9.1)
*
      RETURN
*
      END


