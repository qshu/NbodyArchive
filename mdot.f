      SUBROUTINE MDOT(IM)
*
*
*       Mass loss from evolving stars.
*       ------------------------------
*
      INCLUDE 'common6.h'
      REAL  LUMS(6),TSCLS(6),ZM,TM,TN,RM,TE,AGE
      CHARACTER*8  WHICH1
      DATA  AU  /2.0627E+05/
*
*
*       Define current time (unit of million years) and update counter.
      TPHYS = TIME*TSTAR
      AGE = TPHYS
      NMDOT = NMDOT + 1
      IM = 0
      KS = 0
*
*       Find global index and original name of next star to be checked.
    1 DO 10 J = 1,N
          LJ = NAME(J)
          IF (TEV(LJ).LT.TIME) THEN
              I = J
              LI = LJ
              GO TO 20
          END IF
   10 CONTINUE
*
*       Determine new evolution time (current TMDOT may be escaped star).
      KW = 0
      GO TO 60
*
*       Skip possible c.m. of compact subsystem or zero mass ghost.
   20 IF (NSUB.GT.0.AND.(BODY(I).GT.BODY1.OR.BODY(I).LE.0.0D0)) THEN
          TEV(LI) = TEV(LI) + 0.01*TCR
          GO TO 1
      END IF
*
*       See whether body #I is a ghost member of merger configuration.
      IF (NMERGE.GT.0) THEN
          IF (BODY(I).LE.0.0D0) THEN
              TEV(LI) = TEV(LI) + 0.01*TCR
              GO TO 1
          END IF
      END IF
*
*       Copy the appropriate initial mass (ZAMS or unevolved value).
      IF (LI.LE.500) THEN
          ZM0 = BODY0(LI)
      ELSE
          ZM0 = BODY(I)
      END IF
*
*       Obtain stellar parameters at current epoch.
      ZM = ZM0*ZMBAR
      CALL STAR(ZM,TM,TN,TSCLS,LUMS)
      CALL HRDIAG(ZM,AGE,TM,TN,LUMS,TSCLS,RM,TE,KW)
*
*       Set new mass & radius in scaled units and form mass difference.
      ZMNEW = ZM/ZMBAR
      RNEW = 0.005*RM*(RBAR/AU)
      DM = BODY(I) - ZMNEW
      DMR = ABS(DM/BODY(I))
*
      IF (I.LT.IFIRST.AND.DMR.GT.0.01) THEN
          KSPAIR = KVEC(I)
*       Terminate KS pair or merger configuration and redetermine index.
          IF (NAME(N+KSPAIR).GT.0) THEN
              CALL RESOLV(KSPAIR,3)
              CALL KSTERM
              KS = 1
          ELSE
              CALL RESET
*       Restore phase indicator (NB! IPHASE = -1 skips RESOLVE in KSTERM).
              IPHASE = 0
          END IF
          GO TO 1
      END IF
*
*       Base new time scale for changes in radius & mass on stellar type.
      IF (KW.EQ.1) THEN
          DTM = 0.1*TM
          DTR = TM - AGE
      ELSE IF (KW.EQ.2) THEN
          DTM = 0.1*TSCLS(1)
          DTR = TM + TSCLS(1) - AGE
      ELSE IF (KW.EQ.3) THEN
          DTM = 0.1*TSCLS(2)
          DTR = TM + TSCLS(1) + TSCLS(2) - AGE
      ELSE IF (KW.EQ.4) THEN
          DTM = 0.1*TSCLS(3)
          DTR = TM + TSCLS(1) + TSCLS(2) + TSCLS(3) - AGE
      ELSE IF (KW.EQ.5) THEN
          DTM = 0.1*TSCLS(6)
          DTR = TN + TSCLS(6) - AGE
      ELSE IF (KW.EQ.6) THEN
          DTM = 0.1*(TN - TM)
          DTR = TN - AGE
      ELSE IF (KW.GE.7) THEN
          DTM = 1.0E+10
          DTR = 1.0E+10
      END IF
*
*       Include special treatment of Thorne-Zytkow objects.
      IF (KSTAR(LI).EQ.9) THEN
*
*       Reduce mass by 10 % at constant radius and adopt dm/dt = 1 Sun/Myr.
          KW = 9
          ZMNEW = 0.9*BODY(I)
          RNEW = RADIUS(I)
          DTM = 0.1*BODY(I)*ZMBAR/TSTAR
*
*       Terminate mass loss at white dwarf limit.
          IF (ZMNEW*ZMBAR.LT.1.4) THEN
              KW = 8
              ZMNEW = 1.4/ZMBAR
              RNEW = 5.0E-05*(RBAR/AU)
              DTM = 1.0E+10
          END IF
*
          DM = BODY(I) - ZMNEW
          ZM = ZMNEW*ZMBAR
          DMR = ABS(DM/BODY(I))
          DTR = DTM
      END IF
*
*       Perform neighbour force corrections if mass loss is significant.
      IF (DMR.GT.0.01) THEN
*
*       Set the new mass.
          BODY(I) = ZMNEW
*
*       Accumulate the total mass loss (solar units) and reduce cluster mass.
          DMSUN = DM*ZMBAR
          ZMDOT = ZMDOT + DMSUN
          ZMASS = ZMASS - DM
      IF (KZ(19).GT.3) THEN
          WRITE (6,30)  I, NAME(I), KW, KSTAR(LI), BODY(I)*ZMBAR,
     &                  DMSUN, ZMDOT, TSCALE*TIME
   30     FORMAT (' MDOT:   I NM KW K* MS DMS ZMDOT T6 ',
     &                      4I5,F6.1,F7.2,F7.1,F8.1)
      END IF
*
*       Update event counter & mass loss.
          IF (KW.EQ.3) THEN
              NRG = NRG + 1
              ZMRG = ZMRG + DMSUN
          ELSE IF (KW.EQ.4) THEN
              NHE = NHE + 1
              ZMHE = ZMHE + DMSUN
          ELSE IF (KW.EQ.5.OR.KW.EQ.9) THEN
              NRS = NRS + 1
              ZMRS = ZMRS + DMSUN
          ELSE IF (KW.EQ.7) THEN
              NWD = NWD + 1
              ZMWD = ZMWD + DMSUN
          ELSE IF (KW.EQ.8) THEN
              NSN = NSN + 1
              ZMSN = ZMSN + DMSUN
          END IF
*
*       Predict coordinates & velocities of body #I and its neighbours.
          NNB = LIST(1,I)
          CALL XVPRED(I,0)
          CALL XVPRED(I,NNB)
*
*       Perform local or global force & energy corrections depending on DM/M.
          IF (DMR.LT.0.02) THEN
              CALL FICORR(I,DM)
          ELSE
              CALL FCORR(I,DM)
*
*       Initialize new polynomials of neighbours & #I for neutron star case.
              IF (DMSUN.GT.0.5.OR.KW.EQ.8) THEN
                  NNB2 = LIST(1,I) + 2
                  LIST(NNB2,I) = I
*
*       Obtain new F & FDOT, then F2DOT & F3DOT and time-steps.
                  DO 40 L = 2,NNB2
                      J = LIST(L,I)
                      T0(J) = TIME
                      DO 35 K = 1,3
                          X0DOT(K,J) = XDOT(K,J)
   35                 CONTINUE
                      CALL FPOLY1(J,J,0)
   40             CONTINUE
*
                  DO 45 L = 2,NNB2
                      J = LIST(L,I)
                      CALL FPOLY2(J,J,0)
   45             CONTINUE
              END IF
*
*       Switch on time-step list indicator after significant mass loss.
              IM = 1
          END IF
      END IF
*
*       Choose minimum of 10 % of evolution time and remaining time interval.
      DTM = MIN(DTM,DTR)
*       Impose a lower limit and convert time interval to scaled units.
      DTM = MAX(DTM,1.0D-03)/TSTAR
*
*       Set new time for checking R & M and update R & classification type.
      TEV(LI) = TIME + DTM
      RADIUS(I) = RNEW
      KW0 = KSTAR(LI)
      KSTAR(LI) = KW
*
      IF (KZ(19).GT.3.AND.(KW0.NE.KW.OR.DMR.GT.0.01)) THEN
          IF (KW0.NE.KW) THEN
              WHICH1 = ' TYPE   '
          ELSE
              WHICH1 = ' MASS   '
          END IF
          WRITE (6,50)  WHICH1, TPHYS, I, NAME(I), DMR, KW0, KW,
     &                  ZM0*ZMBAR, ZM, RADIUS(I)/(0.005*(RBAR/AU)),EMDOT
   50     FORMAT (' NEW',A8,' TPHYS I NAM DM/M KW0 KW M0 M R EMD ',
     &                        F7.1,2I5,F6.2,2I3,2F6.1,1P,E9.1,0P,F10.5)
      END IF
*
*       See if former KS pair can be regularized again.
      IF (KS.GT.0) THEN
          ICOMP = IFIRST
          JCOMP = IFIRST + 1
          RIJ2 = (X(1,ICOMP) - X(1,JCOMP))**2 +
     &           (X(2,ICOMP) - X(2,JCOMP))**2 +
     &           (X(3,ICOMP) - X(3,JCOMP))**2
          IF (RIJ2.LT.RMIN2) THEN
              CALL KSREG
          END IF
      END IF
*
*       Determine the time for next stellar evolution check.
   60 TMDOT = 1.0E+10
      DO 70 J = 1,N
          LJ = NAME(J)
          IF (TEV(LJ).LT.TMDOT) THEN
              TMDOT = TEV(LJ)
          END IF
   70 CONTINUE
*
*       Update the maximum single body mass but skip compact subsystems.
      IF (KW.GE.5.AND.NSUB.EQ.0) THEN
          BODY1 = 0.0
          DO 80 J = 1,N
              BODY1 = MAX(BODY1,BODY(J))
   80     CONTINUE
      END IF
*
*       See if any other stars should be considered.
      IF (TMDOT.LT.TIME) GO TO 1
*
      RETURN
*
      END
