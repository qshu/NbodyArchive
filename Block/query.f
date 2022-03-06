      SUBROUTINE QUERY(I,I1,SEMI)
*
*
*       Analysis of large KS perturbation.
*       ----------------------------------
*
      INCLUDE 'common6.h'
*
*
      IPAIR = KVEC(I1)
      GI = GAMMA(IPAIR)
*
      IF (NCH.EQ.0.AND.SEMI.LT.5.0*RMIN.AND.NAME(I).GE.0.AND.
     &    (KZ(30).NE.0.AND.KZ(30).NE.-2)) THEN
*
*       Search the perturbers (singles and binaries) inside 10*SEMI.
          RSEP = MAX(SEMI,RMIN)
          RCR2 = 100.0*RSEP**2
          NNB1 = LIST(1,I1) + 1
          JCL = 0
          JCL2 = 0
          FMAX = 0.0
          FMAX2 = 0.0
          JP = 0
          DO 85 L = 2,NNB1
              J = LIST(L,I1)
              RD = (X(1,I)-X(1,J))*(XDOT(1,I)-XDOT(1,J)) +
     &             (X(2,I)-X(2,J))*(XDOT(2,I)-XDOT(2,J)) +
     &             (X(3,I)-X(3,J))*(XDOT(3,I)-XDOT(3,J))
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
      WRITE (60,444) TIME, J, NAME(J), SEMI, SQRT(RIJ2), RD, GI
  444 FORMAT (' SEARCH    T J NM A RIJ RD G ',F8.3,2I6,1P,5E10.2)
      CALL FLUSH(60)
*             RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
*    &                                      (X(3,I) - X(3,J))**2
*       Save possible strong perturber just outside limit.
              IF (RIJ2.GT.RCR2) THEN
                  IF (RIJ2.LT.4.0*RCR2) THEN
                      JP = J
                      RP2 = RIJ2
                  END IF
                  GO TO 85
              END IF
              FIJ = (BODY(I) + BODY(J))/RIJ2
              IF (FIJ.GT.FMAX) THEN
                  FMAX2 = FMAX
                  JCL2 = JCL
                  FMAX = FIJ
                  JCL = J
                  RRD = RD
                  RJ2 = RIJ2
                  VIJ2 = (XDOT(1,I) - XDOT(1,J))**2 +
     &                   (XDOT(2,I) - XDOT(2,J))**2 +
     &                   (XDOT(3,I) - XDOT(3,J))**2
              ELSE IF (NAME(J).GE.0.AND.FIJ.GT.FMAX2) THEN
                  FMAX2 = FIJ
                  JCL2 = J
              END IF
   85     CONTINUE
*
*       Quit on dominant external perturber.
          IF (JCL.GT.0.AND.JP.GT.0) THEN
              JLIST(1) = JP
              CALL FPERT(I,JCL,1,JLIST,PERT)
              GP = PERT*RP2/(BODY(I) + BODY(JCL))
              IF (GP.GT.0.10) THEN
              WRITE (60,688) TIME, NAME(JP), BODY(JP)*SMU, RRD, GP
  688         FORMAT (' CHAIN ACCEPT    T NM MP RD GP ',
     &                  F8.3, I6,F7.2,1P,2E10.2)
              CALL FLUSH(60)
              END IF
          END IF
      END IF
*
      RETURN
*
      END
