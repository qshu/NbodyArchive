      SUBROUTINE KCPERT(I,I1,FIRR,FD)
*
*
*       Differential force correction on active KS from chain.
*       ------------------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      REAL*8  XI(3),XIDOT(3),FIRR(3),FD(3),FP(3),FPD(3),FPS(3),FDS(3),
     &        DV(3),FCM(3),FCMD(3)
*
*
*       See if perturber list contains chain c.m. body #ICH. 
      NNB1 = LIST(1,I1) + 1
      DO 20 L = 2,NNB1
          J = LIST(L,I1)
*       Finish search if perturber index exceeds chain c.m. index.
          IF (J.GT.ICH) GO TO 30
*
          IF (J.EQ.ICH) THEN
*       Obtain individual c.m. force with single particle approximation.
          A1 = X(1,J) - X(1,I)
          A2 = X(2,J) - X(2,I)
          A3 = X(3,J) - X(3,I)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
          DV(1) = XDOT(1,J) - XDOT(1,I)
          DV(2) = XDOT(2,J) - XDOT(2,I)
          DV(3) = XDOT(3,J) - XDOT(3,I)
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
*       Save c.m. contribution (to be subtracted).
          FCM(1) = A1*DR3I
          FCM(2) = A2*DR3I
          FCM(3) = A3*DR3I
          FCMD(1) = (DV(1) - A1*DRDV)*DR3I
          FCMD(2) = (DV(2) - A2*DRDV)*DR3I
          FCMD(3) = (DV(3) - A3*DRDV)*DR3I
*
*       Evaluate force & derivative on each KS component (saving the first).
              J1 = I1
              DO 10 KCOMP = 1,2
                  DO 5 K = 1,3
                      FP(K) = 0.0D0
                      FPD(K) = 0.0D0
                      XI(K) = X(K,J1)
                      XIDOT(K) = XDOT(K,J1)
    5             CONTINUE
                  CALL FCHAIN(J1,0,XI,XIDOT,FP,FPD)
                  IF (KCOMP.EQ.1) THEN
                      DO 8 K = 1,3
                          FPS(K) = FP(K)
                          FDS(K) = FPD(K)
    8                 CONTINUE
                  END IF
                  J1 = J1 + 1
   10         CONTINUE
*
*       Add differential contributions to irregular force & derivative.
              BODYIN = 1.0/BODY(I)
              DO 15 K = 1,3
                  FIRR(K) = FIRR(K) + ((BODY(I1)*FPS(K) +
     &                                  BODY(J1)*FP(K))*BODYIN - FCM(K))
                  FD(K) = FD(K) + ((BODY(I1)*FDS(K) +
     &                              BODY(J1)*FPD(K))*BODYIN - FCMD(K))
   15         CONTINUE
          END IF
   20 CONTINUE
*
   30 RETURN
*
      END
