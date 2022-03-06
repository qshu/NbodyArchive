      SUBROUTINE FCORR(I,DM)
*
*
*       Global force corrections due to masss loss.
*       -------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  A(6)
*
*
*       Correct potential energy, forces & first derivatives.
      POTJ = 0.0D0
      DO 40 J = IFIRST,NTOT
          IF (J.EQ.I) GO TO 40
          RIJ2 = 0.0D0
          RIJDOT = 0.0D0
*         RDVDOT = 0.0D0
*
          DO 10 K = 1,3
              A(K) = X(K,I) - X(K,J)
              A(K+3) = XDOT(K,I) - XDOT(K,J)
              RIJ2 = RIJ2 + A(K)**2
              RIJDOT = RIJDOT + A(K)*A(K+3)
*             RDVDOT = RDVDOT + A(K)*(XDOT(K,I) - A(K+6))
   10     CONTINUE
*
          RIJ = SQRT(RIJ2)
          POTJ = POTJ + BODY(J)/RIJ
          A3 = 1.0/(RIJ2*RIJ)
*         A4 = BODY(I)*A3
          A5 = DM*A3
          A6 = 3.0*RIJDOT/RIJ2
*         A7 = 3.0*RDVDOT/RIJ2
*
          DO 15 K = 1,3
              A(K+3) = (A(K+3) - A6*A(K))*A5
*             IF (A0.GT.1.0) THEN
*       Include FDOT corrections due to increased velocity.
*                 A(K+3) = A(K+3) + (XDOT(K,I) - A(K+6))*A4
*                 A(K+3) = A(K+3) - A7*A(K)*A4
*             END IF
   15     CONTINUE
*
*       Use neighbour list to distinguish irregular & regular terms.
          NNB = LIST(1,J) + 1
          DO 30 L = 2,NNB
              IF (LIST(L,J).EQ.I) THEN
                  DO 20 K = 1,3
                      F(K,J) = F(K,J) - 0.5*A(K)*A5
                      FI(K,J) = FI(K,J) - A(K)*A5
                      FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
                      D1(K,J) = D1(K,J) - A(K+3)
   20             CONTINUE
*       Reduce the step and see whether body #J should be added to NLIST.
*                 STEP(J) = STEP(J) - 0.5*(T0(J) + STEP(J) - TIME)
*                 IF (T0(J) + STEP(J).LT.TLIST) THEN
*                     CALL NLMOD(J,1)
*                 END IF
                  GO TO 40
              ELSE IF (LIST(L,J).GT.I) THEN
                  DO 25 K = 1,3
                      F(K,J) = F(K,J) - 0.5*A(K)*A5
                      FR(K,J) = FR(K,J) - A(K)*A5
                      FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
                      D1R(K,J) = D1R(K,J) - A(K+3)
   25             CONTINUE
                  GO TO 40
              END IF
   30     CONTINUE
   40 CONTINUE
*
*       Update the potential energy loss.
      EMDOT = EMDOT - DM*POTJ
*
*       Include kinetic energy correction.
      VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
      EMDOT = EMDOT + 0.5*DM*VI2
*
*       See whether tidal terms should be included.
      IF (KZ(14).GT.0) THEN
          EMDOT = EMDOT - 0.5*DM*(TIDAL(1)*X(1,I)**2 +
     &                            TIDAL(3)*X(3,I)**2)
      END IF
*
      RETURN
*
      END
