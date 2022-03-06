      SUBROUTINE FICORR(I,DM)
*
*
*       Local force corrections due to mass loss.
*       -----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  A(6)
*      REAL*8  XBACK(3,LMAX)
*
*       Correct force & first derivative for neighbours of body #I.
      call jpred(I,TIME,TIME)
*     --11/14/13 19:13-lwang-test---------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      POTJ = 0.0
*     --11/14/13 19:13-lwang-end----------------------------------------*
      NNB1 = LIST(1,I) + 1
      DO 20 L = 2,NNB1
          J = LIST(L,I)
          NNB2 = LIST(1,J) + 1
*
*       Ensure that body #I is a neighbour of body #J.
          DO 1 K = 2,NNB2
              IF (LIST(K,J).EQ.I) GO TO 5
    1     CONTINUE
          GO TO 20
*
    5     RIJ2 = 0.0D0
          A7 = 0.0D0
*
          call jpred(j,time,time)
          DO 10 K = 1,3
              A(K) = X(K,I) - X(K,J)
              A(K+3) = XDOT(K,I) - XDOT(K,J)
              RIJ2 = RIJ2 + A(K)**2
              A7 = A7 + A(K)*A(K+3)
   10     CONTINUE
*
          A7 = 3.0*A7/RIJ2
          A8 = DM/(RIJ2*SQRT(RIJ2))
*
          DO 15 K = 1,3
              A(K+3) = (A(K+3) - A7*A(K))*A8
*     --11/14/13 19:12-lwang-test---------------------------------------*
***** Note:------------------------------------------------------------**
*              POTJ = POTJ + BODY(J)/SQRT(RIJ2)
*     --11/14/13 19:12-lwang-end----------------------------------------*
              F(K,J) = F(K,J) - 0.5*A(K)*A8
              FI(K,J) = FI(K,J) - A(K)*A8
              FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
              D1(K,J) = D1(K,J) - A(K+3)
              FIDOT(K,J) = FIDOT(K,J) - A(K+3)
   15     CONTINUE
   20 CONTINUE
*
*     --11/14/13 19:13-lwang-test---------------------------------------*
***** Note:------------------------------------------------------------**
c$$$*       Update the potential energy loss.
c$$$      EMDOT = EMDOT - DM*POTJ
c$$$*
c$$$*       Include kinetic energy correction.
c$$$      VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
c$$$      EMDOT = EMDOT + 0.5*DM*VI2
c$$$
*     --11/14/13 19:13-lwang-end----------------------------------------*
c$$$*       Backup all predicted x and xdot to xback and xdback
c$$$      NNB = NNB1 -1 
c$$$      DO II = 1, NNB
c$$$         IL = LIST(II+1,I)
c$$$         XBACK(1,II) = X0(1,IL)
c$$$         XBACK(2,II) = X0(2,IL)
c$$$         XBACK(3,II) = X0(3,IL)
c$$$         X0(1,IL) = X(1,IL)
c$$$         X0(2,IL) = X(2,IL)         
c$$$         X0(3,IL) = X(3,IL)
c$$$      END DO
c$$$      
c$$$*       Obtain the potential energy due to all particles.
c$$$      POTJ = 0.0D0
c$$$      DO 30 J = IFIRST,NTOT
c$$$          IF (J.EQ.I) GO TO 30
c$$$          RIJ2 = 0.0D0
c$$$          DO 25 K = 1,3
c$$$              RIJ2 = RIJ2 + (X(K,I) - X0(K,J))**2
c$$$   25     CONTINUE
c$$$          POTJ = POTJ + BODY(J)/SQRT(RIJ2)
c$$$   30 CONTINUE
c$$$*       Revert X0 and X0dot to original values
c$$$      DO II = 1, NNB
c$$$         IL = LIST(II+1,I)
c$$$         X0(1,IL) = XBACK(1,II)
c$$$         X0(2,IL) = XBACK(2,II)
c$$$         X0(3,IL) = XBACK(3,II)
c$$$      END DO
c$$$*
c$$$*       Update the potential energy loss.
c$$$      EMDOT = EMDOT - DM*POTJ
c$$$*
c$$$*       Include kinetic energy correction.
c$$$      VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
c$$$      EMDOT = EMDOT + 0.5*DM*VI2
*
*       See whether linearized tidal terms should be included.
      IF (KZ(14).GT.0.AND.KZ(14).LT.3) THEN
          ECDOT = ECDOT - 0.5*DM*(TIDAL(1)*X(1,I)**2 +
     &                            TIDAL(3)*X(3,I)**2)
      END IF
*
*       Check optional Plummer potential.
      IF (KZ(14).EQ.4.OR.KZ(14).EQ.3) THEN
          RI2 = AP2
          DO 50 K = 1,3
              RI2 = RI2 + X(K,I)**2
  50      CONTINUE
          ECDOT = ECDOT - DM*MP/SQRT(RI2)
      END IF
*
      RETURN
*
      END
