      SUBROUTINE FCORR(I,DM,KW)
*
*
*       Total force corrections due to mass loss.
*       -----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  A(9)
      LOGICAL IKICK
*
*
*       Save the velocity components and square velocity.
      VI2 = 0.0
      DO 1 K = 1,3
          A(K+6) = XDOT(K,I)
          VI2 = VI2 + XDOT(K,I)**2
    1 CONTINUE
*
*       Include velocity kick in case of neutron star, BH or mass-less SN.
      IF (KW.EQ.13.OR.KW.EQ.14.OR.KW.EQ.15) THEN
          IKICK = .TRUE.
*       Distinguish between single star and binary.
          IF (I.LE.N) THEN
              CALL KICK(I,1)
          ELSE
              IPAIR = I - N
              CALL KICK(IPAIR,0)
          END IF
      ELSE
          IKICK = .FALSE.
      END IF
*
*       Define consistent c.m. coordinates & velocities for KS mass loss.
      IF (I.GT.N) THEN
          I2 = 2*(I - N)
          I1 = I2 - 1
          VF2 = 0.0
          DV2 = 0.0
          BODYI = BODY(I)
          IF (BODY(I).EQ.0.0D0) BODYI = BODY(I1) + BODY(I2)
          DO 8 K = 1,3
              X(K,I) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/BODYI
              XDOT(K,I) = (BODY(I1)*XDOT(K,I1) + BODY(I2)*XDOT(K,I2))/
     &                                                           BODYI
              X0(K,I) = X(K,I)
              X0DOT(K,I) = XDOT(K,I)
              VF2 = VF2 + XDOT(K,I)**2
              DV2 = DV2 + (XDOT(K,I) - A(K+6))**2
    8     CONTINUE
          VFAC = SQRT(VF2/VI2)
      END IF
*
*       Correct potential energy, all forces & first derivatives.
      POTJ = 0.0D0
      DO 40 J = IFIRST,NTOT
          IF (J.EQ.I) GO TO 40
          RIJ2 = 0.0D0
          RIJDOT = 0.0D0
          RDVDOT = 0.0D0
*
          DO 10 K = 1,3
              A(K) = X(K,I) - X(K,J)
              A(K+3) = A(K+6) - XDOT(K,J)
              RIJ2 = RIJ2 + A(K)**2
              RIJDOT = RIJDOT + A(K)*A(K+3)
              RDVDOT = RDVDOT + A(K)*(XDOT(K,I) - A(K+6))
   10     CONTINUE
*
          RIJ = SQRT(RIJ2)
          POTJ = POTJ + BODY(J)/RIJ
          A3 = 1.0/(RIJ2*RIJ)
          A4 = BODY(I)*A3
          A5 = DM*A3
          A6 = 3.0*RIJDOT/RIJ2
          A7 = 3.0*RDVDOT/RIJ2
*
          DO 15 K = 1,3
              A(K+3) = (A(K+3) - A6*A(K))*A5
              IF (IKICK) THEN
*       Include FDOT corrections due to increased velocity.
                  A(K+3) = A(K+3) + (XDOT(K,I) - A(K+6))*A4
                  A(K+3) = A(K+3) - A7*A(K)*A4
              END IF
   15     CONTINUE
*
*       Use neighbour list to distinguish irregular & regular terms.
          NNB1 = LIST(1,J) + 1
          DO 30 L = 2,NNB1
              IF (LIST(L,J).EQ.I) THEN
                  DO 20 K = 1,3
                      F(K,J) = F(K,J) - 0.5*A(K)*A5
                      FI(K,J) = FI(K,J) - A(K)*A5
                      FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
                      D1(K,J) = D1(K,J) - A(K+3)
                      FIDOT(K,J) = FIDOT(K,J) - A(K+3)
   20             CONTINUE
*       Reduce the step and see whether body #J should be added to NLIST.
*                 STEP(J) = MAX(0.5D0*STEP(J),TIME - T0(J))
                  IF (T0(J) + STEP(J).LT.TLIST) THEN
                      CALL NLMOD(J,1)
                  END IF
                  GO TO 40
              ELSE IF (LIST(L,J).GT.I) THEN
                  DO 25 K = 1,3
                      F(K,J) = F(K,J) - 0.5*A(K)*A5
                      FR(K,J) = FR(K,J) - A(K)*A5
                      FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
                      D1R(K,J) = D1R(K,J) - A(K+3)
                      FRDOT(K,J) = FRDOT(K,J) - A(K+3)
   25             CONTINUE
                  GO TO 40
              END IF
   30     CONTINUE
   40 CONTINUE
*
*       Update the potential and kinetic energy loss.
      ECDOT = ECDOT - DM*POTJ + 0.5*DM*VI2
*
*       Check energy correction due to velocity kick.
      IF (I.GT.N.AND.IKICK) THEN
          ECDOT = ECDOT - 0.5*BODY(I)*VI2*(VFAC**2 - 1.0)
      END IF
*
*       See whether tidal terms should be included.
      IF (KZ(14).GT.0) THEN
          ECDOT = ECDOT - 0.5*DM*(TIDAL(1)*X(1,I)**2 +
     &                            TIDAL(3)*X(3,I)**2)
      END IF
*
      E(12) = EMDOT
*
      RETURN
*
      END
