      SUBROUTINE KSPOLY(IPAIR,IMOD)
*
*
*       Initialization of KS polynomials.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW/  RANGE,ISLOW(10)
      REAL*8  A1(3,4),A(8),FP(4),FPDOT(3)
*
*
*       Specify indices of c.m. & components and perturber membership (+ 1).
      I = N + IPAIR
      I2 = 2*IPAIR
      I1 = I2 - 1
      NNB1 = LIST(1,I1) + 1
*
*       Initialize variables for accumulating contributions.
      DO 10 K = 1,3
          FP(K) = 0.0D0
          FPDOT(K) = 0.0D0
   10 CONTINUE
*       Set dummy index for summation of c.m. or resolved components.
      JDUM = 0
*
*       Obtain the perturbation & first derivative.
      DO 40 L = 2,NNB1
          J = LIST(L,I1)
          IF (J.GT.N) THEN
*       See whether c.m. approximation applies.
              RIJ2 = (X(1,J) - X(1,I))**2 + (X(2,J) - X(2,I))**2 +
     &                                      (X(3,J) - X(3,I))**2
              K = J - N
              IF (RIJ2.LT.CMSEP2*R(K)**2.AND.LIST(1,2*K-1).GT.0) THEN
                  JDUM = 2*K - 1
                  J = JDUM
              END IF
          END IF
*
*       Sum over both components (reverse sign for second component).
   20     II = I1
          DO 30 KCOMP = 1,2
              DO 25 K = 1,3
                  A(K) = X(K,J) - X(K,II)
                  A(K+3) = XDOT(K,J) - XDOT(K,II)
   25         CONTINUE
*       Current velocities are predicted in routines INTGRT, KSMOD, etc.
              RIJ2 = A(1)*A(1) + A(2)*A(2) + A(3)*A(3)
              A8 = BODY(J)/(RIJ2*SQRT(RIJ2))
              A9 = 3.0D0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))/RIJ2
              IF (KCOMP.EQ.2) A8 = -A8
*
              DO 28 K = 1,3
                  FP(K) = FP(K) + A(K)*A8
                  FPDOT(K) = FPDOT(K) + (A(K+3) - A(K)*A9)*A8
   28         CONTINUE
              II = I2
   30     CONTINUE
*
          IF (J.EQ.JDUM) THEN
              J = J + 1
              GO TO 20
          END IF
   40 CONTINUE
*
*       Check for external perturbations.
      IF (KZ(14).NE.0) THEN
          Q1 = X(1,I1) - X(1,I2)
          Q3 = X(3,I1) - X(3,I2)
          CALL XTRNLP(Q1,Q3,FP)
*
*       Use same formalism for the first derivative (omit Coriolis force).
          VX = XDOT(1,I1) - XDOT(1,I2)
          VZ = XDOT(3,I1) - XDOT(3,I2)
          CALL XTRNLP(VX,VZ,FPDOT)
      END IF
*
*       Transform to regularized force derivative using T' = R.
      DO 45 K = 1,3
          FPDOT(K) = R(IPAIR)*FPDOT(K)
   45 CONTINUE
*
*       Save the relative perturbation.
      FP(4) = SQRT(FP(1)**2 + FP(2)**2 + FP(3)**2)
      GAMMA(IPAIR) = FP(4)*R(IPAIR)**2/BODY(I)
*
*       Set current transformation matrix.
      A1(1,1) =  U(1,IPAIR)
      A1(1,2) = -U(2,IPAIR)
      A1(1,3) = -U(3,IPAIR)
      A1(1,4) =  U(4,IPAIR)
      A1(2,1) =  U(2,IPAIR)
      A1(2,2) =  U(1,IPAIR)
      A1(2,3) = -U(4,IPAIR)
      A1(2,4) = -U(3,IPAIR)
      A1(3,1) =  U(3,IPAIR)
      A1(3,2) =  U(4,IPAIR)
      A1(3,3) =  U(1,IPAIR)
      A1(3,4) =  U(2,IPAIR)
*
*       Construct regularized polynomials from explicit derivatives.
      TDOT2(IPAIR) = 0.0D0
      TDOT3(IPAIR) = 0.0D0
      TDOT4 = 0.0D0
      TDOT5 = 0.0D0
      HDOT(IPAIR) = 0.0D0
      D1HDOT(IPAIR) = 0.0D0
      D2HDOT(IPAIR) = 0.0D0
      D3HDOT(IPAIR) = 0.0D0
*
*       Scale perturbing force & first derivative by modification factor.
      IF (IMOD.GT.1) THEN
          ZMOD = FLOAT(ISLOW(IMOD))
          DO 50 K = 1,3
              FP(K) = ZMOD*FP(K)
              FPDOT(K) = ZMOD*FPDOT(K)
   50     CONTINUE
      END IF
*
*       Form regularized force & two derivatives of time & binding energy.
      DO 60 K = 1,4
          A(K) = A1(1,K)*FP(1) + A1(2,K)*FP(2) + A1(3,K)*FP(3)
          A(K+4) = A1(1,K)*FPDOT(1) + A1(2,K)*FPDOT(2) +
     &                                A1(3,K)*FPDOT(3)
          FU(K,IPAIR) = 0.5D0*H(IPAIR)*U(K,IPAIR) + 0.5D0*R(IPAIR)*A(K)
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0D0*U(K,IPAIR)*UDOT(K,IPAIR)
          TDOT3(IPAIR) = TDOT3(IPAIR) + 2.0D0*UDOT(K,IPAIR)**2 +
     &                                      2.0D0*U(K,IPAIR)*FU(K,IPAIR)
          HDOT(IPAIR) = HDOT(IPAIR) + 2.0D0*UDOT(K,IPAIR)*A(K)
          D1HDOT(IPAIR) = D1HDOT(IPAIR) + 2.0D0*FU(K,IPAIR)*A(K)
   60 CONTINUE
*
*       Set regularized velocity matrix (Levi-Civita matrix not required).
      A1(1,1) =  UDOT(1,IPAIR)
      A1(1,2) = -UDOT(2,IPAIR)
      A1(1,3) = -UDOT(3,IPAIR)
      A1(1,4) =  UDOT(4,IPAIR)
      A1(2,1) =  UDOT(2,IPAIR)
      A1(2,2) =  UDOT(1,IPAIR)
      A1(2,3) = -UDOT(4,IPAIR)
      A1(2,4) = -UDOT(3,IPAIR)
      A1(3,1) =  UDOT(3,IPAIR)
      A1(3,2) =  UDOT(4,IPAIR)
      A1(3,3) =  UDOT(1,IPAIR)
      A1(3,4) =  UDOT(2,IPAIR)
*
*       Include the whole (L*F)' term in explicit derivatives of FU & H.
      DO 65 K = 1,4
          A(K+4) = A(K+4) + A1(1,K)*FP(1) + A1(2,K)*FP(2) +
     &                                      A1(3,K)*FP(3)
          D1HDOT(IPAIR) = D1HDOT(IPAIR) + 2.0D0*UDOT(K,IPAIR)*A(K+4)
          FUDOT(K,IPAIR) = 0.5D0*(H(IPAIR)*UDOT(K,IPAIR) +
     &                            HDOT(IPAIR)*U(K,IPAIR) +
     &                            TDOT2(IPAIR)*A(K) + R(IPAIR)*A(K+4))
          D2HDOT(IPAIR) = D2HDOT(IPAIR) + 2.0*FUDOT(K,IPAIR)*A(K) +
     &                                            4.0*FU(K,IPAIR)*A(K+4)
          TDOT4 = TDOT4 + 2.0*FUDOT(K,IPAIR)*U(K,IPAIR) +
     &                                     6.0*FU(K,IPAIR)*UDOT(K,IPAIR)
   65 CONTINUE
*
*       Form higher derivatives by boot-strapping.
      DO 70 K = 1,4
          D2U(K,IPAIR) = 0.5*(H(IPAIR)*FU(K,IPAIR) +
     &                   D1HDOT(IPAIR)*U(K,IPAIR) + TDOT3(IPAIR)*A(K)) +
     &                   TDOT2(IPAIR)*A(K+4) + HDOT(IPAIR)*UDOT(K,IPAIR)
          D3U(K,IPAIR) = 0.5*(H(IPAIR)*FUDOT(K,IPAIR) + TDOT4*A(K) +
     &                        D2HDOT(IPAIR)*U(K,IPAIR)) +
     &                   1.5*(HDOT(IPAIR)*FU(K,IPAIR) +
     &                        D1HDOT(IPAIR)*UDOT(K,IPAIR) +
     &                        TDOT3(IPAIR)*A(K+4))
          D3HDOT(IPAIR) = D3HDOT(IPAIR) + 2.0*D2U(K,IPAIR)*A(K) +
     &                                         6.0*FUDOT(K,IPAIR)*A(K+4)
          TDOT5 = TDOT5 + 2.0*D2U(K,IPAIR)*U(K,IPAIR) +
     &             8.0*FUDOT(K,IPAIR)*UDOT(K,IPAIR) + 6.0*FU(K,IPAIR)**2
   70 CONTINUE
*
*       Check maximum square step (soft binaries & weak hyperbolic pairs).
      IF (ABS(H(IPAIR)).GT.ECLOSE) THEN
          A2 = 0.5/ABS(H(IPAIR))
      ELSE
          A2 = MIN(R(IPAIR)/BODY(I),0.5/ABS(H(IPAIR)))
      END IF
*
*       Assign a conservative initial step of half the standard value.
      DTU = 0.5*ETAU*SQRT(A2)/(1.0 + 1000.0*GAMMA(IPAIR))**0.333
      DTAU(IPAIR) = DTU
*
*       Convert regularized step to physical units and initialize time T0.
      STEP(I1) = ((((TDOT5*DTU/120.0D0 + TDOT4/24.0D0)*DTU +
     &               TDOT3(IPAIR)*ONE6)*DTU + 0.5D0*TDOT2(IPAIR))*DTU +
     &                                                     R(IPAIR))*DTU
      IF (IMOD.GT.1) THEN
          STEP(I1) = ZMOD*STEP(I1)
      END IF
      T0(I1) = TIME
*
*       Set regularized backward times (TAU merely acts as reference).
      T0U(IPAIR) = TAU
      T1U(IPAIR) = TAU - DTU
      T2U(IPAIR) = TAU - 2.0*DTU
      T3U(IPAIR) = TAU - 3.0*DTU
*
*       Convert from derivatives to differences.
      DO 80 K = 1,4
          D1U(K,IPAIR) = (ONE6*D3U(K,IPAIR)*DTU -
     &                          0.5D0*D2U(K,IPAIR))*DTU + FUDOT(K,IPAIR)
          D2U(K,IPAIR) = 0.5D0*D2U(K,IPAIR) - 0.5D0*D3U(K,IPAIR)*DTU
          D3U(K,IPAIR) = ONE6*D3U(K,IPAIR)
          FU(K,IPAIR) = 0.5D0*FU(K,IPAIR)
          FUDOT(K,IPAIR) = ONE6*FUDOT(K,IPAIR)
*       Half regularized force and sixth the derivative for fast predictor.
   80 CONTINUE
*
      D1HDOT(IPAIR) = (ONE6*D3HDOT(IPAIR)*DTU - 0.5D0*D2HDOT(IPAIR))*DTU
     &                                                   + D1HDOT(IPAIR)
      D2HDOT(IPAIR) = 0.5D0*D2HDOT(IPAIR) - 0.5D0*D3HDOT(IPAIR)*DTU
      D3HDOT(IPAIR) = ONE6*D3HDOT(IPAIR)
*
*       See whether to include relative motion in the time-step list.
      IF (T0(I1) + STEP(I1).LT.TLIST) THEN
          CALL NLMOD(I1,1)
      END IF
*
      RETURN
*
      END
