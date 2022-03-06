      SUBROUTINE KSPERT(I1,NNB0,XI,FP)
*
*
*       Perturbation on KS pair.
*       ------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX,1)
      REAL*8  XI(6),XIDOT(3),FP(6),FPD(3)
*
*       Initialize the perturbing force.
      DO 10 K = 1,6
          FP(K) = 0.0D0
   10 CONTINUE
*
*       Set index of the last single perturber.
      NNB2 = NNB0 + 1
   15 IF (LIST(NNB2,I1).LE.N) GO TO 20
      NNB2 = NNB2 - 1
      IF (NNB2.GT.1) GO TO 15
*       Include special case of only c.m. perturbers.
      GO TO 30
*
*       Obtain the perturbation from single particles.
   20 DO 25 L = 2,NNB2
          K = LIST(L,I1)
          A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
          FP(1) = FP(1) + A1*A6
          FP(2) = FP(2) + A2*A6
          FP(3) = FP(3) + A3*A6
*       Perturbation on first component.
*
          A1 = X(1,K) - XI(4)
          A2 = X(2,K) - XI(5)
          A3 = X(3,K) - XI(6)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
          FP(4) = FP(4) + A1*A6
          FP(5) = FP(5) + A2*A6
          FP(6) = FP(6) + A3*A6
*
*       Perturbation on second component.
   25 CONTINUE
*
*       See whether to include any remaining c.m. perturbers.
      IF (NNB2.GT.NNB0) GO TO 40
*
   30 KDUM = 0
*       Dummy index to enable summation of c.m. or resolved components.
      NNB3 = NNB2 + 1
      DO 35 L = NNB3,NNB0+1
          K = LIST(L,I1)
          A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*       See whether c.m. approximation applies (ignore unperturbed case).
          J = K - N
          IF (RIJ2.GT.CMSEP2*R(J)**2) GO TO 33
*
*       Resolve pair #J and sum over individual components.
          CALL KSRES(J,J1,J2,RIJ2)
          KDUM = J1
          K = KDUM
   32     A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
   33     A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
          FP(1) = FP(1) + A1*A6
          FP(2) = FP(2) + A2*A6
          FP(3) = FP(3) + A3*A6
*       Perturbation on first component.
*
          A1 = X(1,K) - XI(4)
          A2 = X(2,K) - XI(5)
          A3 = X(3,K) - XI(6)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
          FP(4) = FP(4) + A1*A6
          FP(5) = FP(5) + A2*A6
          FP(6) = FP(6) + A3*A6
*
*       Perturbation on second component.
*
          IF (K.EQ.KDUM) THEN
              K = K + 1
              GO TO 32
          END IF
   35 CONTINUE
*
*       Check perturbation correction due to regularized chain.
   40 IF (NCH.GT.0) THEN
          DO 45 L = 2,NNB2
              J = LIST(I1,L)
              IF (J.GT.ICH) GO TO 50
              IF (J.EQ.ICH) THEN
*       Set velocity & force derivative (dummies needed by FCHAIN).
                  DO 42 K = 1,3
                      XIDOT(K) = XDOT(K,I1)
                      FPD(K) = 0.0
   42             CONTINUE
                  J1 = I1
                  CALL FCHAIN(J1,1,XI(1),XIDOT(1),FP(1),FPD(1))
                  J1 = J1 + 1
*                 CALL FCHAIN(J1,1,XI(4),XIDOT(1),FP(4),FPD(1))
                  GO TO 50
              END IF
   45     CONTINUE
      END IF
*
*       Set the relative perturbing force.
   50 DO 55 K = 1,3
          FP(K) = FP(K) - FP(K+3)
   55 CONTINUE
*
      RETURN
*
      END
