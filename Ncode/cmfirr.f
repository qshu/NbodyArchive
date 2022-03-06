      SUBROUTINE CMFIRR(I,IPAIR,XI,XIDOT,FIRR,FD)
*
*
*       Irregular c.m. force & derivative.
*       ----------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  XI(3),XIDOT(3),FIRR(3),FD(3),DX(3),DV(3),FP(6),FPD(6)
*
*
*       Initialize the perturbing force & derivative.
      DO 1 K = 1,6
          FP(K) = 0.0D0
          FPD(K) = 0.0D0
    1 CONTINUE
*
      I2 = 2*IPAIR
      I1 = I2 - 1
      KDUM = 0
*       Define indicator for summing over each KS component rather than c.m.
      IFP = 0
      RPERT2 = CMSEP2*R(IPAIR)**2
*
*       Force loop treats case I > N and any other c.m. neighbours.
      NNB0 = LIST(1,I)
      DO 20 LL = 2,NNB0+1
          K = LIST(LL,I)
          A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
*       Decide appropriate summation (c.m. approximation or components).
          IF (K.LE.N) THEN
              IF (RIJ2.GT.RPERT2) GO TO 10
              GO TO 3
          ELSE IF (RIJ2.LT.CMSEP2*R(K-N)**2) THEN
              KDUM = 2*(K - N) - 1
              IF (LIST(1,KDUM).EQ.0) THEN
                  KDUM = 0
              ELSE
                  K = KDUM
              END IF
          END IF
*
*       Check c.m. approximation for current pair.
          IF (RIJ2.GT.RPERT2) GO TO 10
*
*       Evaluate perturbation on first component due to body #K.
    3     IFP = 1
    4     dr2 = 0.0
          drdv = 0.0
          DO 5 L = 1,3
              dx(L) = X(L,K) - X(L,I1)
              dv(L) = XDOT(L,K) - XDOT(L,I1)
              dr2 = dr2 + dx(L)**2
              drdv = drdv + dx(L)*dv(L)
    5     CONTINUE
*
          dr2i = 1.0/dr2
          dr3i = BODY(K)*dr2i*SQRT(dr2i)
          drdv = 3.0*drdv*dr2i
*
          DO 6 L = 1,3
              FP(L) = FP(L) + dx(L)*dr3i
              FPD(L) = FPD(L) + (dv(L) - dx(L)*drdv)*dr3i
    6     CONTINUE
*
*       Evaluate perturbation on second component due to body #K.
          dr2 = 0.0
          drdv = 0.0
          DO 7 L = 1,3
              dx(L) = X(L,K) - X(L,I2)
              dv(L) = XDOT(L,K) - XDOT(L,I2)
              dr2 = dr2 + dx(L)**2
              drdv = drdv + dx(L)*dv(L)
    7     CONTINUE
*
          dr2i = 1.0/dr2
          dr3i = BODY(K)*dr2i*SQRT(dr2i)
          drdv = 3.0*drdv*dr2i
*
          DO 8 L = 1,3
              FP(L+3) = FP(L+3) + dx(L)*dr3i
              FPD(L+3) = FPD(L+3) + (dv(L) - dx(L)*drdv)*dr3i
    8     CONTINUE
*
          IF (K.GT.KDUM) GO TO 20
          K = K + 1
          GO TO 4
*
*       Sum over components of pair #J or its c.m. using c.m. approximation.
   10     dr2 = 0.0
          drdv = 0.0
          DO 12 L = 1,3
              dx(L) = X(L,K) - XI(L)
              dv(L) = XDOT(L,K) - XIDOT(L)
              dr2 = dr2 + dx(L)**2
              drdv = drdv + dx(L)*dv(L)
   12     CONTINUE
*
          dr2i = 1.0/dr2
          dr3i = BODY(K)*dr2i*SQRT(dr2i)
          drdv = 3.0*drdv*dr2i
*
          DO 14 L = 1,3
              FIRR(L) = FIRR(L) + dx(L)*dr3i
              FD(L) = FD(L) + (dv(L) - dx(L)*drdv)*dr3i
   14     CONTINUE
          IF (K.EQ.KDUM) THEN
              K = K + 1
              GO TO 10
          END IF
   20 CONTINUE
*
*       Add mass-weighted perturbations to force & first derivative.
      IF (IFP.GT.0) THEN
          BODYIN = 1.0/BODY(I)
          DO 30 K = 1,3
              FIRR(K) = (BODY(I1)*FP(K) + BODY(I2)*FP(K+3))*BODYIN +
     &                                                           FIRR(K)
              FD(K) = (BODY(I1)*FPD(K) + BODY(I2)*FPD(K+3))*BODYIN +
     &                                                             FD(K)
   30     CONTINUE
      END IF
*
*       Check force correction due to regularized chain.
      IF (NCH.GT.0) THEN
          CALL KCPERT(I,I1,FIRR,FD)
      END IF 
*
      RETURN
*
      END
