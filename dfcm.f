      SUBROUTINE DFCM(I,I1,FIRR,FD)
*
*
*       Irregular force on c.m. body.
*       -----------------------------
*
      INCLUDE 'common4.h'
      REAL*8  FIRR(3),FD(3),DX(3),DV(3),FP(6),FPD(6)
*
*
*       Initialize the perturbing force & derivative.
      DO 1 K = 1,6
          FP(K) = 0.0
          FPD(K) = 0.0
    1 CONTINUE
*
*       Define indicator for summing over each KS component rather than c.m.
      IFP = 0
*     NP = LIST(1,I1)
      RPERT2 = CMSEP2*R(I-N)**2
      I2 = I1 + 1
      KDUM = 0
*
*       Add contributions due to all perturbers.
*     DO 20 LL = 2,NP+1
*         K = LIST(LL,I1)
      DO 20 LL = 1,NSPERT
          K = JPERT(LL)
          IF (K.EQ.I) GO TO 20
*
          dr2 = 0.0
          DO 2 L = 1,3
              dx(L) = X(L,K) - X(L,I)
              dr2 = dr2 + dx(L)**2
    2     CONTINUE
*
          IF (K.LE.N) THEN
*             IF (dr2.GT.RPERT2) GO TO 10
              GO TO 3
          ELSE
              IPERT = 1
              KPERT = 1
              IF (dr2.GT.RPERT2) IPERT = 0
              IF (LIST(1,2*(K-N)-1).EQ.0) KPERT = 0
              IF (IPERT.EQ.0.AND.KPERT.EQ.0) GO TO 10
          END IF
*
*       Include the case of a perturbing KS pair.
          IF (KPERT.GT.0) THEN
*       Assume that components of pair #(K - N) have been resolved.
              KDUM = 2*(K - N) - 1
              K = KDUM
          END IF
*
*       Check c.m. approximation for current pair.
          IF (IPERT.EQ.0) GO TO 10
*
    3     IFP = 1
*
*       Evaluate perturbation on first component due to body #K.
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
*       See whether the second component of pair #J has been included.
          IF (K.GT.KDUM) GO TO 20
          K = K + 1
          GO TO 4
*
*       Sum over components of pair #J or its c.m. using c.m. approximation.
   10     dr2 = 0.0
          drdv = 0.0
          DO 12 L = 1,3
              dx(L) = X(L,K) - X(L,I)
              dv(L) = XDOT(L,K) - XDOT(L,I)
              dr2 = dr2 + dx(L)**2
              drdv = drdv + dx(L)*dv(L)
   12     CONTINUE
*
          dr2i = 1.0/dr2
          dr3i = BODY(K)*dr2i*SQRT(dr2i)
          drdv = 3.*drdv*dr2i
*
          DO 14 L = 1,3
              Firr(L) = Firr(L) + dx(L)*dr3i
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
          DO 25 K = 1,3
              FIRR(K) = FIRR(K) + (BODY(I1)*FP(K) +
     &                             BODY(I2)*FP(K+3))*BODYIN
              FD(K) = FD(K) + (BODY(I1)*FPD(K) +
     &                             BODY(I2)*FPD(K+3))*BODYIN
   25     CONTINUE
      END IF
*
      RETURN
*
      END
