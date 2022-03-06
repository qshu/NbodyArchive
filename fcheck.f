      SUBROUTINE FCHECK(I,FF,FD,FP)
*
*
*       Force & first derivative.
*       -------------------------
*
      INCLUDE 'common4.h'
      REAL*8  A(9),F1(3),F1DOT(3),FF(3),FD(3),FP
*
*
*
*       Initialize force & first derivative of body #I.
      DO 10 K = 1,3
          FF(K) = 0.0D0
          FD(K) = 0.0D0
   10 CONTINUE
      FP = 0.0D0
*
*       Obtain force & first derivative by summing over all bodies.
      KDUM = 0
      DO 30 JDUM = IFIRST,NTOT
          IF (JDUM.EQ.I) GO TO 30
          J = JDUM
          IF (J.GT.N) THEN
              JPAIR = J - N
*       Use c.m. approximation for unperturbed binary.
              IF (LIST(1,2*JPAIR-1).GT.0) THEN
                  KDUM = 2*JPAIR - 1
                  J = KDUM
              END IF
          END IF
*
   12     DO 15 K = 1,3
              A(K) = X(K,J) - X(K,I)
              A(K+3) = XDOT(K,J) - XDOT(K,I)
   15     CONTINUE
*
          A(7) = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
          A(8) = BODY(J)*A(7)*SQRT(A(7))
          A(9) = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A(7)
*
          DO 20 K = 1,3
              F1(K) = A(K)*A(8)
              F1DOT(K) = (A(K+3) - A(K)*A(9))*A(8)
   20     CONTINUE
*
          DO 25 K = 1,3
              FF(K) = FF(K) + F1(K)
              FD(K) = FD(K) + F1DOT(K)
   25     CONTINUE
          FP = FP - BODY(J)/SQRT(A(1)**2 + A(2)**2 + A(3)**2)
*
          IF (J.EQ.KDUM) THEN
              J = J + 1
              GO TO 12
          END IF
   30 CONTINUE
*
      RETURN
*
      END
