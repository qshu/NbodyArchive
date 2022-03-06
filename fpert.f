      SUBROUTINE FPERT(I,J,NP,PERT)
*
*
*       Perturbing force on dominant bodies.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  FP(3)
      REAL*8  XI1(3),VI1(3),FP1(3),FDUM(3)
*
*
      DO 1 K = 1,3
          FP(K) = 0.0D0
    1 CONTINUE
*
*       Obtain perturbation on body #I & J due to NP members of JLIST.
      DO 10 KK = 1,NP
          K = JLIST(KK)
          IF (K.EQ.I.OR.K.EQ.J) GO TO 10
          L = I
    5     A1 = X(1,K) - X(1,L)
          A2 = X(2,K) - X(2,L)
          A3 = X(3,K) - X(3,L)
          RIJ2 = A1**2 + A2**2 + A3**2
          A4 = BODY(K)/(RIJ2*SQRT(RIJ2))
          IF (L.EQ.J) A4 = -A4
          FP(1) = FP(1) + A1*A4
          FP(2) = FP(2) + A2*A4
          FP(3) = FP(3) + A3*A4
          IF(L.EQ.I) THEN
              L = J
              GO TO 5
          END IF
   10 CONTINUE
*
      DO 2 K=1,3
      XI1(K)=X(K,I)
      VI1(K)=XDOT(K,I)
      FP1(K)=0.0D0
      FDUM(K)=0.0D0
    2 CONTINUE
      CALL MASSCENTROBJECT(I,XI1,VI1,FP1,FDUM)
      CALL DRAGFORCE(I,XI1,VI1,FP1,FDUM,FREG,2)
*
      DO 22 K=1,3
  22  FP(K) = FP(K) + FP1(K)
*
      DO 33 K=1,3
      XI1(K)=X(K,J)
      VI1(K)=XDOT(K,J)
      FP1(K)=0.0D0
      FDUM(K)=0.0D0
   33 CONTINUE
      CALL MASSCENTROBJECT(J,XI1,VI1,FP1,FDUM)
      CALL DRAGFORCE(J,XI1,VI1,FP1,FDUM,FREG,2)
*
      DO 44 K=1,3
  44  FP(K) = FP(K) - FP1(K)
*
*
      PERT = SQRT(FP(1)**2 + FP(2)**2 + FP(3)**2)
*
      RETURN
*  
      END
