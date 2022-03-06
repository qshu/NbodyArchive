      IMPLICIT INTEGER(I-K,M-N),REAL*8(A-H,O-Z),LOGICAL(L)
      REAL*8 VK(8)
      PRINT*,' Enter NSEED, NCASE:'
      READ*,NSEED,NCASE

      DISP = 265.0

      TWOPI = 8.D0*DATAN(1.D0)
      PRINT*,' TWOPI=',TWOPI
      IDUM1 = -NSEED
      DO 111 KK = 1,NCASE
*       Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
      DO 2 K = 1,2
          X1 = RAN2(IDUM1)
          X2 = RAN2(IDUM1)
*       Generate two velocities from polar coordinates S & THETA.
          S = DISP*SQRT(-2.0*LOG(1.0 - X1))
          THETA = TWOPI*X2
          VK(2*K-1) = S*COS(THETA) 
          VK(2*K) = S*SIN(THETA) 
    2 CONTINUE

      DO 4 K = 3,4
          X1 = RAN2(IDUM1)
          X2 = RAN2(IDUM1)
*       Generate two velocities from polar coordinates S & THETA.
          IF(K.EQ.1)S = DISP*SQRT(-2.0*LOG(1.0 - X1))
          THETA = TWOPI*X2
          VK(2*K-1) = S*COS(THETA)
          VK(2*K) = S*SIN(THETA)
    4 CONTINUE

      WRITE(6,112) KK,(VK(K),K=1,8)
 112  FORMAT(1X,1P,I10,8E13.5)
 111  CONTINUE

      STOP
      END
