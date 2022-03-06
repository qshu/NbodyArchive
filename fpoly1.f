


      SUBROUTINE FPOLY1(I1,I2)
*
*
*       Force & first derivative.
*       -------------------------
*
      INCLUDE 'common1.h'
      integer i1,i2,i,j,k
      REAL*4  A(6)
      REAL * 8 a7, a8, a9, mrinv
*
*
*       Loop over all bodies or just one single body (I1 = I2).
      DO 50 I = I1,I2
*
*       Initialize force & first derivative.
      DO 2 K = 1,3
          F(K,I) = 0.0
          FDOT(K,I) = 0.0
    2 CONTINUE
      phi(i) = 0.0
*
*       Sum over all the other particles.
      DO 10 J = 1,N
          IF (J.EQ.I) GO TO 10
*
          DO 5 K = 1,3
              A(K) = X(K,J) - X(K,I)
              A(K+3) = XDOT(K,J) - XDOT(K,I)
    5     CONTINUE
*    
          if(i .le. nbh .and. j.le.nbh) then
              A7 = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
          else if(i .le. nbh .or. j.le.nbh) then
              A7 = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3)+meps2)
          else
              A7 = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3) + EPS2)
          endif
          mrinv = BODY(J)*SQRT(A7)
          phi(i) = phi(i) - mrinv
          A8 = mrinv*A7
          A9 = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A7
*
          DO 8 K = 1,3
              F(K,I) = F(K,I) + A(K)*A8
              FDOT(K,I) = FDOT(K,I) + (A(K+3) - A(K)*A9)*A8
    8     CONTINUE
   10 CONTINUE
*     write(6,666) i, phi(i), (f(k,i), k = 1, 3)
 666  format(i5, 4e16.7)
   50 CONTINUE
*
      RETURN
*
      END


      SUBROUTINE FPOLY1a(i)
*
*
*       Force & first derivative.
*       from/to BH      
*       -------------------------
*
      INCLUDE 'common1.h'
      integer i,j,k, nend
      REAL*4  A(6)
      REAL * 8 a7, a8, a9, mrinv
*
*
*       Loop over all bodies or just one single body (I1 = I2).
*
*       Initialize force & first derivative.
*
*       Sum over all the other particles.
      nend = nbh
      if(i .le. nbh) then
         nend = n
         DO  K = 1,3
            F(K,I) = 0.0
            FDOT(K,I) = 0.0
         enddo
         phi(i) = 0.0
      endif
      DO 10 J = 1,Nend
         IF (J.EQ.I) GO TO 10
*        
         DO 5 K = 1,3
            A(K) = X(K,J) - X(K,I)
            A(K+3) = XDOT(K,J) - XDOT(K,I)
 5       CONTINUE
*        
         if(i .le. nbh .and. j.le.nbh) then
            A7 = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
         else if(i .le. nbh .or. j.le.nbh) then
            A7 = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3)+meps2)
         endif
         mrinv = BODY(J)*SQRT(A7)
         phi(i) = phi(i) - mrinv
         A8 = mrinv*A7
         A9 = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A7
*        
         DO 8 K = 1,3
            F(K,I) = F(K,I) + A(K)*A8
            FDOT(K,I) = FDOT(K,I) + (A(K+3) - A(K)*A9)*A8
 8       CONTINUE
 10   CONTINUE
*     
      END


