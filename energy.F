      SUBROUTINE ENERGY(ZKIN,POT)
*
*
*       Total energy.
*       -------------
*
      INCLUDE 'common1.h'
*
*
      integer i,j,k,jmin
      real*8 dt, a1, a2, zkin, pot, potj, a3
*       Predict current coordinates & velocities (XDOT restored later).
      IF (TIME.EQ.0.0D0) GO TO 10
*
      DO 2 I = 1,N
          DT = TIME - T0(I)
          IF (DT.EQ.0.0D0) GO TO 2
c          write(6,*) 'predict', i, time, dt, t0(i), step(i)
          A1 = 0.2*DT
          A2 = DT/24.
          DO 1 K = 1,3
             
c              X(K,I) = ((((D3(K,I)*a1 + D2(k,i))*a2 +
c     &                    FDOT(K,I))*DT + F(K,I))*DT + X0DOT(K,I))*DT +
c     &                                                           X0(K,I)
c              XDOT(K,I)  = (((D3(K,I)*a2 + ONE6*D2(k,i))*DT +
c     &                              3.0*FDOT(K,I))*DT + 2.0*F(K,I))*DT +
c     &                                                        X0DOT(K,I)
              X(K,I) = (((FDOT(K,I))*DT + F(K,I))*DT + X0DOT(K,I))*DT +
     &                                                         X0(K,I)
              XDOT(K,I)  = ((3.0*FDOT(K,I))*DT + 2.0*F(K,I))*DT +
     &                                                     X0DOT(K,I)
    1     CONTINUE
    2 CONTINUE
*
*       Calculate the potential energy.
   10 ZKIN = 0.0
      POT = 0.0
*     i = 0
*     if (i.ne.0) then
#if (! USE_GRAPE      )
      I = 1
   15 JMIN = I + 1
      POTJ = 0.0
      DO 18 J = JMIN,N
          A1 = X(1,I) - X(1,J)
          A2 = X(2,I) - X(2,J)
          A3 = X(3,I) - X(3,J)
          if(i .le. nbh) then
              POTJ = POTJ + BODY(J)/SQRT(A1*A1 + A2*A2 + A3*A3)
          else
              POTJ = POTJ+BODY(J)/SQRT(A1*A1 + A2*A2 + A3*A3 + EPS2)
          endif
   18 CONTINUE
      POT = POT + BODY(I)*POTJ
      I = I + 1
      IF (I.LT.N) GO TO 15
#else
*     else
         do 12 i = 1, n
            pot = pot+body(i)*phi(i)
12       continue
         pot = -pot*0.5
*     endif
#endif
     
*
*       Obtain the kinetic energy.
      DO 20 I = 1,N
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
   20 CONTINUE
*
      ZKIN = 0.5*ZKIN
*       Total energy = ZKIN - POT.
*
      print*,' energy kin,pot=',ZKIN,POT
      RETURN
*
      END
