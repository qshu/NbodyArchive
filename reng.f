* ======================================
* The PN relativistic energy of a binary
* ======================================
*
      SUBROUTINE RENG (ipair,i1,ksx,ksv)
*
      INCLUDE 'common6.h'
*
*      implicit none
*
      integer ipair,i1,i2
      integer i,j,k
*
      real*8 ksx(3),ksv(3),n12(3)
      real*8 r1,r2,r3,r4,r5
      real*8 v1,vv,v3,v4,v5
      real*8 rdot,rdot2,rdot4
      real*8 c,c2,c3,c4,c5
      real*8 M_both,mu,eta
      real*8 E_1PN,E_2PN,E_25PN
*
* Set label for second particle
      i2 = i1 + 1
*
* Zero energies for safety
      E_1PN = 0.0
      E_2PN = 0.0
      E_25PN = 0.0
*
* Powers of c
      c = CLIGHT
      c2 = c*c
      c3 = c2*c
      c4 = c2*c2
      c5 = c4*c
*
* Masses, separation, velocity and unit vector
      M_both = BODY(i1) + BODY(i2)
      M_2 = M_both*M_both
      M_3 = M_2*M
      mu = BODY(i1)*BODY(i2)/M_both
      eta = mu/M_both
*
      r2 = ksx(1)*ksx(1) + ksx(2)*ksx(2) + ksx(3)*ksx(3)
      r1 = sqrt(r2)
      r3 = r2*r1
      r4 = r2*r2
      r5 = r4*r1
*
      vv = ksv(1)*ksv(1) + ksv(2)*ksv(2) + ksv(3)*ksv(3)
      v = sqrt(vv)
      v3 = vv*v
      v4 = vv*vv
      v5 = v4*v
      v6 = v3*v3
*
      DO i = 1,3
         n12(i) = ksx(i)/r1
      END DO
      rdot = n12(1)*ksv(1) + n12(2)*ksv(2) + n12(3)*ksv(3)
      rdot2 = rdot*rdot
      rdot4 = rdot2*rdot2
*
* =========================================================
* Calculate relativistic energy per unit mass at each order
* =========================================================
*
* ----------------------
* The 1PN energy - 1/c^2
* ----------------------
*
      E_1PN = (1.0/c2)*(M_2/(2.0*r2) + 3.0*(1.0 - 3.0*eta)*v4/8.0
     &     + (3.0 + eta)*vv*M_both/(2.0*r1) + eta*M_both*rdot2/(2.0*r1))
*
* ----------------------
* The 2PN energy - 1/c^4
* ----------------------
*
      E_2PN = (1.0/c4)*(-(2.0 + 15.0*eta)*M_3/(4.0*r3)
     &     + 5.0*(1.0 - 7.0*eta + 13.0*eta*eta)*v6/16.0
     &     + (14.0 - 55.0*eta + 4.0*eta*eta)*M_2*vv/(8.0*r2)
     &     + (4.0 + 69.0*eta + 12.0*eta*eta)*M_2*rdot2/(8.0*r2)
     &     + (21.0 - 23.0*eta - 27.0*eta*eta)*M_both*v4/(8.0*r1)
     &     + eta*(1.0 - 15.0*eta)*M_both*vv*rdot2/(4.0*r1)
     &     - 3.0*eta*(1.0 - 3.0*eta)*M_both*rdot4/(8.0*r1))
*
* ------------------------
* The 2.5PN energy - 1/c^5
* ------------------------
*
      E_25PN = (1.0/c5)*(8.0*eta*M_2*rdot*vv/(5.0*r2))
*
* =====================================
* Apply relativistic energy corrections
* =====================================
*
      RELENERGY1(ipair) = mu*E_1PN
      RELENERGY2(ipair) = mu*E_2PN
      RELENERGY2_5(ipair) = mu*E_25PN
*
      RETURN
*
      END
