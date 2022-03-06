      SUBROUTINE STUMPF(IPAIR,Z)
*
*       Modified Stumpff functions.
*       ---------------------------
*
*       Developed by Seppo Mikkola 10/1997.
*       ns is the max order included in z^n.
*
      INCLUDE 'common6.h'
      parameter (ns=12)
      real*8 c(12)
*
*
*       Generate Stumpff functions for argument Z (SCOEFF set by START).
      c(ns)=1-z*SCOEFF(12)
      c(ns-1)=1-z*SCOEFF(11)
      do i=ns-2,1,-1
          c(i)=1-z*c(i+2)*SCOEFF(i)
      end do
*
*       Copy the five first coefficients (for predictor & corrector).
      DO 2 I = 1,5
          SF(I,IPAIR) = C(I)
*         SF(I,IPAIR) = 1.0
    2 CONTINUE
*
*       Generate Stumpff functions for argument 4*Z.
      Z4 = 4.0D0*Z
      c(ns)=1-z4*SCOEFF(12)
      c(ns-1)=1-z4*SCOEFF(11)
      do i=ns-2,1,-1
          c(i)=1-z4*c(i+2)*SCOEFF(i)
      end do
*
*       Copy c_4 & c_5 to sixth & seventh location (for time integration).
      DO 3 I = 6,7
          SF(I,IPAIR) = C(I-2)
*         SF(I,IPAIR) = 1.0
    3 CONTINUE
*
      RETURN
*
      END
