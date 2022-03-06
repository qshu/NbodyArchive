      SUBROUTINE XTRNLF(I,FIRR,FD)
*
*
*       External force & first derivative.
*       ----------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  FIRR(3),FD(3),XI(3),VI(3)
*
*
*       Choose the relevant galactic tidal force (standard or point mass).
      IF (KZ(14).LE.2) THEN
          FIRR(1) = FIRR(1) + TIDAL(1)*X(1,I) + TIDAL(4)*XDOT(2,I)
          FIRR(2) = FIRR(2) - TIDAL(4)*XDOT(1,I)
          FIRR(3) = FIRR(3) + TIDAL(3)*X(3,I)
          FD(1) = FD(1) + TIDAL(1)*XDOT(1,I) + TIDAL(4)*FIRR(2)
          FD(2) = FD(2) - TIDAL(4)*FIRR(1)
          FD(3) = FD(3) + TIDAL(3)*XDOT(3,I)
      END IF
*
      RETURN
*
      END
