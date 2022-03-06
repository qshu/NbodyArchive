      SUBROUTINE XTRNLF(XI,XIDOT,FIRR,FREG,FD,FDR,KCASE)
*
*
*       External force & first derivative.
*       ----------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  XI(3),XIDOT(3),FIRR(3),FREG(3),FD(3),FDR(3)
*
*
*       See whether to include the galactic tidal force.
      IF (TIDAL(1).GT.0.0) THEN
          FIRR(1) = FIRR(1) + TIDAL(4)*XIDOT(2)
          FIRR(2) = FIRR(2) - TIDAL(4)*XIDOT(1)
          FD(1) = FD(1) + TIDAL(4)*(FIRR(2) + FREG(2))
          FD(2) = FD(2) - TIDAL(4)*(FIRR(1) + FREG(1))
*       Add smooth part to the regular components (KCASE = 1 in REGINT).
          IF (KCASE.GT.0) THEN
              FREG(1) = FREG(1) + TIDAL(1)*XI(1)
              FREG(3) = FREG(3) + TIDAL(3)*XI(3)
              FDR(1) = FDR(1) + TIDAL(1)*XIDOT(1)
              FDR(3) = FDR(3) + TIDAL(3)*XIDOT(3)
          END IF
      END IF
*
      RETURN
*
      END
