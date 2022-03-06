      SUBROUTINE DECORR(DE)
*
*
*       Collision/coalescence energy correction.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
*
*       Accumulate PN energy change for conservation purposes.
      ECOLL = ECOLL + DE
*
      WRITE (17,1)  TIME, NCH, DE, ECOLL
    1 FORMAT (' DECORR    T NCH DE EC  ',F9.4,I4,F11.7,F10.6)
      CALL FLUSH(17)
      RETURN
*
      END
