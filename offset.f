      SUBROUTINE OFFSET(DTOFF)
*
*
*       Offset of global times.
*       -----------------------
*
      INCLUDE 'common6.h'
*
*
*       Update the global offset time.
    1 TOFF = TOFF + DTOFF
*
*       Reduce all individual times and epochs by offset interval.
      DO 10 I = 1,NTOT
          T0(I) = T0(I) - DTOFF
          T0R(I) = T0R(I) - DTOFF
          TEV(I) = TEV(I) - DTOFF
          TEV0(I) = TEV0(I) - DTOFF
          EPOCH(I) = EPOCH(I) - DTOFF*TSTAR
   10 CONTINUE
*
*       Set new global times.
      TIME = TIME - DTOFF
      TADJ = TADJ - DTOFF
      TNEXT = TNEXT - DTOFF
      TPREV = TPREV - DTOFF
      IF (KZ(19).GT.2) THEN
          TPLOT = TPLOT - DTOFF
          TMDOT = TMDOT - DTOFF
      END IF
*
*       See whether more reductions are needed.
      IF (TIME.GE.TOFF) GO TO 1
*
*       Activate control indicator for new scheduling.
      IPHASE = -1
*
      RETURN
*
      END
