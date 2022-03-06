      SUBROUTINE RESORT(LENGTH,NXTLST,TM)
*
*     determine next block particles without assuming
*     sorted time step list (R.Sp.)
*     (advantage for large N - but see correction in routine INTGRT)
*     tm minimum time
      INCLUDE 'common6.h'
      INTEGER LENGTH,J,NXTLST(*)
*         find next minimum t+dt
         TM = TIMENW(IFIRST)
         DO 5 J = IFIRST, NTOT
            IF (TIMENW(J).LT.TM)TM=TIMENW(J)
  5      CONTINUE
*         find particles on level TM
      INXT = 0
      DO 10 J = IFIRST, NTOT
         IF(TIMENW(J).GT.TM) THEN
            CONTINUE
         ELSE
            INXT = INXT + 1
            NXTLST(INXT) = J
         END IF
 10   CONTINUE
      LENGTH = INXT
*
      RETURN
      END
