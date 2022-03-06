      SUBROUTINE RESORT(NXTLEN,NXTLST)
*
*     determine next block particles without assuming
*     sorted time step list (R.Sp.)
*     (advantage for large N - but see correction in routine INTGRT)
*     tm minimum time
      INCLUDE 'common6.h'
      INTEGER NXTLEN,NXTLST(NMAX)
*
*         find particles on level TIME
      NXTLEN = 0
      DO 10 J = IFIRST, NTOT
         IF(TIMENW(J).LE.TIME) THEN
            NXTLEN = NXTLEN + 1
            NXTLST(NXTLEN) = J
         END IF
 10   CONTINUE
*
      RETURN
      END
