      SUBROUTINE RSSORT(IMIN,NXTLEN,NXTLST)
!VAMPINST include
*
      INCLUDE 'common6.h'
      INTEGER NXTLEN,NXTLST(NMAX)
*
!VAMPINST subroutine_start
      TMIN = TIMENW(IFIRST)
      IMIN = IFIRST
*
      DO 3 J = IFIRST+1, NTOT
         IF(TIMENW(J).LT.TMIN)THEN
         TMIN = TIMENW(J)
         IMIN = J
         END IF
  3   CONTINUE
*
      DO 5 J = IFIRST, NTOT
         IF(DABS(TIMENW(J)-TMIN).LT.DTK(40)) THEN
            NXTLEN = NXTLEN + 1
            NXTLST(NXTLEN) = J
         END IF
  5   CONTINUE
*
!VAMPINST return
          RETURN
      RETURN
!VAMPINST subroutine_end
      END
