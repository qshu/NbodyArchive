      SUBROUTINE SHORT(LENGTH,NXTLST)
!VAMPINST include
*
*
*       Short time-step list.
*       ---------------------
*
      INCLUDE 'common6.h'
      INTEGER  NXTLST(NMAX)
*
*
*       Decide between modifying list or increasing membership.
!VAMPINST subroutine_start
    1 NNB = LSHORT(1)
      IF (NNB.GT.10) THEN
          K1 = NNB - 9
          K2 = NNB - (K1 - 1)
*       Reduce membership below 10 by removing old entries.
          DO 5 K = K1,K2
              LSHORT(K) = LSHORT(K+K1)
    5     CONTINUE
          NNB = NNB - K1
*       Add new KS candidates and c.m. particles unless already present.
          DO 10 K = 1,LENGTH
              J = NXTLST(K)
              DO 8 L = 2,NNB+1
                  IF (LSHORT(L).EQ.J) GO TO 10
    8         CONTINUE
              IF (STEP(J).LT.SMIN) THEN
                  NNB = NNB + 1
                  LSHORT(NNB+1) = J
              END IF
   10     CONTINUE
          LSHORT(1) = NNB
      ELSE
*       Check existing members before adding current small step particles.
          DO 20 K = 1,LENGTH
              J = NXTLST(K)
              DO 15 L = 2,NNB+1
                  IF (LSHORT(L).EQ.J) GO TO 20
   15         CONTINUE
              IF (STEP(J).LT.SMIN) THEN
                  NNB = NNB + 1
                  LSHORT(NNB+1) = J
*                 LSHORT(1) = NNB
              END IF
   20     CONTINUE
          LSHORT(1) = NNB
      END IF
*
!VAMPINST return
          RETURN
      RETURN
*
!VAMPINST subroutine_end
      END
