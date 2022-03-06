      SUBROUTINE LAGR(C)
*
*
*       Lagrangian radii.
*       -----------------
*
      INCLUDE 'common6.h'
      REAL  R2
      COMMON/WORK1/  R2(NMAX)
      REAL*8  C(3),FLAGR(10)
      DATA FLAGR/-1.9,-1.7,-1.5,-1.3,-1.1,-.9,-.7,-.5,-.3,-.1/
*
*
*       Set square radii of single particles & c.m. bodies.
      NP = 0
      DO 10 I = IFIRST,NTOT
          NP = NP + 1
          R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                                  (X(3,I) - C(3))**2
          JLIST(NP) = I
   10 CONTINUE
*
*       Sort square distances with respect to the centre C.
      CALL SORT1(NP,R2,JLIST)
*
*       Determine the half-mass radius.
      ZM = 0.0D0
      ZMH = 0.5*ZMASS
      I = 0
   20 I = I + 1
      IM = JLIST(I)
      ZM = ZM + BODY(IM)
      IF (ZM.LT.ZMH) GO TO 20
*
*       Replace approximate half-mass radius by actual value.
      RSCALE = SQRT(R2(I))
*
      if(rank.eq.0)then
*       Check output options (line printer or unit 7 or both).
      IF (KZ(7).EQ.2.OR.KZ(7).EQ.4) THEN
          WRITE (6,40)
     &             (.5*LOG10(R2(MAX(NINT(10.0**FLAGR(K)*NP),1))),K=1,10)
   40     FORMAT ('  LAGR:  ',10F7.3)
      END IF
      end if
*
#ifdef PARALLEL
      if(rank.eq.0)then
#endif
      IF (KZ(7).GE.3) THEN
          WRITE (70+rank)  TIME,
     &            (0.5*LOG10(R2(MAX(NINT(10.0**FLAGR(K)*NP),1))),K=1,10)
      END IF
#ifdef PARALLEL
      end if
#endif
*
      RETURN
*
      END
