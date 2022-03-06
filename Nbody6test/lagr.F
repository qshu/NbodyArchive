      SUBROUTINE LAGR(C)
*
*
*       Lagrangian radii.
*       -----------------
*
      INCLUDE 'common6.h'
      REAL*8  R2(NMAX)
      REAL*8  C(3),FLAGR(11),RLAGR(11),AVMASS(11)
*       Lagrangian radii at 1,2,5,7.5,10,15,20,50,75,85 % of total mass.
      DATA FLAGR/0.01,0.02,0.05,0.10,0.20,0.30,0.40,0.50,0.75,0.90,1.00/
*  correspond to: 1% , 2% , 5% , 10%, 20%, 30% , 40% ,50%, 75% ,90%,100%
*
*
*       Set square radii of all single particles & c.m. bodies.
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
*  Determine the Lagrangian radii for specified mass fractions.
*     RLAGR = Lagrangian radius
*     AVMASS = average mass of a spherical shell with radius R2(I)
*     AVMRC = average mass of particles within core radius RC
*     ISHELL = particle counter within a shell
      AVMRC = 0.D0
      RC2 = RC*RC
      NC = 0
      ZM = 0.0D0
      ZMH = 0.5D0*ZMASS
      I = 0
*
      DO 15 J = 1, 11
         AVMASS(J) = 0.D0
         ISHELL = 0
 20      I = I + 1
         IM = JLIST(I)
         ZM = ZM + BODY(IM)
         AVMASS(J) = AVMASS(J) + BODY(IM)
         IF(R2(I).LT.RC2)THEN
            AVMRC = AVMRC + BODY(IM)
            NC = NC + 1
         END IF
         ISHELL = ISHELL + 1
         IF (I.LT.NTOT.AND.ZM.LT.FLAGR(J)*ZMASS) GOTO 20
         RLAGR(J) = SQRT(R2(I))
         AVMASS(J) = AVMASS(J)/ISHELL
 15   CONTINUE
      AVMRC = AVMRC/NC
*
*       Determine half-mass radius separately.
      ZM = 0.0
      ZMH = 0.5D0*ZMASS
      I = 0
   30 I = I + 1
      IM = JLIST(I)
      ZM = ZM + BODY(IM)
      IF (ZM.LT.ZMH) GO TO 30
*
*       Replace approximate half-mass radius by actual value.
      RSCALE = SQRT(R2(I))
*
#ifdef PARALLEL
      if(rank.eq.0)then
#endif
*       Check output options (line printer or unit 7 or both).
      IF (KZ(7).EQ.2.OR.KZ(7).EQ.4) THEN
          WRITE (6,40)  (FLAGR(K),K=1,11),' <RC'
   40     FORMAT (/,3X,'M/MT:  ',1P,11D9.2,A4)
          WRITE (6,41)  (RLAGR(K),K=1,11),RC
   41     FORMAT (/,3X,'RLAGR: ',1P,12D9.2)
          WRITE (6,42)  (AVMASS(K),K=1,11),AVMRC
   42     FORMAT (/,3X,'AVMASS:',1P,12D9.2)
      END IF
*       Write on File
      IF (KZ(7).GE.3) THEN
            WRITE (7,100) TIME+TOFF, (RLAGR(K),K=1,11),
     &           (AVMASS(K),K=1,11),RC,AVMRC
 100         FORMAT(1X,1P,25(1X,D15.7))
      END IF
#ifdef PARALLEL
      end if
#endif
*
      RETURN
*
      END
