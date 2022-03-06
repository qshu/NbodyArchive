      SUBROUTINE NBCHECK(NREG,IREG)
      INCLUDE 'common6.h'
      DIMENSION IREG(NREG)
      INTEGER NBREG(NMAX)
      REAL R2(NMAX)
*
      IOVR=0
      IDIS=0
      DO 1001 L=1,NREG-1
*
      I = IREG(L)
      NNBI = LIST(1,I)
      RS2I = RS(IREG(I))**2
*
      DO 1002 LL=L+1,NREG
*
      II = IREG(LL)
      NNBII = LIST(1,II)
      RR = DSQRT((X(1,I)-X(1,II))**2 + (X(2,I)-X(2,II))**2 +
     *        (X(3,I)-X(3,II))**2)
*
      IF(RR.LE.RS(I)+RS(II))THEN
*      Find neighbours of other particles (II) which are not nnb of (I)
*
          DO 50 J = 2,NNBI+1
*
              DO 55 JJ = 2,NNBII+1
*
                  IF(LIST(JJ,II).EQ.LIST(J,I))GOTO 55


      IOVR = IOVR + 1
      ELSE
*      Neighbour spheres are distinct
      IDIS = IDIS + 1
      END IF
*
 1002 CONTINUE
 1001 CONTINUE
      PRINT*,' NREG,IOVR,IDIS=',NREG,IOVR,IDIS
      END
