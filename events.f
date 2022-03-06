      SUBROUTINE EVENTS
*
*
*       Output of mass loss or tidal capture events.
*       --------------------------------------------
*
      INCLUDE 'common6.h'
      DATA  AU  /2.0627E+05/
*
*
*       Check counter for stellar evolution events.
      IF (NMDOT.GT.0) THEN
          DO 5 J = 1,10
              NTYPE(J) = 0
    5     CONTINUE
*
          KM = 1
          DO 10 J = 1,N
              LJ = NAME(J)
              KW = KSTAR(LJ)
              KW = MIN(KW,10)
              KW = MAX(KW,1)
              NTYPE(KW) = NTYPE(KW) + 1
              KM = MAX(KM,KW)
   10     CONTINUE
*
          WRITE (6,15)
   15     FORMAT (/,5X,'NMDOT  NRG  NHE   NRS  NWD  NSN  NTZ  NBS',
     &               '  ZMRG  ZMHE  ZMRS  ZMWD  ZMSN  ZMDOT  NTYPE')
          WRITE (6,20)  NMDOT, NRG, NHE, NRS, NWD, NSN, NTZ, NBS,
     &                  ZMRG, ZMHE, ZMRS, ZMWD, ZMSN, ZMDOT,
     &                  (NTYPE(J),J=1,KM)
   20     FORMAT (' #4',I7,2I5,I6,4I5,5F6.0,F7.0,I7,9I4)
      END IF
*
*       Check output for tidal capture or collisions.
      IF (NDISS + NCOLL.GT.0) THEN
          WRITE (6,25)
   25     FORMAT (/,5X,'NDISS  NTIDE  NSYNC  NCOLL   EBIN  ESYNC',
     &               '  ECOLL')
          WRITE (6,30)  NDISS, NTIDE, NSYNC, NCOLL,  EBIN, ESYNC, ECOLL
   30     FORMAT (' #5',4I7,3F7.3)
      END IF
*
      RETURN
*
      END
