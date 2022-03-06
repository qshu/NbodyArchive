      SUBROUTINE EVENTS
*
*
*       Output of mass loss or tidal capture events.
*       --------------------------------------------
*
      INCLUDE 'common6.h'
      INTEGER  NTYPE(17)
      real*8 thookf,tbgbf
      external thookf,tbgbf
*
*
*       Check counter for stellar evolution events.
      IF (NMDOT.GT.0) THEN
          DO 5 J = 1,16
              NTYPE(J) = 0
    5     CONTINUE
*
          KM = 1
          DO 10 J = 1,N
              KW = KSTAR(J) + 1
              KW = MIN(KW,16)
              KW = MAX(KW,1)
              NTYPE(KW) = NTYPE(KW) + 1
              KM = MAX(KM,KW)
   10     CONTINUE
*
          WRITE (6,15)
   15     FORMAT (/,5X,'NMDOT  NRG  NHE   NRS  NNH  NWD  NSN  NBH  NBS',
     &               '  ZMRG  ZMHE  ZMRS  ZMNH  ZMWD  ZMSN  ZMDOT',
     &               '  NTYPE')
          WRITE (6,20)  NMDOT, NRG, NHE, NRS, NNH, NWD, NSN, NBH, NBS,
     &                  ZMRG, ZMHE, ZMRS, ZMNH, ZMWD, ZMSN, ZMDOT,
     &                  (NTYPE(J),J=1,KM)
   20     FORMAT (' #4',I7,2I5,I6,5I5,6F6.1,F7.1,I7,I6,14I4)
      END IF
*
*       Determine turnoff mass at current cluster age (cf. routine STAR).
      IF (TIME.LE.0.0D0) THEN
          TURN = BODY1*ZMBAR
      ELSE
          TPHYS = TIME*TSTAR
          TURN = BODY1*ZMBAR
          TURN2 = 2.0*TURN
   25     TM = MAX(zpars(8),thookf(turn))*tbgbf(turn)
          IF (TM.GT.TPHYS) THEN
              TURN = 1.01*TURN
          ELSE
              TURN = 0.985*TURN
          END IF
          IF (ABS(TM - TPHYS).GT.1.0.AND.TURN.LT.TURN2) GO TO 25
      END IF
*
*       Check output for tidal capture or collisions.
      IF (NDISS + NCOLL.GT.0) THEN
          EESC = E(5) + E(7)
          WRITE (6,30)
   30     FORMAT (/,5X,'NDISS  NTIDE  NSYNC  NCOLL    EBIN   ECOLL',
     &               '   EMDOT  ECDOT  EKICK   EESC')
          WRITE (6,35)  NDISS, NTIDE, NSYNC, NCOLL, EBIN, ECOLL,
     &                  EMDOT, ECDOT, EKICK, EESC
   35     FORMAT (' #5',4I7,3F8.2,3F7.2)
      END IF
*
      RETURN
*
      END
