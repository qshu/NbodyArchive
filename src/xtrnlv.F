      SUBROUTINE XTRNLV(I1,I2)
!VAMPINST include
*
*
*       External potential.
*       -------------------
*
      INCLUDE 'common6.h'
*
*
*       See whether to include the galactic tidal force.
!VAMPINST subroutine_start
      ET = 0.0D0
      IF (TIDAL(1).GT.0.0) THEN
          DO 10 I = I1,I2
              ET = ET - 0.5D0*BODY(I)*(TIDAL(1)*X(1,I)**2 +
     &                                 TIDAL(3)*X(3,I)**2)
   10     CONTINUE
      END IF
*
*       Place sum in ETIDE and single particle contribution in HT.
      IF (I2.GT.I1) THEN
          ETIDE = ET
      ELSE
          HT = ET
      END IF
*
!VAMPINST return
          RETURN
      RETURN
*
!VAMPINST subroutine_end
      END
