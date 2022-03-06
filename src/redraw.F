      SUBROUTINE REDRAW(IC,LMIN,LMAX,IJ,LI,SUC,LOOP)
!VAMPINST include
*
*
*       Check chain connections.
*       -----------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      INTEGER IC(1),ICC(2),IJ(NMXM,2)
      LOGICAL SUC,LOOP
*
*
!VAMPINST subroutine_start
      SUC=.FALSE.
      LOOP=.FALSE.
      ICC(1)=IC(LMIN)
      ICC(2)=IC(LMAX)
      DO I=1,2
      DO J=1,2
      IF(ICC(I).EQ.IJ(LI,J))THEN
      JC=3-J
      LOOP=.TRUE.
      DO L=LMIN,LMAX
      IF(IC(L).EQ.IJ(LI,JC)) Then
!VAMPINST return
          RETURN
      RETURN
      END IF
      END DO
      SUC=.TRUE.
      LOOP=.FALSE.
      IF(I.EQ.1)THEN
      LMIN=LMIN-1
      IC(LMIN)=IJ(LI,JC)
!VAMPINST return
          RETURN
      RETURN
      ELSE
      LMAX=LMAX+1
      IC(LMAX)=IJ(LI,JC)
!VAMPINST return
          RETURN
      RETURN
      END IF
      END IF
      END DO
      END DO
*
!VAMPINST return
          RETURN
      RETURN
!VAMPINST subroutine_end
      END
