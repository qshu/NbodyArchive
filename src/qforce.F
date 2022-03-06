      SUBROUTINE QFORCE(Q,F,QF)
!VAMPINST include
*
*
*       KS transpose of Q times 2F.
*       ---------------------------
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Q(4),F(3),QF(4)
*
*
!VAMPINST subroutine_start
      QF(1)=2.D0*(+Q(1)*F(1)+Q(2)*F(2)+Q(3)*F(3))
      QF(2)=2.D0*(-Q(2)*F(1)+Q(1)*F(2)+Q(4)*F(3))
      QF(3)=2.D0*(-Q(3)*F(1)-Q(4)*F(2)+Q(1)*F(3))
      QF(4)=2.D0*(+Q(4)*F(1)-Q(3)*F(2)+Q(2)*F(3))
*
!VAMPINST return
          RETURN
      RETURN
!VAMPINST subroutine_end
      END
