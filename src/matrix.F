      SUBROUTINE MATRIX(U,A)
!VAMPINST include
*
*
*       Levi-Civita matrix.
*       -------------------
*
      REAL*8  U(4),A(3,4)
*
*
*       Set current transformation matrix.
!VAMPINST subroutine_start
      A(1,1) =  U(1)
      A(1,2) = -U(2)
      A(1,3) = -U(3)
      A(1,4) =  U(4)
      A(2,1) =  U(2)
      A(2,2) =  U(1)
      A(2,3) = -U(4)
      A(2,4) = -U(3)
      A(3,1) =  U(3)
      A(3,2) =  U(4)
      A(3,3) =  U(1)
      A(3,4) =  U(2)
*
!VAMPINST return
          RETURN
      RETURN
*
!VAMPINST subroutine_end
      END
