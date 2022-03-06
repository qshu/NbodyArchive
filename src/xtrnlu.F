      SUBROUTINE XTRNLU(X,M,NB,UG)
!VAMPINST include
*
*
*       Harmonic test potential.
*       ------------------------
*
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8 X(1),M(1)
*
*
!VAMPINST subroutine_start
      UG=0.0
      L=0
      DO I=1,NB
      DO K=1,3
      L=L+1
      UG=UG-5.0*X(L)**2*M(I)
      END DO
      END DO
*
!VAMPINST return
          RETURN
      RETURN
!VAMPINST subroutine_end
      END
