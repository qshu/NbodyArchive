      SUBROUTINE PERMIT(RPERT,IGO)
!VAMPINST include
*
*
*       Permission for unperturbed triple or quad.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
*
*
*       Search any existing subsystems.
!VAMPINST subroutine_start
      ISUB = 0
      DO 10 L = 1,NSUB
*       Distinguish between triple, quad & chain case.
          IF (JCOMP.LE.N.AND.ISYS(L).EQ.1) THEN
              ISUB = L
          ELSE IF (JCOMP.GT.N.AND.ISYS(L).GE.2) THEN
              ISUB = L
          END IF
   10 CONTINUE
*
*       Do not allow a second regularization of the same type.
      IF (ISUB.GT.0) THEN
          IGO = 1
*       Enforce termination at next extension if new system is smaller.
          IF (RPERT.LT.RMAXS(ISUB)) THEN
              STEPS(ISUB) = 0.0D0
          END IF
      END IF
*
!VAMPINST return
          RETURN
      RETURN
*
!VAMPINST subroutine_end
      END
