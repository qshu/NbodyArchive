      SUBROUTINE PERMIT(RPERT,IGO)
*
*
*       Permission for unperturbed triple or quad.
*       ------------------------------------------
*
      INCLUDE 'common4.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
*
*
*       Search any existing subsystem.
      ISUB = 0
      DO 10 L = 1,NSUB
*       Distinguish between triple & quad case.
          IF (JCOMP.LE.N.AND.NAMES(4,L).EQ.0) THEN
              ISUB = L
      ELSE IF (JCOMP.GT.N.AND.NAMES(4,L).GT.0) THEN
              ISUB = L
          END IF
   10 CONTINUE
*
*       Do not allow a second regularization of the same type.
      IF (ISUB.GT.0) THEN
          IGO = 1
*       Enforce termination at next extension if new system is smaller.
          IF (RPERT.LT.RMAXS(ISUB)) THEN
              STEPS(ISUB) = 0.0
          END IF
      END IF
*
      RETURN
*
      END
