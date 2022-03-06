      SUBROUTINE DRAGMBHINIT(I1,I2)
*
*
*       Force & first derivative.
*       -------------------------
*
      INCLUDE 'common4.h'
*     Call Drag Force
      DO 90 I=I1,I2
      CALL DRAGBLCKHL1(I)
      CALL DRAGFORCE1(I)
 90   CONTINUE
*
      RETURN     
      END
