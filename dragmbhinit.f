      SUBROUTINE DRAGMBHINIT
*
*
*       Force & first derivative.
*       -------------------------
*
      INCLUDE 'common6.h'
*     Call Drag Force
      DO 90 I=1,N
      CALL DRAGBLCKHL1(I)
      CALL DRAGFORCE1(I)
      DO 95 K = 1,3
*
        F(K,I)=FI(K,I)+FR(K,I)
        FDOT(K,I)=D1(K,I)+D1R(K,I)
*
        D0(K,I) = FI(K,I)
        D0R(K,I) = FR(K,I)
        FIDOT(K,I)=D1(K,I)
        FRDOT(K,I)=D1R(K,I)
   95 CONTINUE
   90 CONTINUE
*
      RETURN     
      END
