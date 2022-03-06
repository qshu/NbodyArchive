      SUBROUTINE DRAGMCO_INIT(I1,I2,KCASE)
*
*
*       Force & first derivative.
*       -------------------------
*
      INCLUDE 'common6.h'
*     Call Drag Force & Massive Central Object Force
      DO 90 I=I1,I2
      CALL MASSCENTROBJECT1(I)
      CALL DRAGFORCE1(I)
*      DO 95 K = 1,3    - This loop is done fploly1.f
*
*        F(K,I)=FI(K,I)+FR(K,I)
*        FDOT(K,I)=D1(K,I)+D1R(K,I)
*
*        D0(K,I) = FI(K,I)
*        D0R(K,I) = FR(K,I)
*        FIDOT(K,I)=D1(K,I)
*        FRDOT(K,I)=D1R(K,I)
*   95 CONTINUE
   90 CONTINUE
*
      RETURN     
      END
