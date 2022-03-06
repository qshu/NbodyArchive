        SUBROUTINE DRAGSTARPCM(I,FC,FCD)
*
*
*       Pertubation on C.M. due to Central Force
*       ----------------------------------------
*
*     Types of BODY(I),X(K,I),XDOT(K,I) is Defined in "comon6.h"
      INCLUDE 'common6.h'
      INTEGER I,I1,I2
      REAL*8  FC(3),FCD(3)
      REAL*8  FCM(3),FCMD(3),F1(3),F2(3),F1D(3),F2D(3)
*
*       KS Pair Components (resolved in intgrt)
      I2 = 2*(I-N)
      I1 = I2-1
      DO K=1,3
       FCM(K) = 0.0
       FCMD(K) = 0.0
       F1(K) = 0.0
       F1D(K) = 0.0
       F2(K) = 0.0
       F2D(K) = 0.0
      END DO
      CALL MASSCENTROBJECT(I,X(1,I),XDOT(1,I),FCM,FCMD)
      CALL DRAGFORCE(I,X(1,I),XDOT(1,I),FCM,FCMD,FR,2)
      CALL MASSCENTROBJECT(I1,X(1,I1),XDOT(1,I1),F1,F1D)
      CALL DRAGFORCE(I1,X(1,I1),XDOT(1,I1),F1,F1D,FR,2)
      CALL MASSCENTROBJECT(I2,X(1,I2),XDOT(1,I2),F2,F2D)
      CALL DRAGFORCE(I2,X(1,I2),XDOT(1,I2),F2,F2D,FR,2)
*     
      DO K=1,3
       FC(K) =FC(K) -FCM(K) +(BODY(I1)*F1(K) +BODY(I2)*F2(K) )/BODY(I)
       FCD(K)=FCD(K)-FCMD(K)+(BODY(I1)*F1D(K)+BODY(I2)*F2D(K))/BODY(I)
      END DO
*
      RETURN
*
      END
