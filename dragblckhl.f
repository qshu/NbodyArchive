        SUBROUTINE DRAGBLCKHL(I,FIRR,FD,XI,VI)
*
*
*       Drag force & first derivative.
*       -------------------------
*
*              Types of BODY(I),X(K,I),XDOT(K,I) is Defined in "comon6.h" 
      INCLUDE 'common6.h'
*     COMMON /BLACKFORCE/FDBLACKH(3)
      REAL*8  FIRR(3),FD(3),XI(3),VI(3),AA1(9),FBLACKH(3),FDBLACKH(3)
*
*            Auxiliary Calculations for Black Hole Interaction
      DO K=1,3
      AA1(K)=-XI(K)
      AA1(K+3)=-VI(K)
      END DO
      AA1(7)=1.0/(AA1(1)*AA1(1)+AA1(2)*AA1(2)+AA1(3)*AA1(3)+EPS1*EPS1)
      AA1(8)=CMBLHOLE*AA1(7)*SQRT(AA1(7))
      AA1(9)=3.0*(AA1(1)*AA1(4)+AA1(2)*AA1(5)+AA1(3)*AA1(6))*AA1(7)      
*
      DO K=1,3
      FBLACKH(K)=AA1(K)*AA1(8)
      FDBLACKH(K)=(AA1(K+3)-AA1(K)*AA1(9))*AA1(8)
      END DO
*      EBLCKHL=EBLCKHL-BODY(I)*(FBLACKH(1)*XDOT(1,I)+
*     &        FBLACKH(2)*XDOT(2,I)+ FBLACKH(3)*XDOT(3,I))*STEP(I)
*              
*              Total Force Acting on a Star
*
       DO K=1,3
         FIRR(K)=FIRR(K)+FBLACKH(K)
         FD(K)=FD(K)+FDBLACKH(K)
       END DO 
*
      RETURN
*
      END
      
