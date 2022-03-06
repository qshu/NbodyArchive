        SUBROUTINE DRAGBLCKHL1(I)
*
*
*       Drag force & first derivative.
*       -------------------------
*
*              Types of BODY(I),X(K,I),XDOT(K,I) is Defined in "comon6.h" 
      INCLUDE 'common4.h'
      REAL*8  AA1(9),FBLACKH(3),FDBLACKH(3)
*
*            Auxiliary Calculations for Black Hole Interaction
      DO K=1,3
      AA1(K)=-X(K,I)
      AA1(K+3)=-XDOT(K,I)
      END DO
      AA1(7)=1.0/(AA1(1)*AA1(1)+AA1(2)*AA1(2)+AA1(3)*AA1(3)+EPS1*EPS1)
      AA1(8)=CMBLHOLE*AA1(7)*SQRT(AA1(7))
      AA1(9)=3.0*(AA1(1)*AA1(4)+AA1(2)*AA1(5)+AA1(3)*AA1(6))*AA1(7)      
*
      DO K=1,3
      FBLACKH(K)=AA1(K)*AA1(8)
      FDBLACKH(K)=(AA1(K+3)-AA1(K)*AA1(9))*AA1(8)
      END DO
*
*              Total Force Acting on a Star
*
       DO K=1,3
         F(K,I)=F(K,I)+FBLACKH(K)/2.D0
         FDOT(K,I)=FDOT(K,I)+FDBLACKH(K)/6.D0
       END DO 
*
      RETURN
*
      END
      
