      SUBROUTINE ORBITELEMENTS(I)
*
*
*       Drag force & first derivative.
*       -------------------------
*
      INCLUDE 'common6.h'
      REAL*8  ECC_3P(NMAX),
     &        INCLIN(NMAX),ANGLE1,INCLIN_L(NMAX),
     &        EXCENTR2(NMAX),CENER,STARNUMBER1,CEXC,
     &        RP_1(NMAX),RP_2(NMAX),ANGL_P1P2,
     &        ANGL_P2P3,ANGLE_C,V1(NMAX),V2(NMAX),V3(NMAX),
     &        V12(NMAX),V23(NMAX),
     &        RP_3(NMAX),CONST1,CONST3,
     &        P_ORBIT,PR1,PR3,R1R3,RTEK1,
     &        ANGLE2(NMAX),RTEK,XI(3),XIDOT(3),
     &        INCLDISTR1(NMAX),EXCEDISTR2(NMAX),
     &        INCLDISTR10(NMAX),EXCEDISTR20(NMAX),RTEK_C,
     &        A4,RAD11(NMAX),POT11,EBIND1,ETOTLPR, UGOL 
     &        EBIND11     
      INTEGER STARNUMBER, IANGLE, K111, K1111, NOFPARTCL
      SAVE NOFPARTCL
      DATA NOFPARTCL /0/
      LOGICAL XCOROT(NMAX),COUNTER1(NMAX)
      LOGICAL P123(NMAX),P1_FIXED(NMAX),P2_FIXED(NMAX),
     &        P3_FIXED(NMAX),LPR,E_BIND(NMAX)
      COMMON/ENERG/ E_BIND
*
      DO K=1,3
      XI(K)=X(K,I)
      XIDOT(K)=XDOT(K,I)
      END DO
*      
*      A4=0.1/DSQRT(XI(1)*XI(1)+XI(2)*XI(2)+XI(3)*XI(3))
*      POT11=-A4
*      EBIND1=0.5D0*(XIDOT(1)**2.0+XIDOT(2)*2.0
*     &                +XIDOT(3)*2.0)+POT11
*      
      LPR=NAME(I).EQ.1621
*     IF(LPR) THEN
*     IF(LPR) GOTO 556
*     IF(LPR) CALL ENERGY
      EBIND1=BODY(I)*(0.5D0*(XDOT(1,I)**2.0+XDOT(2,I)**2.0+ 
     &                       XDOT(3,I)**2.0)- 
     &       SMBH/DSQRT(X(1,I)**2.0+X(2,I)**2.0+X(3,I)**2.0))
      EBIND11=BODY(I)*0.5D0*(XIDOT(1)**2.0+XIDOT(2)*2.0
     &                +XIDOT(3)*2.0)+PHIDBL(I)
*     END IF
*       Recalling the Name of the Star
      J=NAME(I)
      IF(J.GT.N) THEN
*      GO TO 556
      J=J-N
      END IF
*    
*       Memorizing previous positions
      IF(ITHIRD_POINT(J).EQ.0.0) THEN
      ITHIRD_POINT(J)=1.0
      RTEK=SQRT(XI(1)**2.0+XI(2)**2.0+XI(3)**2.0)
      DO K=1,3
      P_1(K,J)=P_2(K,J)
      P_1DOT(K,J)=P_2DOT(K,J)
      END DO
      P123(J)=.FALSE.
      P1_FIXED(J)=.TRUE.
      P2_FIXED(J)=.TRUE.
      P3_FIXED(J)=.FALSE.
      ANGLE1= (XI(1)*XIDOT(2)-XI(2)*XIDOT(1))/
     &        (SQRT((XI(2)*XIDOT(3)-XI(3)*XIDOT(2))**2.0+
     &              (XI(3)*XIDOT(1)-XI(1)*XIDOT(3))**2.0+
     &              (XI(1)*XIDOT(2)-XI(2)*XIDOT(1))**2.0))
      INCLIN_L(J)=(360/TWOPI)*ACOS(ANGLE1)
         IF(INCLIN_L(J).LE.135.0) THEN
         XCOROT(J)=.TRUE.
         ELSEIF(INCLIN_L(J).GT.135.0) THEN
         XCOROT(J)=.FALSE.
         END IF
      COUNTER1(J)=.TRUE.
      GOTO 555
      END IF
*
*
      RTEK=SQRT(XI(1)**2.0+XI(2)**2.0+XI(3)**2.0)
      RTEK_C=SQRT(XI(1)**2.0+XI(2)**2.0)
*
 551  CONTINUE
*
      IF(P123(J).AND.COUNTER1(J)) THEN 
      DO K=1,3
      P_3(K,J)=XI(K)
      P_3DOT(K,J)=XIDOT(K)
      END DO
*       Laplace Inclination  
      ANGLE1=(P_3(1,J)*XIDOT(2)-P_3(2,J)*XIDOT(1))/
     &       (SQRT((P_3(2,J)*P_3DOT(3,J)-P_3(3,J)*P_3DOT(2,J))**2.0+
     &             (P_3(3,J)*P_3DOT(1,J)-P_3(1,J)*P_3DOT(3,J))**2.0+
     &             (P_3(1,J)*P_3DOT(2,J)-P_3(2,J)*P_3DOT(1,J))**2.0))
      INCLIN_L(J)=(360/TWOPI)*ACOS(ANGLE1)
*
*       Three-point's method eccentricity 'ECC_3P'
*
      CONST1=SQRT((P_2(2,J)*P_3(3,J)-P_2(3,J)*P_3(2,J))**2.0+
     &            (P_2(3,J)*P_3(1,J)-P_2(1,J)*P_3(3,J))**2.0+
     &            (P_2(1,J)*P_3(2,J)-P_2(2,J)*P_3(1,J))**2.0)/
     &       SQRT((P_1(2,J)*P_3(3,J)-P_1(3,J)*P_3(2,J))**2.0+
     &            (P_1(3,J)*P_3(1,J)-P_1(1,J)*P_3(3,J))**2.0+
     &            (P_1(1,J)*P_3(2,J)-P_1(2,J)*P_3(1,J))**2.0)
      CONST3=SQRT((P_2(2,J)*P_1(3,J)-P_2(3,J)*P_1(2,J))**2.0+
     &            (P_2(3,J)*P_1(1,J)-P_2(1,J)*P_1(3,J))**2.0+
     &            (P_2(1,J)*P_1(2,J)-P_2(2,J)*P_1(1,J))**2.0)/
     &       SQRT((P_3(2,J)*P_1(3,J)-P_3(3,J)*P_1(2,J))**2.0+
     &            (P_3(3,J)*P_1(1,J)-P_3(1,J)*P_1(3,J))**2.0+
     &            (P_3(1,J)*P_1(2,J)-P_3(2,J)*P_1(1,J))**2.0)
      RP_1(J)=SQRT(P_1(1,J)**2.0+P_1(2,J)**2.0+P_1(3,J)**2.0) 
      RP_2(J)=SQRT(P_2(1,J)**2.0+P_2(2,J)**2.0+P_2(3,J)**2.0) 
      RP_3(J)=SQRT(P_3(1,J)**2.0+P_3(2,J)**2.0+P_3(3,J)**2.0) 
      P_ORBIT=(CONST1*RP_1(J)+CONST3*RP_3(J)-RP_2(J))/
     &        (CONST1+CONST3-1.0)
      PR1=P_ORBIT-RP_1(J)
      PR3=P_ORBIT-RP_3(J)
      R1R3=SQRT((P_1(2,J)*P_3(3,J)-P_1(3,J)*P_3(2,J))**2.0+
     &          (P_1(3,J)*P_3(1,J)-P_1(1,J)*P_3(3,J))**2.0+
     &          (P_1(1,J)*P_3(2,J)-P_1(2,J)*P_3(1,J))**2.0)
      ECC_3P(J)=(SQRT((PR1*P_3(1,J)-PR3*P_1(1,J))**2.0+
     &                (PR1*P_3(2,J)-PR3*P_1(2,J))**2.0+
     &                (PR1*P_3(3,J)-PR3*P_1(3,J))**2.0))/R1R3
      INCLIN(J)=(360.0/TWOPI)*ACOS((P_2(1,J)*P_3(2,J)-
     &                              P_2(2,J)*P_3(1,J))/
     &        SQRT((P_2(2,J)*P_3(3,J)-P_2(3,J)*P_3(2,J))**2.0+
     &             (P_2(3,J)*P_3(1,J)-P_2(1,J)*P_3(3,J))**2.0+
     &             (P_2(1,J)*P_3(2,J)-P_2(2,J)*P_3(1,J))**2.0))
*
      DO K =1,3 
      P_1(K,J)=P_2(K,J)
      P_2(K,J)=P_3(K,J)
      P_2DOT(K,J)=P_3DOT(K,J)
      END DO
*
      ETOTLPR=EBIND1+EDISS11
      if(rank.eq.0) then
      IF(LPR)THEN
*      IF(E_BIND(J)) THEN
*      GOTO 556
*      ELSE
*
      IF(XI(2).GT.0.0) THEN
      UGOL=ACOS(XI(1)/RTEK_C)
      ELSE
      UGOL=-ACOS(XI(1)/RTEK_C)
      END IF 
*      PRINT*,'N,NUMBER(I),Z:',NAME(I),NTOT,XI(3)
*      PRINT*,'RTEK_XY:',RTEK_C
*      IF(LPR) PRINT*,'INCLIN BY ANG.M  TIME :',INCLIN_L(J),TTOT,EDISS11
*      IF(LPR) PRINT*,'INCLIN BY 3 POINTS EN:',INCLIN(J),NOFPARTCL,EBIND1
*      PRINT*,'Eccentricity & A-Axis angl,ugol:',ECC_3P(J),UGOL
*      PRINT*,'****************************************'
*
      WRITE(29,524) TIME, INCLIN_L(J),EDISS11,EBIND11,RTEK, EBIND1
 524  FORMAT(6(D15.5))
      WRITE(36,613) UGOL,RTEK_C
 613  FORMAT(2(D15.5))
      WRITE(76,614) TIME,XI(1),XI(2),XI(3),XIDOT(1),XIDOT(2),XIDOT(3)
 614  FORMAT(7(D15.5))
      END IF
      end if
      END IF
*
*
 555  ITHIRD_POINT(J)=1.0
      IF(P1_FIXED(J)) THEN
      RTEK1=SQRT(P_1(1,J)**2.0+P_2(2,J)**2.0+P_3(3,J)**2.0)	
      ANGL_P1P2=(360.0/TWOPI)*
     &       ACOS((P_1(1,J)*XI(1)+P_1(2,J)*XI(2)+P_1(3,J)*XI(3))/
     &     (RTEK*SQRT(P_1(1,J)**2.0+P_1(2,J)**2.0+P_1(3,J)**2.0)))
      END IF
* 
      IF(ANGL_P1P2.GE.2.0D0.AND.P2_FIXED(J)) THEN
      DO K=1,3
      P_2(K,J)=XI(K)
      P_2DOT(K,J)=XIDOT(K)
      P_3(K,J)=XI(K)
      END DO
      P1_FIXED(J)=.FALSE.
      P2_FIXED(J)=.FALSE.
      P3_FIXED(J)=.TRUE.
      END IF
*      
*     
      IF(P3_FIXED(J)) THEN
      ANGL_P2P3=(360.0/TWOPI)*
     &ACOS((P_3(1,J)*XI(1)+P_3(2,J)*XI(2)+
     &P_3(3,J)*XI(3))/         
     &(RTEK*SQRT(P_3(1,J)**2.0+P_3(2,J)**2.0+P_3(3,J)**2.0)))      
      V1(J)=P_1(1,J)*P_1DOT(1,J)+P_1(2,J)*P_1DOT(2,J)+
     &   P_1(3,J)*P_1DOT(3,J)
      V2(J)=P_2(1,J)*P_2DOT(1,J)+P_2(2,J)*P_2DOT(2,J)+
     &   P_2(3,J)*P_2DOT(3,J)
      V3(J)=XI(1)*XIDOT(1)+XI(2)*XIDOT(2)+
     &   XI(3)*XIDOT(3)
      V12(J)=V1(J)*V2(J)
      V23(J)=V2(J)*V3(J)
      END IF
*
      IF(ANGL_P2P3.GE.2.0D0.AND.P3_FIXED(J)) THEN
      ANGLE_P2P3=0.0
      P123(J)=.TRUE.
      END IF	
*
 556  RETURN
*
*
      END
      





























