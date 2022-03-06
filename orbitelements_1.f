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
     &        INCLDISTR10(NMAX),EXCEDISTR20(NMAX)       
*       
      INTEGER STARNUMBER, IANGLE, K111, K1111, NOFPARTCL
      SAVE NOFPARTCL
      DATA NOFPARTCL /0/
      LOGICAL XCOROT(NMAX),COUNTER1(NMAX)
      LOGICAL P123(NMAX),P1_FIXED(NMAX),P2_FIXED(NMAX),
     &        P3_FIXED(NMAX)
*
      DO K=1,3
      XI(K)=X(K,I)
      XIDOT(K)=XDOT(K,I)
      END DO
*      
*       Recalling the Name of the Star
      J=NAME(I)
      IF(J.GT.N) THEN
      GO TO 556
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
         IF(INCLIN_L(J).LE.90.0) THEN
         XCOROT(J)=.TRUE.
         ELSEIF(INCLIN_L(J).GT.90.0) THEN
         XCOROT(J)=.FALSE.
         END IF
      COUNTER1(J)=.TRUE.
      GOTO 555
      END IF
*
*
      RTEK=SQRT(XI(1)**2.0+XI(2)**2.0+XI(3)**2.0)
*
 551  CONTINUE
*
      IF(P123(J).AND.COUNTER1(J)) THEN 
      P123(J)=.FALSE.
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
      IF(V23(J).LT.0.0.AND.V2(J).GT.0.0.AND.TTOT.GT.19.0
     &.OR.V12(J).LT.0.0.AND.V1(J).GT.0.0.AND.TTOT.GT.19.0) THEN
*
      NOFPARTCL=NOFPARTCL+1
      PRINT*,'N,NUMBER(I):',NAME(I),NTOT
      PRINT*,'INCLIN BY ANG.MOM  :',INCLIN_L(J)
      PRINT*,'INCLIN BY 3 POINTS :',INCLIN(J),NOFPARTCL
      PRINT*,'Eccentricity & A-Axis angl        :',ECC_3P(J)
      PRINT*,'****************************************'
      COUNTER1(J)=.FALSE.
*      END IF
*      IF(XCOROT(J)) THEN
      INCLDISTR1(J)=INCLIN(J)
      EXCEDISTR2(J)=ECC_3P(J)
*      ELSE
      INCLDISTR10(J)=INCLIN(J)
      EXCEDISTR20(J)=ECC_3P(J)
      END IF
*
*      END IF
*
      END IF
*
      IF(NOFPARTCL.EQ.4950) THEN
      NOFPARTCL=0
      PRINT*,'It is done!'
*
      DO IANGLE=0,180,1
      STARNUMBER=0
      DO I1=1,N
      IF(INCLDISTR1(I1).GT.IANGLE.AND.INCLDISTR1(I1).LE.IANGLE+1) THEN
      STARNUMBER=STARNUMBER+1
      END IF
      END DO
      WRITE(18,524) IANGLE,STARNUMBER
 524  FORMAT(2(I5))
      END DO
*
      DO CEXC=0.0,3.99,0.01
      STARNUMBER=0
      DO I1=1,N
      IF(EXCEDISTR2(I1).GT.CEXC.AND.EXCEDISTR2(I1).LE.CEXC+0.01) THEN
      STARNUMBER=STARNUMBER+1
      END IF
      END DO
      WRITE(39,531) CEXC,STARNUMBER
 531  FORMAT(D15.5,I5)
      END DO
*
      END IF
*
*     END IF
*
 555  ITHIRD_POINT(J)=1.0
      IF(P1_FIXED(J)) THEN
      RTEK1=SQRT(P_1(1,J)**2.0+P_2(2,J)**2.0+P_3(3,J)**2.0)	
      ANGL_P1P2=(360.0/TWOPI)*
     &       ACOS((P_1(1,J)*XI(1)+P_1(2,J)*XI(2)+P_1(3,J)*XI(3))/
     &     (RTEK*SQRT(P_1(1,J)**2.0+P_1(2,J)**2.0+P_1(3,J)**2.0)))
      END IF
* 
      IF(ANGL_P1P2.GE.20.0D0.AND.P2_FIXED(J)) THEN
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
      IF(ANGL_P2P3.GE.20.0D0.AND.P3_FIXED(J)) THEN
      ANGLE_P2P3=0.0
      P123(J)=.TRUE.
      END IF	
*
 556  RETURN
*
*
      END
      





























