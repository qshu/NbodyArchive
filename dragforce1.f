       SUBROUTINE DRAGFORCE1(I)
*
*
*       Drag force & first derivative.
*       -------------------------
*
*              Types of BODY(I),X(K,I),XDOT(K,I) is Defined in "comon6.h" 
      INCLUDE 'common6.h'
      REAL*8  FIRR(3),FD(3),FDRAG(3),FDDRAG(3),
     &        R2DISC,R2DENS,RDISC,DDENS,DDOTDENS,VXRE,VYRE,VZRE,VRE,
     &        CKFDRAG,DVXRE,DVYRE,DVZRE,DVRE,KFDOT,ECCENTR1(NMAX),
     &        AORBIT(NMAX),BORBIT(NMAX),CORBIT(NMAX),INCLIN1(NMAX)
*
      LOGICAL LPR
      LPR=NAME(I).EQ.250
*
*              Auxiliary Calculations
      R2DISC=X(1,I)*X(1,I)+X(2,I)*X(2,I)
      R2DENS=R2DISC+0.01
      RDISC=SQRT(R2DISC)
*
      DDENS=R2DENS**(-3.0/8.0)*EXP(-X(3,I)*X(3,I)/(2.0*HZ))
      DDOTDENS=(-3.0/4.0)*(1.0/R2DENS)*
     &         (X(1,I)*XDOT(1,I)+X(2,I)*XDOT(2,I))-
     &         (1.0/HZ)*X(3,I)*XDOT(3,I)
*
      VXDISC=-C3*X(2,I)/((R2DISC)**(3.0/4.0))
      VYDISC=C3*X(1,I)/((R2DISC)**(3.0/4.0))
      VXRE=VXDISC-XDOT(1,I)
      VYRE=VYDISC-XDOT(2,I)
      VZRE=-XDOT(3,I)
      VRE=SQRT(VXRE*VXRE + VYRE*VYRE + VZRE*VZRE)
      CKFDRAG=Q_DRAG*DDENS*VRE
*
*              Mean of Drag Force
      FDRAG(1)=CKFDRAG*VXRE
      FDRAG(2)=CKFDRAG*VYRE
      FDRAG(3)=CKFDRAG*VZRE
*      	 
      EDISS1=EDISS1+BODY(I)*(FDRAG(1)*XDOT(1,I)+FDRAG(2)*XDOT(2,I)+
     &       FDRAG(3)*XDOT(3,I))*STEP(I) 
      
*       Turning-Points of Orbits
      IF(LPR) THEN
      V1T1(I)=X(1,I)*XDOT(1,I)+X(2,I)*XDOT(2,I)+X(3,I)*XDOT(3,I)
      V2T2(I)=XOLD(1,I)*VOLD(1,I)+XOLD(2,I)*VOLD(2,I)+
     &        XOLD(3,I)*VOLD(3,I)
      V12T(I)=V1T1(I)*V2T2(I)
*
      IF(V12T(I).LE.0.0 .AND. V2T2(I).GE.0.0) THEN
      RSTARMAX(I)=SQRT(R2DISC+X(3,I)*X(3,I))
      COUNTER1(I)=COUNTER1(I)+1
      XOLD1(1,I)=X(1,I)
      XOLD1(2,I)=X(2,I)
      XOLD1(3,I)=X(3,I)
      END IF
*      
      IF(V12T(I).LE.0.0 .AND. V2T2(I).LT.0.0) THEN   
      RSTARMIN(I)=SQRT(R2DISC+X(3,I)*X(3,I))
      COUNTER1(I)=COUNTER1(I)+1
      XOLD2(1,I)=X(1,I)		
      XOLD2(2,I)=X(2,I)
      XOLD2(3,I)=X(3,I)
      END IF
*
      IF(COUNTER1(I).EQ.1) THEN
      IF(XOLD(3,I)*X(3,I).LE.0.0) GO TO 516
      IF(XOLD(2,I)*X(2,I).LE.0.0) GO TO 516
      IF(XOLD(1,I)*X(1,I).LE.0.0) GO TO 516
      GO TO 515
 516  XINTER(1,I)=X(1,I)
      XINTER(2,I)=X(2,I)
      XINTER(3,I)=X(3,I)
      END IF
*
 515  CONTINUE
      XOLD(1,I)=X(1,I)
      XOLD(2,I)=X(2,I)
      XOLD(3,I)=X(3,I)
      VOLD(1,I)=XDOT(1,I)
      VOLD(2,I)=XDOT(2,I)
      VOLD(3,I)=XDOT(3,I)
      END IF
*
*      Excenricity & Inclination of Orbits
      IF(COUNTER1(I).EQ.2) THEN
      ECCENTR1(I)=(RSTARMAX(I)-RSTARMIN(I))/(RSTARMAX(I)+RSTARMIN(I))
*     PRINT*,'ECCENTRICITY=',ECCENTR1(I)
      AORBIT(I)=XINTER(2,I)*XOLD1(3,I)+XINTER(3,I)*XOLD2(2,I)+
     &           XOLD1(2,I)*XOLD2(3,I)-XOLD1(3,I)*XOLD2(2,I)-
     &           XOLD1(2,I)*XINTER(3,I)-XOLD2(3,I)*XINTER(2,I)
      BORBIT(I)=-(XINTER(1,I)*XOLD1(3,I)+XINTER(3,I)*XOLD2(1,I)+
     &           XOLD1(1,I)*XOLD2(3,I)-XOLD1(3,I)*XOLD2(1,I)-
     &           XOLD1(1,I)*XINTER(3,I)-XOLD2(3,I)*XINTER(1,I))
      CORBIT(I)=XINTER(1,I)*XOLD1(2,I)+XINTER(2,I)*XOLD2(1,I)+
     &          XOLD1(1,I)*XOLD2(2,I)-XOLD1(2,I)*XOLD2(1,I)-
     &          XOLD2(2,I)*XINTER(1,I)-XOLD1(1,I)*XINTER(2,I)
*      
      IF(VXDISC*XDOT(1,I).GT.0.0 .AND. VYDISC*X(2,I).GT.0.0) THEN
      INCLIN1(I)=(360/TWOPI)*ACOS(ABS(CORBIT(I))/SQRT(AORBIT(I)**2+
     &            BORBIT(I)**2+CORBIT(I)**2))
      ELSE
      INCLIN1(I)=180-(360/TWOPI)*ACOS(ABS(CORBIT(I))/
     &           SQRT(AORBIT(I)**2+BORBIT(I)**2+CORBIT(I)**2))
      END IF
*
      COUNTER1(I)=O
*     PRINT*,'INCLINATION=',INCLIN1(I)
      END IF
*
*     IF(LPR) THEN
*     WRITE(51,511) TIME,(FDRAG(K), FDDRAG(K),X(K,I),XDOT(K,I),K=1,3)
*     WRITE(52,512) TIME,INCLIN1(I),ECCENTR1(I)
*511  FORMAT(1X,1P,13(D15.5,1X))
*512  FORMAT(1X,1P,3(D15.5,1X))
*     END IF
*      
*              Auxiliary Calculations
      DVXRE=-(C3**2.0)*X(1,I)/(R2DISC**(3.0/2.0))-
     &       (FI(1,I)+FR(1,I)+FDRAG(1))
      DVYRE=-(C3**2.0)*X(2,I)/(R2DISC**(3.0/2.0))-
     &       (FI(2,I)+FR(2,I)+FDRAG(2))
      DVZRE=-(FI(3,I)+FR(3,I)+FDRAG(3))
      DVRE=(VXRE*DVXRE+VYRE*DVYRE+VZRE*DVZRE)/(VRE*VRE)
      KFDOT=Q_DRAG*DDENS*VRE
*   
*              Mean of Drag Force Derivative
      FDDRAG(1)=KFDOT*(DDOTDENS*VXRE+DVRE*VXRE+DVXRE)
      FDDRAG(2)=KFDOT*(DDOTDENS*VYRE+DVRE*VYRE+DVYRE)
      FDDRAG(3)=KFDOT*(DDOTDENS*VZRE+DVRE*VZRE+DVZRE)
*     PRINT*,FDDRAG(1),FDDRAG(2),FDDRAG(3)
*   
*     END IF
*    
*              Total Force Acting on a Star
      DO K=1,3
         FI(K,I)=FI(K,I)+FDRAG(K)
         D1(K,I)=D1(K,I)+FDDRAG(K)
      END DO 
*
      RETURN
*
      END
      
