       SUBROUTINE DRAGFORCE(I,FIRR,FD,XI,VI,KDR)
*
*
*       Drag force & first derivative.
*       -------------------------
*
*              Types of BODY(I),X(K,I),XDOT(K,I) is Defined in "comon6.h" 
      INCLUDE 'common6.h'
      REAL*8  FIRR(3),FD(3),FDRAG(3),FDDRAG(3),
     &        R2DISC,R2DENS,RDISC,DDENS,DDOTDENS,VXRE,VYRE,VZRE,VRE,
     &        CKFDRAG,DVXRE,DVYRE,DVZRE,DVRE,KFDOT,XI(3),VI(3)
*
*
      DO K=1,3
      X(K,I)=XI(K)
      XDOT(K,I)=VI(K)
      END DO
*               Auxiliary Calculations
      R2DISC=X(1,I)*X(1,I)+X(2,I)*X(2,I)
      R2DENS=R2DISC+0.01
      RDISC=SQRT(R2DISC)
*
      DDENS=R2DENS**(-3.0/8.0)*EXP(-X(3,I)*X(3,I)/(2.0*HZ))
      DDOTDENS=(-3.0/4.0)*(1.0/R2DENS)*
     &         (X(1,I)*XDOT(1,I)+X(2,I)*XDOT(2,I))-
     &         (1.0/HZ)*X(3,I)*XDOT(3,I)
*
      VXDISC=-C3*X(2,I)/((R2DENS)**(3.0/4.0))
      VYDISC=C3*X(1,I)/((R2DENS)**(3.0/4.0))
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
      IF(KDR.EQ.0)
     &EDISS1=EDISS1+BODY(I)*(FDRAG(1)*XDOT(1,I)+FDRAG(2)*XDOT(2,I)+
     &       FDRAG(3)*XDOT(3,I))*STEP(I)      	  
*
*              Auxiliary Calculations
      DTR = TIME - T0R(I)
      DVXRE=-(C3**2.0)*X(1,I)/(R2DENS**(3.0/2.0))-
     &       (FIRR(1)+FR(1,I)+FRDOT(1,I)*DTR+FDRAG(1))
      DVYRE=-(C3**2.0)*X(2,I)/(R2DENS**(3.0/2.0))-
     &       (FIRR(2)+FR(2,I)+FRDOT(2,I)*DTR+FDRAG(2))
      DVZRE=-(FIRR(3)+FR(3,I)+FRDOT(3,I)*DTR+FDRAG(3))
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
      fdx = dsqrt(fdrag(1)**2+fdrag(2)**2+fdrag(3)**2)
      fddx = dsqrt(fddrag(1)**2+fddrag(2)**2+fddrag(3)**2)
      stepdd = fdx/fddx
*     if(stepdd.lt.step(i))print*,' i,dt,dtd=',i,step(i),stepdd
*              Total Force Acting on a Star
      DO K=1,3
         FIRR(K)=FIRR(K)+FDRAG(K)
         FD(K)=FD(K)+FDDRAG(K)
      END DO 
*
      RETURN
*
      END
