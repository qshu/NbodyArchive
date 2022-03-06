       SUBROUTINE DRAGFORCE(I,XI,XIDOT,FIRR,FD,FREG,ICASE)
*
*
*       Drag force & first derivative.
*       -------------------------
*
*              Types of BODY(I),X(K,I),XDOT(K,I) is Defined in "comon6.h" 
      INCLUDE 'common6.h'
      REAL*8  FIRR(3),FD(3),FDRAG(3),FDDRAG(3),FTOT(3),FREG(3),
     &        R2DISC,R2DENS,RDISC,DDENS,DDOTDENS,VXRE,VYRE,VZRE,VRE,
     &        CKFDRAG,DVXRE,DVYRE,DVZRE,DVRE,KFDOT,XI(3),XIDOT(3),
     &        R2DOT,MPLI
*
*        Calling Parameter: ICASE=0 from nbint.f
*                                =1 from regint.f
*                                =2 from kspert.f, fpert.f
*
      DT = TIME - T0(I)
      DTR = TIME - T0R(I)
*
      DO K=1,3
         IF(ICASE.EQ.0) FTOT(K) = FIRR(K) + FR(K,I) + DTR*FRDOT(K,I)
         IF(ICASE.EQ.1) FTOT(K) = FIRR(K) + FREG(K)
         IF(ICASE.GT.1) FTOT(K) = FI(K,I) + DT*FIDOT(K,I) +
     &               FR(K,I) + DTR*FRDOT(K,I)
*        It is possible that higher precision FTOT using D2,D3,D2R,D3R
*         may be needed in the future.
      END DO
*
*               Auxiliary Calculations
      R2DISC=XI(1)*XI(1)+XI(2)*XI(2)
      R2DENS=R2DISC+0.01
      RDISC=SQRT(R2DISC)
*
*
*               Mass of the cluster according to Plummer distr. inside of R(I)
      C3=0.1
      IF(RDISC.LT.1.0E-03.OR.RDISC.GT.0.8) GOTO 20
      MPLI=(RDISC/(1.5*TWOPI/16.0))**3.0/
     &     (1+(RDISC/(1.5*TWOPI/16.0))**2.0)**1.5
      C3=SQRT(C3+MPLI)
 20   CONTINUE
*
*
      DDENS=R2DENS**(-3.0/8.0)*EXP(-XI(3)*XI(3)/(2.0*HZ))
      DDOTDENS=(-3.0/4.0)*(1.0/R2DENS)*
     &         (XI(1)*XIDOT(1)+XI(2)*XIDOT(2))-
     &         (1.0/HZ)*XI(3)*XIDOT(3)
*
      VXDISC=-C3*XI(2)/((R2DISC)**(3.0/4.0))
      VYDISC=C3*XI(1)/((R2DISC)**(3.0/4.0))
      VXRE=VXDISC-XIDOT(1)
      VYRE=VYDISC-XIDOT(2)
      VZRE=-XIDOT(3)
      VRE=SQRT(VXRE*VXRE + VYRE*VYRE + VZRE*VZRE)
      CKFDRAG=Q_DRAG*DDENS*VRE
*
*              Mean of Drag Force
      FDRAG(1)=CKFDRAG*VXRE
      FDRAG(2)=CKFDRAG*VYRE
      FDRAG(3)=CKFDRAG*VZRE  
* 
      IF(ICASE.EQ.0) EDISS1_LOC=EDISS1_LOC+BODY(I)*(FDRAG(1)*XIDOT(1)+
     &FDRAG(2)*XIDOT(2)+FDRAG(3)*XIDOT(3))*STEP(I)
*      PRINT*,'DRAG FORCE: EDISS1_LOC & rank',rank,EDISS1_LOC
*
*              Auxiliary Calculations
      R2DOT = 2.0*(XI(1)*XIDOT(1)+XI(2)*XIDOT(2))
      DVXRE = -C3*XIDOT(2)*R2DENS**(-0.75)
     &        +C3*0.75*XI(2)*R2DOT*R2DENS**(-1.75)-(FTOT(1)+FDRAG(1))
      DVYRE = +C3*XIDOT(1)*R2DENS**(-0.75)
     &        -C3*0.75*XI(1)*R2DOT*R2DENS**(-1.75)-(FTOT(2)+FDRAG(2))
      DVZRE = -(FTOT(3)+FDRAG(3))
*
*      DVXRE=-(C3**2.0)*XI(1)/(R2DENS**(3.0/2.0))-
*     &       (FTOT(1)+FDRAG(1))
*      DVYRE=-(C3**2.0)*XI(2)/(R2DENS**(3.0/2.0))-
*     &       (FTOT(2)+FDRAG(2))
*      DVZRE=-(FTOT(3)+FDRAG(3))
*      DVRE=(VXRE*DVXRE+VYRE*DVYRE+VZRE*DVZRE)/(VRE*VRE)
      DVRE=(VXRE*DVXRE+VYRE*DVYRE+VZRE*DVZRE)/VRE
      KFDOT=Q_DRAG*DDENS*VRE
*   
*              Mean of Drag Force Derivative
      FDDRAG(1)=KFDOT*(DDOTDENS*VXRE+DVRE*VXRE/VRE+DVXRE)
      FDDRAG(2)=KFDOT*(DDOTDENS*VYRE+DVRE*VYRE/VRE+DVYRE)
      FDDRAG(3)=KFDOT*(DDOTDENS*VZRE+DVRE*VZRE/VRE+DVZRE)
*     PRINT*,FDDRAG(1),FDDRAG(2),FDDRAG(3)
*   
*     END IF
*    
      fdx = dsqrt(fdrag(1)**2+fdrag(2)**2+fdrag(3)**2)
      fddx = dsqrt(fddrag(1)**2+fddrag(2)**2+fddrag(3)**2)
      stepdd = fdx/fddx
*     if(stepdd.lt.step(i))print*,' i,dt,dtd=',i,step(i),stepdd
*      if(name(i).eq.191)print*,' !!!!!!!i,dt,dtd=',fdx,step(i),stepdd
*              Total Force Acting on a Star
      DO K=1,3
         FIRR(K)=FIRR(K)+FDRAG(K)
         FD(K)=FD(K)+FDDRAG(K)
      END DO 
*
      RETURN
*
      END


