       SUBROUTINE DRAGFORCE1(I)
*
*
*       Drag force & first derivative.
*       -------------------------
*
*              Types of BODY(I),X(K,I),XDOT(K,I) is Defined in "comon6.h" 
      INCLUDE 'common6.h'
      REAL*8  FIRR(3),FD(3),FDRAG(3),FDDRAG(3),R2DENS,
     &        RDISC,DDENS,DDENSDOT,VXRE,VYRE,VZRE,VRE,R2_DELTA
     &        CKFDRAG,DVXRE,DVYRE,DVZRE,DVRE,KFDOT,R_DOT 
*
      LOGICAL LPR1,E_BIND(NMAX)
      COMMON/ENERG/ E_BIND
*
*
*              Auxiliary Calculations
      L=NAME(I)      
*      IF(TIME.LT.1.0) THEN
*      IF(L.NE.2759) GOTO 10
      RDISC=(X(1,I)*X(1,I)+X(2,I)*X(2,I))**0.5
      R_DOT=(X(1,I)*XDOT(1,I)+X(2,I)*XDOT(2,I))/RDISC
      R2_DELTA=(RDISC**2.0+(D_DELTA/RDISC)**6.0)**-0.375*
     &         EXP(-(RDISC**4.0))            
*     
   20 CONTINUE
*
      DDENS=D_SIGM*R2_DELTA*EXP(-(X(3,I)**2.0)/(2.0*(HZ**2.0)))
      DDOTDENS=EXP(-(RDISC)**4.0)*
     & (RDISC**1.25*(RDISC**8.0+D_DELTA**6.0)**-0.375*
     & (-0.75*(RDISC**8.0-3.0*D_DELTA**6.0)/
     & (RDISC**8.0+D_DELTA**6.0)-
     & 4.0*RDISC**2.0))*R_DOT-DDENS*X(3,I)*XDOT(3,I)/HZ**2.0
*
      VXDISC=-(SMBH**0.5)*X(2,I)/((RDISC)**1.5)
      VYDISC=+(SMBH**0.5)*X(1,I)/((RDISC)**1.5)
      VXRE=VXDISC-XDOT(1,I)
      VYRE=VYDISC-XDOT(2,I)
      VZRE=-XDOT(3,I)
      VRE=DSQRT(VXRE*VXRE + VYRE*VYRE + VZRE*VZRE)
      CKFDRAG=Q_DRAG*DDENS*VRE
*
*              Mean of Drag Force
      FDRAG(1)=CKFDRAG*VXRE
      FDRAG(2)=CKFDRAG*VYRE
      FDRAG(3)=CKFDRAG*VZRE
*      	 
      IF(TIME.EQ.0) EDISS1_LOC=EDISS1_LOC+BODY(I)*(FDRAG(1)*XDOT(1,I)+
     &FDRAG(2)*XDOT(2,I)+FDRAG(3)*XDOT(3,I))*STEP(I)
      LPR1=NAME(I).EQ.2074
      IF(LPR1) THEN
      IF(TTOT.EQ.0) EDISS1_LOC1=EDISS1_LOC1+BODY(I)*(FDRAG(1)*XDOT(1,I)+
     &FDRAG(2)*XDOT(2,I)+FDRAG(3)*XDOT(3,I))*STEP(I)
      END IF 
*              Auxiliary Calculations
      DVXRE = -(SMBH**0.5)*XDOT(2,I)*(RDISC**(-1.5))
     &        +(SMBH**0.5)*1.5*X(2,I)*R_DOT*(RDISC**(-2.5))
     &        -(FI(1,I)+FDRAG(1))
      DVYRE = +(SMBH**0.5)*XDOT(1,I)*(RDISC**(-1.5))
     &        -(SMBH**0.5)*1.5*X(1,I)*R_DOT*(RDISC**(-2.5))
     &        -(FI(2,I)+FDRAG(2))
      DVZRE = -(FI(3,I)+FDRAG(3))
*
      DVRE=(VXRE*DVXRE+VYRE*DVYRE+VZRE*DVZRE)/VRE
      KFDOT=Q_DRAG*DDENS
*
*   
*              Mean of Drag Force Derivative
      FDDRAG(1)=KFDOT*(DDOTDENS*VXRE*VRE+DVRE*VXRE+DVXRE*VRE)
      FDDRAG(2)=KFDOT*(DDOTDENS*VYRE*VRE+DVRE*VYRE+DVYRE*VRE)
      FDDRAG(3)=KFDOT*(DDOTDENS*VZRE*VRE+DVRE*VZRE+DVZRE*VRE)
*   
*     END IF
*    
*              Total Force Acting on a Star
      DO K=1,3
         FI(K,I)=FI(K,I)+FDRAG(K)
         D1(K,I)=D1(K,I)+FDDRAG(K)
      END DO 
*
  10  CONTINUE
      RETURN
*
      END
      




