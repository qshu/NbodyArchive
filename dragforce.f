       SUBROUTINE DRAGFORCE(I,XI,XIDOT,FIRR,FD,FREG,ICASE)
*
*
*       Drag force & first derivative.
*       -------------------------
*
*              Types of BODY(I),X(K,I),XDOT(K,I) is Defined in "comon6.h" 
      INCLUDE 'common6.h'
      REAL*8  FIRR(3),FD(3),FDRAG(3),FDDRAG(3),FTOT(3),FREG(3),
     &        RDISC,DDENS,DDOTDENS,VXRE,VYRE,VZRE,VRE,
     &        CKFDRAG,DVXRE,DVYRE,DVZRE,DVRE,KFDOT,XI(3),XIDOT(3),
     &        R_DOT, R_DOT1,R2_DELTA,R_MINI, V_TI,V_XI
      LOGICAL LPR1, E_BIND(NMAX)
      COMMON/ENERG/ E_BIND
*
*        Calling Parameter: ICASE=0 from nbint.f
*                                =1 from regint.f
*                                =2 from kspert.f, fpert.f
*
      L=NAME(I)
*     IF(TIME.LT.1.0) THEN
*     IF(L.NE.2759) GOTO 10
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
      DO K=1,3
      XI(K)=X(K,I)
      XIDOT(K)=XDOT(K,I)
      END DO
*               Auxiliary Calculations
      R_MINI=SQRT(XI(1)*XI(1)+XI(2)*XI(2)+XI(3)*XI(3))
      V_XI=SQRT(XIDOT(1)*XIDOT(1)+XIDOT(2)*XIDOT(2)+XIDOT(3)*XIDOT(3))/
     &     SQRT(SMBH/R_MINI)
*      IF(R_MINI.LT.9.0E-04.AND.V_XI.LE.1.4) THEN
      IF(STEP(I).LT.9.0E-07.AND.NAME(I).NE.NAME(ICOMP).AND.
     &        V_XI.LE.1.5) THEN
*      IF(rank.eq.0)
*     &PRINT*,'RMINI=, name, stepi, v',R_MINI,NAME(I),STEP(I),V_XI
*      PRINT*,'STEP I',NAME(I),STEP(I),R_MINI,NAME(ICOMP),NAME(JCOMP)
*      PRINT*,'STEP I==============================================='
*      NTOT=NTOT-1
*      N=N-1
*      CALL REMOVE(I,1)
*      PRINT*,'remove is done',NTOT,N
      END IF
*
      RDISC=(XI(1)*XI(1)+XI(2)*XI(2))**0.5
      R_DOT=(XI(1)*XIDOT(1)+XI(2)*XIDOT(2))/RDISC
      R2_DELTA=(RDISC**2.0+(D_DELTA/RDISC)**6.0)**-0.375*
     &         EXP(-(RDISC**4.0)) 
*
   20 CONTINUE
* 
      DDENS=D_SIGM*R2_DELTA*EXP(-(XI(3)**2.0)/(2.0*(HZ**2.0)))
      DDOTDENS=EXP(-(RDISC)**4.0)*
     & (RDISC**1.25*(RDISC**8.0+D_DELTA**6.0)**-0.375*
     & (-0.75*(RDISC**8.0-3.0*D_DELTA**6.0)/(RDISC**8.0+D_DELTA**6.0)-
     &  4.0*RDISC**2.0))*R_DOT-DDENS*XI(3)*XI(3)/HZ**2.0
*
      VXDISC=-(SMBH**0.5)*XI(2)/((RDISC)**1.5)
      VYDISC=+(SMBH**0.5)*XI(1)/((RDISC)**1.5)
      VXRE=VXDISC-XIDOT(1)
      VYRE=VYDISC-XIDOT(2)
      VZRE=-XIDOT(3)
      VRE=(VXRE*VXRE + VYRE*VYRE + VZRE*VZRE)**0.5
      CKFDRAG=Q_DRAG*DDENS*VRE
*
*              Mean of Drag Force
      FDRAG(1)=CKFDRAG*VXRE
      FDRAG(2)=CKFDRAG*VYRE
      FDRAG(3)=CKFDRAG*VZRE  
*      IF(NAME(I).EQ.3162)
*     &PRINT*,'FDRAG=',(FDRAG(1)**2.0+FDRAG(2)**2.0+FDRAG(3)**2.0)**0.5
* 
      IF(ICASE.EQ.0) EDISS1_LOC=EDISS1_LOC+BODY(I)*(FDRAG(1)*XIDOT(1)+
     &FDRAG(2)*XIDOT(2)+FDRAG(3)*XIDOT(3))*STEP(I)
      LPR1=NAME(I).EQ.2074
      IF(LPR1) THEN
      IF(ICASE.EQ.0) EDISS1_LOC1=EDISS1_LOC1+BODY(I)*(FDRAG(1)*XIDOT(1)+
     &FDRAG(2)*XIDOT(2)+FDRAG(3)*XIDOT(3))*STEP(I)
      END IF

*      PRINT*,'DRAG FORCE: EDISS1_LOC & rank',rank,EDISS1_LOC
*
*              Auxiliary Calculations
      DVXRE = -(SMBH**0.5)*XIDOT(2)*(RDISC**(-1.5))
     &        +(SMBH**0.5)*1.5*XI(2)*R_DOT*(RDISC**(-2.5))
     &        -(FTOT(1)+FDRAG(1))
      DVYRE = +(SMBH**0.5)*XIDOT(1)*(RDISC**(-1.5))
     &        -(SMBH**0.5)*1.5*XI(1)*R_DOT*(RDISC**(-2.5))
     &        -(FTOT(2)+FDRAG(2))
      DVZRE = -(FTOT(3)+FDRAG(3))
*
      DVRE=(VXRE*DVXRE+VYRE*DVYRE+VZRE*DVZRE)/VRE
      KFDOT=Q_DRAG*DDENS
*
*   
*              Mean of Drag Force Derivative
      FDDRAG(1)=KFDOT*(DDOTDENS*VXRE*VRE+DVRE*VXRE+DVXRE*VRE)
      FDDRAG(2)=KFDOT*(DDOTDENS*VYRE*VRE+DVRE*VYRE+DVYRE*VRE)
      FDDRAG(3)=KFDOT*(DDOTDENS*VZRE*VRE+DVRE*VZRE+DVZRE*VRE)
*     PRINT*,FDDRAG(1),FDDRAG(2),FDDRAG(3)
*   
*     END IF

*  
      fdx = dsqrt(fdrag(1)**2+fdrag(2)**2+fdrag(3)**2)
      fddx = dsqrt(fddrag(1)**2+fddrag(2)**2+fddrag(3)**2)
      stepdd = fdx/fddx
*      if(stepdd.lt.1.0E-08)print*,
*     &         ' i,dt,dtd=',i,name(i),step(i),stepdd
*      if(stepdd.lt.step(i).and.NAME(I).EQ.324)print*,
*     &         ' i,fdx,fddx,R,Z=',i,fdx,fddx,RDISC,XI(3)
*      if(NAME(I).EQ.2074)print*,
*     &  'dr=',
*      if(name(i).eq.3547)print*,' !!!!!!!i,dt,dtd=',fdx,step(i),stepdd
*              Total Force Acting on a Star
      DO K=1,3
         FIRR(K)=FIRR(K)+FDRAG(K)
         FD(K)=FD(K)+FDDRAG(K)
      END DO 
*
  10  CONTINUE
      DO K=1,3
      XI(K)=X(K,I)
      XIDOT(K)=XDOT(K,I)
      END DO
*               Auxiliary Calculations
      RDISC=(XI(1)*XI(1)+XI(2)*XI(2)+XI(3)*XI(3))**0.5
      V_TI=step(i)*SQRT(XIDOT(1)**2.0+XIDOT(2)**2.0+XIDOT(3)**2.0)
*      if(rank.eq.0.and.RDISC.LE.0.2.and.V_TI.GE.5.0E-3.AND.
*     &DABS(XI(3)).LE.HZ) 
*     &PRINT*,',i, r, shag_i=',I, NAME(I),NAME(JCOMP),
*     &                NAME(ICOMP), JCOMP, ICOMP, RDISC, V_TI
      RETURN
*
      END


