      SUBROUTINE INFLSMBH
*
*
*       Parameters of particles inside of influence radius ob SMBH.
*       --------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8 R_XYZ(NMAX),R_MINIM,STEP_M,VEL_S,RSIII
      REAL*8 R_MINIM_1,STEP_MR,VEL_S1,VEL_IN(NMAX),R_DIS
      INTEGER IMYA,IMYA_1,NUMBER1,
     &        NUMBER2,NUMBER3,NUMBER4,NUMBER5,NUMBER6,
     &        NUMBER7,NUMBER8
*
*
      R_MINIM=10.0
      RMAKC=1.0
      STEP_MR=1.0
      NUMBER1=0
      NUMBER2=0
      NUMBER3=0
      NUMBER4=0
      NUMBER5=0
      NUMBER6=0
*
      R_DIS=0.0
*
      DO I = IFIRST,NTOT
      R_XYZ(I)=(X(1,I)**2.0+X(2,I)**2.0+X(3,I)**2.0)**0.5 
*      IF(TIME.EQ.1.0) PRINT*,'list',LIST(1,I),LIST(1,N+I),LIST(1,3010)
*      IF(TIME.EQ.1.0) PRINT*,'radi',IFIRST,NTOT,ICOMP,JCOMP,NPAIRS
*
*
      IF(R_XYZ(I).LT.R_MINIM) THEN
      R_MINIM=R_XYZ(I)
      STEP_M=STEPR(I)
      IMYA=NAME(I)
      VEL_S=SQRT(XDOT(1,I)**2.0+XDOT(2,I)**2.0+XDOT(3,I)**2.0)/
     &      SQRT(SMBH/R_MINIM)
      NUMBER7=LIST(1,I)
*     RSIII=RS(I)
      END IF
      IF(STEPR(I).LT.STEP_MR) THEN
      STEP_MR=STEPR(I)
      R_MINIM_1=R_XYZ(I)
      IMYA_1=NAME(I)
      VEL_S1=SQRT(XDOT(1,I)**2.0+XDOT(2,I)**2.0+XDOT(3,I)**2.0)/
     &      SQRT(SMBH/R_MINIM)
      NUMBER8=LIST(1,I)
      END IF
      END DO
*
      NNB=0
      DO I = 1,N
      R_XYZ(I)=(X(1,I)**2.0+X(2,I)**2.0+X(3,I)**2.0)**0.5
      IF(R_XYZ(I).LT.D_SCALE/3.0) THEN
      NNB=NNB+1
      LIST(1,NTOT+1)=NNB
      LIST(I+1,NTOT+1)=I
      END IF
      END DO    
      PRINT*,'NNB==============================',NNB,D_SCALE,D_DELTA
*
* 
*
        DO I = 1,N
        J=NAME(I)   
        VEL_IN(I)=(XDOT(1,I)**2.0+XDOT(2,I)**2.0+XDOT(3,I)**2.0)**0.5/
     &           SQRT(SMBH/R_XYZ(I))
        IF(VEL_IN(I).LE.1.4.AND.VEL_IN(I).GT.1.0) NUMBER2=NUMBER2+1
        IF(VEL_IN(I).LE.1.0) NUMBER3=NUMBER3+1
*
        IF(R_XYZ(I).LE.D_DELTA) THEN
        NUMBER4=NUMBER4+1
        IF(VEL_IN(I).LE.1.4.AND.VEL_IN(J).GT.1.0) NUMBER5=NUMBER5+1
*      
        IF(VEL_IN(I).LE.1.0) THEN
        NUMBER6=NUMBER6+1
        END IF
*
        END IF
*
*         To remember names of particles since this list is updated
*         only each adjust time()
*
        NUMBER1=NUMBER1+1
        IF(I.EQ.INT(N*SMBH)) THEN
*       LIST(1,NTOT+1)=I
*       D_SCALE=R_XYZ(I)
        GO TO 75
        END IF
        END DO 
   75 CONTINUE
*
      PRINT*,'Num.of Stars < Circ & Parab.vel(R-infl):',NUMBER3,NUMBER2
      PRINT*,'Num.of Stars < Circ & Parab.vel(R-innc):',NUMBER6,NUMBER5
      PRINT*,'Num.of Stars inside R-innc & R-infl    :',NUMBER4,NUMBER1
      PRINT*,'R-incut & R-infl:',D_DELTA,D_SCALE
      PRINT*,'R-min &    STEPR:',R_MINIM,STEP_M
      PRINT*,'Name of R-min to SMBH & its Vel        :',IMYA, VEL_S
      PRINT*,'R for min STEPR :',R_MINIM_1,STEP_MR
      PRINT*,'Name of STEPR-min & its Vel.           :',IMYA_1,VEL_S1
      PRINT*,'List1,List2    :',NUMBER7,NUMBER8,SMBH
*
      RETURN
*
      END
