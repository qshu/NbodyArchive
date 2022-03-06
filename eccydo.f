      SUBROUTINE ECCYDO
*
*     Black Hole output & orbit parameter computation
*     called in nbint  line 122
*
      INCLUDE 'common6.h'
*
*     Defining Variables for following computations
      REAL*8  XFBH(1:3),XSBH(1:3),VFBH(1:3),VSBH(1:3)
      REAL*8  RELRBH(1:3),RELVBH(1:3),LREL(1:3)
      REAL*8  NOVEC(1:3),RLVEC(1:3),CMBH(1:3),CMVBH(1:3)
      REAL*8  NORMFH(1:3),NORMSH(1:3)
      REAL*8  DENSR(1:7),DENS(1:6),DIST(1:3)
      REAL*8  ZAVZ,TIMESA,OUTST
      INTEGER FNAME,SNAME,NUM,COUNT(1:6)
      REAL*8  MFBH,MSBH,MRED
      REAL*8  RBHN,RBHDOT,LRELN,HENERG
      REAL*8  ECC,INC,SEMIMA,LOMEGA,SOMEGA,SUMOME
      REAL*8  COSI,COSLOM,COSSOM
      REAL*8  RLVECN,NOVECN,PROJRL
      REAL*8  NNORFH,NNORSH,COSIFH,COSISH,INCFBH,INCSBH
      REAL*8  DISTN
*
      INTEGER NUMB,ZAEHL
      REAL*8  CMVBH2,NVCOM2,MVCOM2,SVCOM2,BRACK,BRACKT
      REAL*8  NECC,NECY,MECC,MECY,SECC,SECY
      REAL*8  BRECK,BRECKT,BRYCK,BRYCKT
      REAL*8  VCOM2(1:30000),ECCAR(1:30000),ECYAR(1:30000)
*
*     Defining Variables of compare calculation
      REAL*8  EI2,XLI2,CHECK,VREL2,RDOT,XMBIN,ECY,ETWO,RITWO
      REAL*8  XIARG,XINCL,XOARG,XOMEG,SEMI,EBIND
      INTEGER KP1,KP2
      REAL*8  EI(1:3),XLI(1:3)
*
      REAL*8  RELAB,SAVEAB,MINAB
      INTEGER M1NAME,M1SAVE,M2NAME,M2SAVE,NCHBH1,NCHBH2,PLACE
*
      REAL*8  COSAVE(1:6,1:2000),DESAVE(1:6,1:2000)
      REAL*8  COSUM(1:6),DESUM(1:6),MCOUNT(1:6),MDENS(1:6)
      REAL*8  BRCOU(1:6),BRCOUT(1:6),BRDEN(1:6),BRDENT(1:6)
      REAL*8  SCOUNT(1:6),SDENS(1:6),ERRTIM
      REAL*8  XNBBH1(1:3),XNBBH2(1:3)
      INTEGER NUMOV,ERROUT,SCHBH1,SCHBH2
*
*     Look for changes in neighbour list
      M1SAVE=M1NAME
      M2SAVE=M2NAME
      SAVEAB=50000000
      DO 10 J=1,NTOT
        IF(NAME(J).EQ.(NZERO-1))THEN
          DO 12 K=2,(LIST(1,J)+1)
            PLACE=LIST(K,J)
            RELAB=DSQRT((X(1,J)-X(1,PLACE))**2+
     &	              (X(2,J)-X(2,PLACE))**2+(X(3,J)-X(3,PLACE))**2)
	    IF(RELAB.LE.SAVEAB)THEN
	      M1NAME=NAME(PLACE)
	      SAVEAB=RELAB
	    ENDIF
   12     CONTINUE
        ENDIF
   10 CONTINUE
      SAVEAB=5000000
      DO 14 J=1,NTOT
        IF(NAME(J).EQ.(NZERO))THEN
          DO 16 K=2,(LIST(1,J)+1)
            PLACE=LIST(K,J)
            RELAB=DSQRT((X(1,J)-X(1,PLACE))**2+
     &	              (X(2,J)-X(2,PLACE))**2+(X(3,J)-X(3,PLACE))**2)
	    IF(RELAB.LE.SAVEAB)THEN
	      M2NAME=NAME(PLACE)
	      SAVEAB=RELAB
	    ENDIF
   16     CONTINUE
        ENDIF
   14 CONTINUE
      IF(M1NAME.NE.M1SAVE.AND.M1NAME.NE.NZERO)THEN
        NCHBH1=NCHBH1+1
        SCHBH1=SCHBH1+1
      ENDIF
      IF(M2NAME.NE.M2SAVE.AND.M2NAME.NE.(NZERO-1))THEN
        NCHBH2=NCHBH2+1
        SCHBH2=SCHBH2+1
      ENDIF
*
*     Output step is given by OUTST
      OUTST=0.002
*
*     Look whether new output is necessary
      ZAVZ=TIME-TIMESA
      IF(ZAVZ.LT.OUTST)THEN
        RETURN
      ENDIF
      TIMESA=TIME
*
*     Intervals for error output
      ERROUT=496
      IF(NUMOV.EQ.ERROUT)THEN
      NUMOV=0
      DO 20 J=1,6
         COSUM(J)=0.D0
	 DESUM(J)=0.D0
	 BRCOUT(J)=0.D0
	 BRDENT(J)=0.D0
	 MCOUNT(J)=0.D0
	 MDENS(J)=0.D0
	 BRCOU(J)=0.D0
	 BRDEN(J)=0.D0
	 SCOUNT(J)=0.D0
	 SDENS(J)=0.D0
   20 CONTINUE
      ENDIF
*
*     Save special coordinates (Black Holes on initial positions
*                                NZERO-1 AND NZERO)
      DO 47 J=1,NTOT
        IF(NAME(J).EQ.(NZERO-1))THEN
          XFBH(1)=X(1,J)
          XFBH(2)=X(2,J)
          XFBH(3)=X(3,J)
          VFBH(1)=XDOT(1,J)
          VFBH(2)=XDOT(2,J)
          VFBH(3)=XDOT(3,J)
          FNAME=NAME(J)
          MFBH=BODY(J)
        ELSE IF(NAME(J).EQ.M1SAVE)THEN
	  XNBBH1(1)=X(1,J)
	  XNBBH1(2)=X(2,J)
	  XNBBH1(3)=X(3,J)
        ENDIF
   47 CONTINUE
      DO 48 J=1,NTOT
        IF(NAME(J).EQ.NZERO)THEN
          XSBH(1)=X(1,J)
          XSBH(2)=X(2,J)
          XSBH(3)=X(3,J)
          VSBH(1)=XDOT(1,J)
          VSBH(2)=XDOT(2,J)
          VSBH(3)=XDOT(3,J)
          SNAME=NAME(J)
          MSBH=BODY(J)
        ELSE IF(NAME(J).EQ.M2SAVE)THEN
	  XNBBH2(1)=X(1,J)
	  XNBBH2(2)=X(2,J)
	  XNBBH2(3)=X(3,J)
        ENDIF
   48 CONTINUE
*
*     Compute reduced mass
      MRED=(MFBH*MSBH)/(MFBH+MSBH)
*
*     Compute relative space/velocity coordinates
      RELRBH(1)=XFBH(1)-XSBH(1)
      RELRBH(2)=XFBH(2)-XSBH(2)
      RELRBH(3)=XFBH(3)-XSBH(3)
      RELVBH(1)=VFBH(1)-VSBH(1)
      RELVBH(2)=VFBH(2)-VSBH(2)
      RELVBH(3)=VFBH(3)-VSBH(3)
*
*     Compute mean of relative space coordinates
      RBHN=DSQRT((XFBH(1)-XSBH(1))**2+(XFBH(2)-XSBH(2))**2
     &                              +(XFBH(3)-XSBH(3))**2)
*
*     Compute RBHN dot
      RBHDOT=(RELVBH(1)*RELRBH(1)+RELVBH(2)*RELRBH(2)
     &                           +RELVBH(3)*RELRBH(3))/RBHN
*
*     Compute relative angular momentum
      LREL(1)=MRED*(RELVBH(3)*RELRBH(2)-RELVBH(2)*RELRBH(3))
      LREL(2)=MRED*(RELVBH(1)*RELRBH(3)-RELVBH(3)*RELRBH(1))
      LREL(3)=MRED*(RELVBH(2)*RELRBH(1)-RELVBH(1)*RELRBH(2))
      LRELN=DSQRT(LREL(1)**2+LREL(2)**2+LREL(3)**2)
*
*     Compute energy
      HENERG=(0.5*MRED*RBHDOT**2)+(LRELN**2/(2*MRED*RBHN**2))
     &                           -(MFBH*MSBH/RBHN)
*
*     Compute relative eccentricity ECC
      ECC=DSQRT(1+(2*HENERG*LRELN**2)/(MRED*(MFBH*MSBH)**2))
*
*     Compute semi-major axis SEMIMA in relative System
      SEMIMA=-MFBH*MSBH/(2*HENERG)
*
*     Compute relative inclination INC
      COSI=LREL(3)/LRELN
      INC=ACOS(COSI)
*
*     Compute relative node-vector NOVEC
      NOVEC(1)=-LREL(2)
      NOVEC(2)=LREL(1)
      NOVEC(3)=0.0D0
      NOVECN=DSQRT(NOVEC(1)**2+NOVEC(2)**2)
*
*     Compute relative longitude of the ascending node LOMEGA
      IF(NOVECN.EQ.0)THEN
         LOMEGA=0.0D0
      ELSE
         COSLOM= -LREL(2)/NOVECN
         IF(NOVEC(2).GT.0) LOMEGA=ACOS(COSLOM)
         IF(NOVEC(2).LT.0) LOMEGA=TWOPI - ACOS(COSLOM)
      ENDIF
*
*     Compute relative Runge-Lenz Vector RLVEC
      RLVEC(1)=(LREL(3)*RELVBH(2)-LREL(2)*RELVBH(3))/(MFBH*MSBH)
     &                           -RELRBH(1)/RBHN
      RLVEC(2)=(LREL(1)*RELVBH(3)-LREL(3)*RELVBH(1))/(MFBH*MSBH)
     &                           -RELRBH(2)/RBHN
      RLVEC(3)=(LREL(2)*RELVBH(1)-LREL(1)*RELVBH(2))/(MFBH*MSBH)
     &                           -RELRBH(3)/RBHN
      RLVECN=DSQRT(RLVEC(1)**2+RLVEC(2)**2+RLVEC(3)**2)
      PROJRL=DSQRT(RLVEC(1)**2+RLVEC(2)**2)
*
*     Compute relative perihelion distance SOMEGA
*     Attention: Angel between projected RL Vector
*                and x-unityvector calculated
      IF(NOVECN.EQ.0)THEN
         COSSOM=RLVEC(1)/RLVECN
         SOMEGA=ACOS(COSSOM)
      ELSE
         COSSOM=(RLVEC(1)/PROJRL)
         IF(RLVEC(2).GT.0) SOMEGA=ACOS(COSSOM)
         IF(RLVEC(2).LT.0) SOMEGA=TWOPI - ACOS(COSSOM)
      ENDIF
*
*
*     Compute inclination in absolute coordinates
      NORMFH(1)=VFBH(3)*XFBH(2)-VFBH(2)*XFBH(3)
      NORMFH(2)=VFBH(1)*XFBH(3)-VFBH(3)*XFBH(1)
      NORMFH(3)=VFBH(2)*XFBH(1)-VFBH(1)*XFBH(2)
      NORMSH(1)=VSBH(3)*XSBH(2)-VSBH(2)*XSBH(3)
      NORMSH(2)=VSBH(1)*XSBH(3)-VSBH(3)*XSBH(1)
      NORMSH(3)=VSBH(2)*XSBH(1)-VSBH(1)*XSBH(2)
      NNORFH=DSQRT(NORMFH(1)**2+NORMFH(2)**2+NORMFH(3)**2)
      NNORSH=DSQRT(NORMSH(1)**2+NORMSH(2)**2+NORMSH(3)**2)
      COSIFH=NORMFH(3)/NNORFH
      COSISH=NORMSH(3)/NNORSH
      INCFBH=ACOS(COSIFH)
      INCSBH=ACOS(COSISH)
*
*     Compute C.M. motion of BH binary
      DO 70 K=1,3
      CMBH(K)=(MFBH*XFBH(K)+MSBH*XSBH(K))/(MFBH+MSBH)
   70 CONTINUE
      DO 71 K=1,3
      CMVBH(K)=(MFBH*VFBH(K)+MSBH*VSBH(K))/(MFBH+MSBH)
   71 CONTINUE
*
*     Alternative calculation for comparison
*     Variables
      VREL2=RELVBH(1)**2+RELVBH(2)**2+RELVBH(3)**2
      RDOT=RBHDOT*RBHN
      XMBIN=MFBH+MSBH
      RITWO=RELRBH(1)**2+RELRBH(2)**2+RELRBH(3)**2
*
*       Construct Runge-Lenz vector (Heggie & Rasio, 1995, IAU174, Eq.(5)).
          EI2 = 0.0D0
          XLI2 = 0.0D0
          CHECK = 0.0D0
          DO 76 K = 1,3
              KP1 = 1 + MOD(K,3)
              KP2 = 1 + MOD(K+1,3)
              EI(K) = (VREL2*RELRBH(K) - RDOT*RELVBH(K))/XMBIN -
     &                                        RELRBH(K)/DSQRT(RITWO)
              XLI(K)=RELRBH(KP1)*RELVBH(KP2)-RELRBH(KP2)*RELVBH(KP1)
              EI2 = EI2 + EI(K)**2
              XLI2 = XLI2 + XLI(K)**2
              CHECK = CHECK + EI(K)*XLI(K)
   76     CONTINUE
*
          XIARG = DSQRT((XLI(1)**2+XLI(2)**2)/XLI2)
          IF (XLI(3).GT.0.D0) THEN
              XINCL = ASIN(XIARG)
              ELSE
              XINCL = TWOPI/2.D0 - ASIN(XIARG)
          END IF
*
          XOARG = EI(2)/DSQRT(EI(1)**2+EI(2)**2)
*         XOARG = SQRT((EI(3)*XLI(1))**2+(EI(3)*XLI(2))**2+
*    &                 (EI(1)*XLI(1)+EI(2)*XLI(2))**2)/
*    &           (DSQRT(EI(1)**2+EI(2)**2)*SQRT(XLI(1)**2+XLI(2)**2))
          IF (EI(1).GT.0.D0) THEN
          IF (EI(2).GT.0.D0) XOMEG = ASIN(XOARG)
          IF (EI(2).LT.0.D0) XOMEG = TWOPI - ASIN(-XOARG)
          ELSE
          IF (EI(2).LT.0.D0) XOMEG = TWOPI/2.D0 + ASIN(-XOARG)
          IF (EI(2).GT.0.D0) XOMEG = TWOPI/2.D0 - ASIN(XOARG)
          END IF
*
*       Determine semi-major axis & eccentricity of inner binary.
      SEMI = 2.0D0/RBHN - VREL2/XMBIN
      SEMI = 1.0/SEMI
*       Save new binary data for gaseous model.
      EBIND = 0.5D0*MSBH*MFBH/SEMI
*
      ETWO = (1.0D0 - RBHN/SEMI)**2 + RDOT**2/(SEMI*XMBIN)
      IF(ETWO.GT.0.D0)THEN
      ECY = SQRT(ETWO)
      ELSE
      ECY = SQRT(-ETWO)
      END IF
*
*     Compute density around first BH
*     The following radii are used for calculations
      DENSR(7)=0.00
      DENSR(6)=0.01
      DENSR(5)=0.02
      DENSR(4)=0.04
      DENSR(3)=0.08
      DENSR(2)=0.16
      DENSR(1)=0.32
*     Set COUNT(NUM) to zero
      DO 78 NUM=1,6
      COUNT(NUM)=0
   78 CONTINUE
*
      DO 80 J=1,NTOT
        IF(NAME(J).NE.(NZERO-1))THEN
            DO 82 K=1,3
            DIST(K)=X(K,J)-XFBH(K)
   82       CONTINUE
            DISTN=DSQRT(DIST(1)**2+DIST(2)**2+DIST(3)**2)
            IF (DISTN.LE.DENSR(6)) THEN
              COUNT(6)=COUNT(6)+1
            ELSE IF (DISTN.LE.DENSR(5)) THEN
              COUNT(5)=COUNT(5)+1
            ELSE IF (DISTN.LE.DENSR(4)) THEN
              COUNT(4)=COUNT(4)+1
            ELSE IF (DISTN.LE.DENSR(3)) THEN
              COUNT(3)=COUNT(3)+1
            ELSE IF (DISTN.LE.DENSR(2)) THEN
              COUNT(2)=COUNT(2)+1
            ELSE IF (DISTN.LE.DENSR(1)) THEN
              COUNT(1)=COUNT(1)+1
            ENDIF
        ENDIF
   80 CONTINUE
*     Compute densities
      DO 84 NUM=1,6
        DENS(NUM)=COUNT(NUM)/
     &          (2*TWOPI*(DENSR(NUM)**3-DENSR(NUM+1)**3)/3)
   84 CONTINUE
*     Compute mean densities/particle numbers and errors
*     Attention: The number-of-values-variable NUMOV is set to
*                zero as like as the mean and error values in
*    	         lines 93-108.
      NUMOV=NUMOV+1
      DO 85 NUM=1,6
        COSAVE(NUM,NUMOV)=COUNT(NUM)
	DESAVE(NUM,NUMOV)=DENS(NUM)
	COSUM(NUM)=COSUM(NUM)+COSAVE(NUM,NUMOV)
	DESUM(NUM)=DESUM(NUM)+DESAVE(NUM,NUMOV)
	MCOUNT(NUM)=COSUM(NUM)/NUMOV
	MDENS(NUM)=DESUM(NUM)/NUMOV
	DO 86 ZAEHL=1,NUMOV
	  BRCOU(NUM)=(MCOUNT(NUM)-COSAVE(NUM,ZAEHL))**2
	  BRCOUT(NUM)=BRCOUT(NUM)+BRCOU(NUM)
	  BRDEN(NUM)=(MDENS(NUM)-DESAVE(NUM,ZAEHL))**2
	  BRDENT(NUM)=BRDENT(NUM)+BRDEN(NUM)
   86   CONTINUE
        SCOUNT(NUM)=DSQRT(BRCOUT(NUM)/(NUMOV*(NUMOV-1)))
	SDENS(NUM)=DSQRT(BRDENT(NUM)/(NUMOV*(NUMOV-1)))
   85 CONTINUE
      DO 87 NUM=1,6
        BRCOUT(NUM)=0.D0
	BRDENT(NUM)=0.D0
   87 CONTINUE	   	  
*
*     Time for error output
      ERRTIM=TIME-0.5
*      
*     Compute <VCOM^2>,<ECC>,<ECY>
      IF(TIME.GT.22)THEN
        CMVBH2=CMVBH(1)**2+CMVBH(2)**2+CMVBH(3)**2
        NUMB=NUMB+1
        VCOM2(NUMB)=CMVBH2
        ECCAR(NUMB)=ECC
        ECYAR(NUMB)=ECY
        NVCOM2=NVCOM2+CMVBH2
        MVCOM2=NVCOM2/NUMB
        NECC=NECC+ECC
        NECY=NECY+ECY
        MECC=NECC/NUMB
        MECY=NECY/NUMB
        DO 88 ZAEHL=1,NUMB
          BRACK=(VCOM2(ZAEHL)-MVCOM2)**2
          BRACKT=BRACKT+BRACK
          BRECK=(ECCAR(ZAEHL)-MECC)**2
          BRECKT=BRECKT+BRECK
          BRYCK=(ECYAR(ZAEHL)-MECY)**2
          BRYCKT=BRYCKT+BRYCK
   88   CONTINUE
        SVCOM2=DSQRT(BRACKT/(NUMB*(NUMB-1)))
        SECC=DSQRT(BRECKT/(NUMB*(NUMB-1)))
        SECY=DSQRT(BRYCKT/(NUMB*(NUMB-1)))
	BRACKT=0.D0
	BRECKT=0.D0
	BRYCKT=0.D0
      ENDIF
*
*     Data output
      WRITE(*,90) TIME,FNAME,MFBH,XFBH(1),XFBH(2),XFBH(3),
     &                            VFBH(1),VFBH(2),VFBH(3),
     &                 SNAME,MSBH,XSBH(1),XSBH(2),XSBH(3),
     &                            VSBH(1),VSBH(2),VSBH(3),
     &                      RELRBH(1),RELRBH(2),RELRBH(3),
     &                      RELVBH(1),RELVBH(2),RELVBH(3),
     &                            CMBH(1),CMBH(2),CMBH(3),
     &                         CMVBH(1),CMVBH(2),CMVBH(3),
     &                         RLVEC(1),RLVEC(2),RLVEC(3),
     &                 RBHN,HENERG,SEMIMA,LOMEGA,SOMEGA,
     &                 RLVECN,ECC,INCFBH,INCSBH,INC,TIMESA,
     &                        MVCOM2,SVCOM2,MECC,SECC,NUMB,
     &                        COUNT(6),MCOUNT(6),SCOUNT(6),
     &                        COUNT(5),MCOUNT(5),SCOUNT(5),
     &                        COUNT(4),MCOUNT(4),SCOUNT(4),
     &                        COUNT(3),MCOUNT(3),SCOUNT(3),
     &                        COUNT(2),MCOUNT(2),SCOUNT(2),
     &                        COUNT(1),MCOUNT(1),SCOUNT(1),
     &                        DENS(6),MDENS(6),SDENS(6),
     &                        DENS(5),MDENS(5),SDENS(5),
     &                        DENS(4),MDENS(4),SDENS(4),
     &                        DENS(3),MDENS(3),SDENS(3),
     &                        DENS(2),MDENS(2),SDENS(2),
     &                        DENS(1),MDENS(1),SDENS(1),
     &                        NCHBH1,NCHBH2,SCHBH1,SCHBH2,
     &                        NUMOV,ERRTIM,
     &                     M1SAVE,XNBBH1(1),XNBBH1(2),XNBBH1(3),
     &                     M2SAVE,XNBBH2(1),XNBBH2(2),XNBBH2(3),
     &                     SEMI,EBIND,ECY,MECY,SECY,XINCL,XOMEG
   90 FORMAT('     BHOrbitdata',F15.6,I9,E15.4,E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                I9,E15.4,E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E20.9,E20.9,E20.9,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                             E15.4,E15.4,E15.4,E15.4,E15.4,
     &                       E15.4,E15.4,E15.4,E15.4,E15.4,E15.4,
     &                                E15.4,E15.4,E15.4,E15.4,I8,
     &                                            I6,E15.4,E15.4,
     &                                            I6,E15.4,E15.4,
     &                                            I6,E15.4,E15.4,
     &                                            I6,E15.4,E15.4,
     &                                            I6,E15.4,E15.4,
     &                                            I6,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                           I12,I12,I12,I12,
     &                                                 I10,E15.4,
     &                                      I9,E15.4,E15.4,E15.4,
     &                                      I9,E15.4,E15.4,E15.4,
     &                 E15.4,E15.4,E15.4,E15.4,E15.4,E15.4,E15.4)
*
*     Set intervals for neighbour change counters
      NCHBH1=0
      NCHBH2=0
*
      IF(NUMOV.EQ.ERROUT)THEN
         SCHBH1=0
         SCHBH2=0
      ENDIF
*
      RETURN
*
      END













