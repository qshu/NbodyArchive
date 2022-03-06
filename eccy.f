      SUBROUTINE ECCY
*
*     Black Hole Output & Eccentricity Calculation
*     called by lagr.F
*
      INCLUDE 'common6.h'
*
*     Defining Variables for following computations
      REAL*8  XFBH(1:3),XSBH(1:3),VFBH(1:3),VSBH(1:3)
      REAL*8  RELRBH(1:3),RELVBH(1:3),LREL(1:3)
      REAL*8  NOVEC(1:3),RLVEC(1:3),CMBH(1:3),CMVBH(1:3)
      REAL*8  NORMFH(1:3),NORMSH(1:3)
      REAL*8  DENSR(1:7),DENS(1:6),DIST(1:3)
      INTEGER FNAME,SNAME,NUM,COUNT(1:6)
      REAL*8  MFBH,MSBH,MRED
      REAL*8  RBHN,RBHDOT,LRELN,HENERG
      REAL*8  ECC,INC,SEMIMA,LOMEGA,SOMEGA
      REAL*8  COSI,COSLOM,COSSOM
      REAL*8  RLVECN,NOVECN,PROJRL
      REAL*8  NNORFH,NNORSH,COSIFH,COSISH,INCFBH,INCSBH
      REAL*8  DISTN,TIMESA
*
*     Defining Variables of compare calculation
      REAL*8  EI2,XLI2,CHECK,VREL2,RDOT,XMBIN,ECY,ETWO,RITWO
      REAL*8  XIARG,XINCL,XOARG,XOMEG,SEMI,EBIND
      INTEGER KP1,KP2
      REAL*8  EI(1:3),XLI(1:3)
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
        ELSE
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
        ELSE
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
         IF(NOVEC(2).LT.0) LOMEGA=TWOPI - ACOS(COSSOM)
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
      DENSR(6)=0.05
      DENSR(5)=0.10
      DENSR(4)=0.15
      DENSR(3)=0.20
      DENSR(2)=0.25
      DENSR(1)=0.30
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
*
*     Data output
      PRINT*,''
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
     &      COUNT(6),DENS(6),COUNT(5),DENS(5),COUNT(4),DENS(4),
     &      COUNT(3),DENS(3),COUNT(2),DENS(2),COUNT(1),DENS(1)
   90 FORMAT('     BHAnalysis ',F15.6,I6,E15.4,E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                I6,E15.4,E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E20.9,E20.9,E20.9,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                             E15.4,E15.4,E15.4,E15.4,E15.4,
     &                       E15.4,E15.4,E15.4,E15.4,E15.4,E15.4,
     &                                I6,E15.4,I6,E15.4,I6,E15.4,
     &                                I6,E15.4,I6,E15.4,I6,E15.4)
*
*     Data output (Comparison)
      WRITE(*,95) TIME,FNAME,MFBH,XFBH(1),XFBH(2),XFBH(3),
     &                            VFBH(1),VFBH(2),VFBH(3),
     &                 SNAME,MSBH,XSBH(1),XSBH(2),XSBH(3),
     &                            VSBH(1),VSBH(2),VSBH(3),
     &                      RELRBH(1),RELRBH(2),RELRBH(3),
     &                      RELVBH(1),RELVBH(2),RELVBH(3),
     &                            CMBH(1),CMBH(2),CMBH(3),
     &                         CMVBH(1),CMVBH(2),CMVBH(3),
     &                                  EI(1),EI(2),EI(3),
     &                 RBHN,EBIND,SEMI,LOMEGA,XOMEG,
     &                 RLVECN,ECY,INCFBH,INCSBH,XINCL,TIMESA,
     &      COUNT(6),DENS(6),COUNT(5),DENS(5),COUNT(4),DENS(4),
     &      COUNT(3),DENS(3),COUNT(2),DENS(2),COUNT(1),DENS(1)
   95 FORMAT('     BHComAnalys',F15.6,I6,E15.4,E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                I6,E15.4,E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                                         E20.9,E20.9,E20.9,
     &                                         E15.4,E15.4,E15.4,
     &                                         E15.4,E15.4,E15.4,
     &                             E15.4,E15.4,E15.4,E15.4,E15.4,
     &                       E15.4,E15.4,E15.4,E15.4,E15.4,E15.4,
     &                                I6,E15.4,I6,E15.4,I6,E15.4,
     &                                I6,E15.4,I6,E15.4,I6,E15.4)
*
      RETURN
*
      END



