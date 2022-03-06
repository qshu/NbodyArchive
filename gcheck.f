      SUBROUTINE GCHECK
*
*      Galactic Centre and Bulge Forces
*      Dynamical Friction
*      Predictor/Corrector for Apparent Motion of Galactic Centre
*      R.Sp./ O.G. Nov. 2001
*
      INCLUDE 'common6.h'
      COMMON /GC/ XG0(3),VG0(3),AG0(3),ADOTG0(3),
     *      XG(3),VG(3),AG(3),ADOTG(3),AG2(3),AG3(3),
     *      XGI(3),VGI(3),AGI(3),ADOTGI(3),FTOT(3),FDTOT(3),
     *      ADF0(3),ADOTF0(3),ADF(3),ADOTF(3),
     *      XMBH,CC,CFAC,XCOUL,XMCL,TCEN0,STEPCC,EG0,
     *      DEDFF,ECLS,ECLS0,ETAG,RGAL,VGAL,XMGAL,IMEM(NMAX)
      REAL*8  W0(4),W1(4),W2(4),W3(4)
      COMMON/ETT/DEDF1,DEDF2,DEDF3,DEDF4
*
*     Predictor for galactic centre
*     Note that XG,VG contain gravitation plus dynamical friction
*
      S = TIME - TCEN0
      S1 = 1.5D0*S
      S2 = 2.0D0*S
      DO 1 K = 1,3
      FTOT(K) = (AG0(K) + ADF0(K))/2.D0
      FDTOT(K) = (ADOTG0(K) + ADOTF0(K))/6.D0
      XG(K) = ((FDTOT(K)*S + FTOT(K))*S +VG0(K))*S + XG0(K)
      VG(K) = (FDTOT(K)*S1 + FTOT(K))*S2 + VG0(K)
      AG(K) = ADOTG0(K)*S + AG0(K)
      ADOTG(K) = ADOTG0(K)
      ADF(K) = ADOTF0(K)*S + ADF0(K)
      ADOTF(K) = ADOTF0(K)
 1    CONTINUE
*
      PRINT*,rank,' GC: XG=',(XG(K),K=1,3),' ||=',
     *   DSQRT(XG(1)**2+XG(2)**2+XG(3)**2),' S=',S
      PRINT*,rank,'     VG=',(VG(K),K=1,3),' ||=',
     *   DSQRT(VG(1)**2+VG(2)**2+VG(3)**2)
      PRINT*,rank,'    CMR=',(CMR(K),K=1,3),' ||=',CMR(4)
      PRINT*,rank,' CMRDOT=',(CMRDOT(K),K=1,3),' ||=',CMRDOT(4)
      PRINT*,rank,'   FTOT=',(FTOT(K),K=1,3)
      PRINT*,rank,'  FDTOT=',(FDTOT(K),K=1,3)
      PRINT*,rank,' AG=',(AG(K),K=1,3),' ||=',
     *                  DSQRT(AG(1)**2+AG(2)**2+AG(3)**2)
      PRINT*,rank,' ADOTG=',(ADOTG(K),K=1,3),' ||=',
     *                  DSQRT(ADOTG(1)**2+ADOTG(2)**2+ADOTG(3)**2)
      PRINT*,rank,' ADF=',(ADF(K),K=1,3),' ||=',
     *                  DSQRT(ADF(1)**2+ADF(2)**2+ADF(3)**2)
      PRINT*,rank,' ADOTF=',(ADOTF(K),K=1,3),' ||=',
     *                  DSQRT(ADOTF(1)**2+ADOTF(2)**2+ADOTF(3)**2)
*
      ET0 = ET
      ET = 0.0D0
*
      DEDF1 = 0.0D0
      DEDF2 = 0.0D0
      DEDF3 = 0.0D0
      DEDF4 = 0.0D0
*
      DO 55 I = IFIRST,NTOT
         DO 56 K = 1,3
*      Start Energy Checks
 56      XGI(K) = XG(K) + X(K,I)
         RIJ2 = XGI(1)**2 + XGI(2)**2 + XGI(3)**2
         DR2I = 1.0/RIJ2
         DR3I = XMBH*DR2I*DSQRT(DR2I)
          PHIEXT1 = -DR3I*RIJ2
          PHIEXT2 = CC*DLOG(DSQRT(RIJ2)/RGAL)
*
          DEDF3 = DEDF3 + BODY(I)*PHIEXT1
          DEDF4 = DEDF4 + BODY(I)*PHIEXT2
          DEDF1 = DEDF1 + BODY(I)*(VG(1)**2+VG(2)**2+VG(3)**2)/2.D0
          DEDF2 = DEDF2 + BODY(I)*(VG(1)*XDOT(1,I) + VG(2)*XDOT(2,I) +
     *            VG(3)*XDOT(3,I))
 55    CONTINUE
*
       RGCM=DSQRT((XG(1)+CMR(1))**2+(XG(2)+CMR(2))**2+(XG(3)+CMR(3))**2)
       ECM = XMCL*CMRDOT(4)**2/2.D0
       ECLK = XMCL*(VG(1)**2+VG(2)**2+VG(3)**2)/2.D0
       RIJ2 = XG(1)**2 + XG(2)**2 + XG(3)**2
         DR2I = 1.0/RIJ2
         DR3I = XMBH*DR2I*DSQRT(DR2I)
         PHICLS = -DR3I*RIJ2 + CC*DLOG(DSQRT(RIJ2)/RGAL)
*
       ECLS = ECLK + XMCL*PHICLS
       ETIDE = DEDF1 + DEDF2 + DEDF3 + DEDF4 - EG0 - ECLS + ECLS0
       ETOT = EBIN + ZKIN - POT + ETIDE
       WRITE(6,1111)rank,EBIN,ZKIN,-POT,DEDF1,DEDF2,DEDF3,DEDF4,EG0
 1111  FORMAT(I5,' ETOT,i=',1P,9D15.7)
       PRINT*,rank,' ETOTx=',ETOT,EBIN+ZKIN-POT,ETIDE
       PRINT*,rank,' ETOT+.25=',ETOT+0.25D0,' ETIDE=',ETIDE,' ECM=',ECM
       PRINT*,rank,' ECLK,POTCL=',ECLK,XMCL*PHICLS,' ECLS,0=',
     *   ECLS,ECLS0,' diff=',ECLS-ECLS0,' DEDFF=',DEDFF
       CALL FLUSH(6)
*
 100  CONTINUE
*
      RETURN
*
      END
