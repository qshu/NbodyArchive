      SUBROUTINE OFFSET(DTOFF)
*
*
*       Offset of global times.
*       -----------------------
*
      INCLUDE 'common6.h'
      COMMON /GC/ XG0(3),VG0(3),AG0(3),ADOTG0(3),
     *      XG(3),VG(3),AG(3),ADOTG(3),AG2(3),AG3(3),
     *      XGI(3),VGI(3),AGI(3),ADOTGI(3),FTOT(3),FDTOT(3),
     *      ADF0(3),ADOTF0(3),ADF(3),ADOTF(3),
     *      XMBH,CC,CFAC,XCOUL,XMCL,TCEN0,STEPCC,EG0,TGCL,APW1,APW3,
     *      DEDFF,ECLS,ECLS0,ETAG,RGAL,VGAL,XMGAL,XTIDAL,RTIDAL,
     *      IMEM(NMAX),ICLUST,IDUMMY
*
*
*       Update the global offset time.
    1 TOFF = TOFF + DTOFF
*
*       Reduce all individual times and epochs by offset interval.
      DO 10 I = 1,NTOT
          T0(I) = T0(I) - DTOFF
          T0R(I) = T0R(I) - DTOFF
          TEV(I) = TEV(I) - DTOFF
          TEV0(I) = TEV0(I) - DTOFF
          EPOCH(I) = EPOCH(I) - DTOFF*TSTAR
   10 CONTINUE
*       Reduce global timers for galactic centre
      TCEN0 = TCEN0 - DTOFF
      TGCL = TGCL - DTOFF
*
*       Set new global times.
      TIME = TIME - DTOFF
      TADJ = TADJ - DTOFF
      TNEXT = TNEXT - DTOFF
      TPREV = TPREV - DTOFF
      IF (KZ(19).GT.2) THEN
          TPLOT = TPLOT - DTOFF
          TMDOT = TMDOT - DTOFF
      END IF
*
*       See whether more reductions are needed.
      IF (TIME.GE.TOFF) GO TO 1
*
*       Activate control indicator for new scheduling.
      IPHASE = -1
*
      RETURN
*
      END
