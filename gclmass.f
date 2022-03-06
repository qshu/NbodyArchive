      SUBROUTINE GCLMASS(C)
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
     *      XMBH,CC,CFAC,XCOUL,XMCL,TCEN0,STEPCC,EG0,TGCL,APW1,APW3,
     *      DEDFF,ECLS,ECLS0,ETAG,RGAL,VGAL,XMGAL,XTIDAL,RTIDAL,
     *      IMEM(NMAX),ICLUST,IDUMMY
*
      REAL*8  R2(NMAX),C(3)
*
      PI43 = TWOPI*2.D0/3.D0
      ALP = APW1 + 1.D0
      RGAL = DSQRT(XG(1)**2 + XG(2)**2 + XG(3)**2)
      XMGAL = CC*RGAL**ALP
      RHOGAL = (XMGAL+XMBH)/PI43/RGAL**3
      XMCLTT = XMCL
*
*       Set square radii of all single particles & c.m. bodies.
      NP = 0
      DO 10 I = 1,N
          NP = NP + 1
          R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                                  (X(3,I) - C(3))**2
          JLIST(NP) = I
   10 CONTINUE
*       Sort square distances with respect to the centre C.
      CALL SORT1(NP,R2,JLIST)
*
      ZM = 0.0D0
*
         DO 20 I = 1,N
         IM = JLIST(I)
         ZM = ZM + BODY(IM)
         RHOCL = ZM/PI43/R2(I)**3
*
*        Check whether mass within Langrangian radius is complete.
         IF (RHOCL.GT.RHOGAL) THEN
             RTIDAL = DSQRT(R2(I))
         ELSE
             IF (R2(I).LT.XTIDAL*RTIDAL)THEN
                XMCLTT = ZM
                GO TO 25
             END IF
         END IF
*     
 20   CONTINUE
*
 25   CONTINUE
*        Only apply if ICLUST GE 1
      IF(ICLUST.GE.1) XMCL = XMCLTT
*
      RETURN
*
      END


