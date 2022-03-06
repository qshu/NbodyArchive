      SUBROUTINE GALCEN(I,FREG,FDR)
*
*      Galactic Centre and Bulge Forces
*      Dynamical Friction
*      External Forces on stars
*      R.Sp./ O.G. Nov. 2001
*
      INCLUDE 'common6.h'
      COMMON /GC/ XG0(3),VG0(3),AG0(3),ADOTG0(3),
     *      XG(3),VG(3),AG(3),ADOTG(3),AG2(3),AG3(3),
     *      XGI(3),VGI(3),AGI(3),ADOTGI(3),FTOT(3),FDTOT(3),
     *      ADF0(3),ADOTF0(3),ADF(3),ADOTF(3),
     *      XMBH,CC,CFAC,XCOUL,XMCL,TCEN0,STEPCC,EG0,TGCL,
     *      DEDFF,ECLS,ECLS0,ETAG,RGAL,VGAL,XMGAL,IMEM(NMAX)
      COMMON/POTENT/PHII(NMAX),PHIR(NMAX),PHIR1(NMAX)
      COMMON/XXX/IOLD,TOLD
      REAL*8 FREG(3),FDR(3),DER(3)
*
      IF(I.EQ.IOLD.AND.TIME.EQ.TOLD)THEN
      PRINT*,' Warning T,I,IOLD=',TIME,I,IOLD
      RETURN
      END IF
*
      IOLD = I
*
*      Relative positions and velocities to g.c. using pred. pos. of particle
*
*     Predictor for galactic centre
*     Note that XG,VG contain gravitation plus dynamical friction
*
*      Relative positions and velocities to g.c. using pred. pos.
*
       DO 1 K = 1,3
          XGI(K) = XG(K) + X(K,I)
          VGI(K) = VG(K) + XDOT(K,I)
 1     CONTINUE
*
       RIJ2 = XGI(1)**2 + XGI(2)**2 + XGI(3)**2
       DR2I = 1.0/RIJ2
       DRDV = XGI(1)*VGI(1) + XGI(2)*VGI(2) + XGI(3)*VGI(3)
*
*      First: central black hole
*
       DR3I = XMBH*DR2I*DSQRT(DR2I)
       DRDP = 3.0D0*DRDV*DR2I
*
*      Force and derivative at particle position due to g.c.
*
       DO 2 K = 1,3
          AGI(K) = -XGI(K)*DR3I
          ADOTGI(K) = -(VGI(K) - XGI(K)*DRDP)*DR3I
 2     CONTINUE
*
*      Second: bulge     
*
       DR3I = CC*DR2I
       DRDP = 2.0D0*DRDV*DR2I
*
*      Force and derivative at particle position due to g.c.
*
       DO 3 K = 1,3
          AGI(K) = AGI(K) - XGI(K)*DR3I
          ADOTGI(K) = ADOTGI(K) - (VGI(K) - XGI(K)*DRDP)*DR3I
 3     CONTINUE
*
*      Subtraction of apparent (pred.) force on cluster centre
*
       DO 4 K = 1,3
          AGI(K) = AGI(K) - AG(K)
          ADOTGI(K) = ADOTGI(K) - ADOTG(K)
 4     CONTINUE
*
       DO 5 K = 1,3
*      Subtract external forces (pred.) for non-cluster members
*         IF(IMEM(I).EQ.0)THEN
*            AGI(K) = AGI(K) - ADF(K)
*            ADOTGI(K) = ADOTGI(K) - ADOTF(K)
*         END IF
*
          FREG(K) = FREG(K) + AGI(K)
          FDR(K) = FDR(K) + ADOTGI(K)
*
 5     CONTINUE
*
*        DO 8 K = 1,3
*        DEDFF = DEDFF + BODY(I)*ADF(K)*VGI(K)*STEPR(I)
*8       CONTINUE
*
       RETURN
*
       END
