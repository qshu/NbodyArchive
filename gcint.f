      SUBROUTINE GCINT
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
*
*     Predictor for galactic centre
*     Note that XG,VG contain gravitation plus dynamical friction
*
      IFLAG = 0
      S = TIME - TCEN0
      S1 = 1.5D0*S
      S2 = 2.0D0*S
      DO 1 K = 1,3
      FTOT(K) = (AG0(K) + ADF0(K))/2.D0
      FDTOT(K) = (ADOTG0(K) + ADOTF0(K))/6.D0
      XG(K) = ((FDTOT(K)*S + FTOT(K))*S +VG0(K))*S + XG0(K)
      VG(K) = (FDTOT(K)*S1 + FTOT(K))*S2 + VG0(K)
 1    CONTINUE
*
*     High order prediction of galactic centre if no correction
*
*     IF (TIME.LT.TCEN0+STEPCC) THEN
*
*     DO 11 K = 1,3
*     XG(K) = XG(K) + S**2*(AG2(K) + S*AG3(K)/5.D0)*S**2/24.D0
*     VG(K) = VG(K) + S**2*(AG2(K) + S*AG3(K)/4.D0)*S/6.D0
* 11  CONTINUE
*
*     END IF
*
 555   CONTINUE
*
*      Force and derivative at cluster centre due to gravity of galactic centre
*
       RIJ2 = XG(1)**2 + XG(2)**2 + XG(3)**2
       DR2I = 1.0D0/RIJ2
       DRDV = XG(1)*VG(1) + XG(2)*VG(2) + XG(3)*VG(3)
*
*      First: central black hole
*
       DR3I = XMBH*DR2I*DSQRT(DR2I)
       DRDP = 3.0D0*DRDV*DR2I
*
       DO 2 K = 1,3
          AG(K) = -XG(K)*DR3I
          ADOTG(K) = -(VG(K) - XG(K)*DRDP)*DR3I
 2     CONTINUE
*
*      Second: bulge
*
       DR3I = CC*DR2I
       DRDP = 2.0D0*DRDV*DR2I
*
       DO 3 K = 1,3
          AG(K) = AG(K) - XG(K)*DR3I
          ADOTG(K) = ADOTG(K) - (VG(K) - XG(K)*DRDP)*DR3I
 3     CONTINUE
*
*      Force and derivative cluster centre due to dynamical friction
*
       RIJ2 = XG(1)**2 + XG(2)**2 + XG(3)**2
       VIJ2 = VG(1)**2 + VG(2)**2 + VG(3)**2
       DR2I = 1.0D0/RIJ2
       DV2I = 1.0D0/VIJ2
       DV3I = DV2I*DSQRT(DV2I)
       DFAC = -CFAC*XCOUL*XMCL*DR2I*DV3I*VIJ2
       DRDV = XG(1)*VG(1) + XG(2)*VG(2) + XG(3)*VG(3)
       DADV = AG(1)*VG(1) + AG(2)*VG(2) + AG(3)*VG(3)
*
       DRDP = DADV*DV2I
*
       DO 4 K = 1,3
*
       ADF(K) = DFAC*VG(K)
       ADOTF(K) = DFAC*(AG(K) - VG(K)*DADV*DV2I 
     *                                  - 2.D0*VG(K)*DRDV*DR2I)
 4     CONTINUE
*
       IF(IFLAG.EQ.1)THEN
       DO 55 K =1,3
       AG0(K) = AG(K)
       ADF0(K) = ADF(K)
       ADOTG0(K) = ADOTG(K)
       ADOTF0(K) = ADOTF(K)
 55    CONTINUE
       END IF
*     Correction of galactic centre
*
      IF (TIME.GE.TCEN0+STEPCC) THEN
*
      DTR = TIME - TCEN0
      DTSQ = DTR**2
      DT6 = 6.0/(DTR*DTSQ)
      DT2 = 2.0/DTSQ
      DTSQ12 = ONE12*DTSQ
      DTR13 = ONE3*DTR
*
        DO 6 K = 1,3
*
          DFR = AG0(K) + ADF0(K) - AG(K) - ADF(K)
          FDR0 = ADOTG(K) + ADOTF(K)
          FRD = ADOTG0(K) + ADOTF0(K)
      SUM = FRD + FDR0
      AT3 = 2.0*DFR + DTR*SUM
      BT2 = -3.0*DFR - DTR*(SUM + FRD)
*       Use here new variables for consistency in parallel execution (R.Sp.)
          XG0(K) = XG(K) + (0.6*AT3 + BT2)*DTSQ12
          VG0(K) = VG(K) + (0.75*AT3 + BT2)*DTR13
*
          AG2(K) = BT2*DT2
*         AG2(K) = (3.0*AT3 + BT2)*DT2
          AG3(K) = AT3*DT6
          XG(K) = XG0(K)
          VG(K) = VG0(K)
*
 6    CONTINUE
*
      TCEN0 = TIME
*
*       Obtain new regular integration step using composite expression.
      DO 90 K = 1,3
          W1(K) = ADOTG(K) + ADOTF(K)
          W2(K) = AG2(K)
          W3(K) = AG3(K)
  90  CONTINUE
*
      W0(4) = AG(1)**2 + AG(2)**2 + AG(3)**2 +
     *        ADF(1)**2 + ADF(2)**2 + ADF(3)**2
      W1(4) = W1(1)**2 + W1(2)**2 + W1(3)**2
      W2(4) = W2(1)**2 + W2(2)**2 + W2(3)**2
      W3(4) = W3(1)**2 + W3(2)**2 + W3(3)**2
*
*       Form new step by relative criterion (extra SQRT for large F3DOT).
      IF (W3(4).LT.1.0E+20) THEN
          W0(1) = (SQRT(W0(4)*W2(4)) + W1(4))/
     &                                       (SQRT(W1(4)*W3(4)) + W2(4))
      ELSE
          W0(1) = (SQRT(W0(4)*W2(4)) + W1(4))/
     &                                 (SQRT(W1(4))*SQRT(W3(4)) + W2(4))
      END IF
      W0(1) = ETAG*W0(1)
      TTMP = SQRT(W0(1))
*       Obtain regular force change using twice the predicted step.
      DTC = 2.0*TTMP
      S2 = 0.5*DTC
      S3 = ONE3*DTC
      W0(1) = 0.0
      DO 105 K = 1,3
          W0(2) = ((W3(K)*S3 + W2(K))*S2 + W1(K))*DTC
          W0(1) = W0(1) + W0(2)**2
  105 CONTINUE
*
*       See whether regular step can be increased by factor 2.
      IF (W0(1).LT.W0(4)) THEN
          TTMP = DTC
      END IF
*
*       Select discrete value (increased by 2, decreased by 2 or unchanged).
        IF (TTMP .GT. 2.0*STEPCC) THEN
            IF (DMOD(TIME,2.0*STEPCC) .EQ. 0.0D0) THEN
              TTMP = MIN(2.0*STEPCC,1.D0)
          ELSE
              TTMP = STEPCC
          END IF
      ELSE IF (TTMP .LT. STEPCC) THEN
          TTMP = 0.5*STEPCC
      ELSE
          TTMP = STEPCC
      END IF
*
      STEPCC = TTMP
*
      XXG = DSQRT(XG(1)**2 + XG(2)**2 + XG(3)**2)
      XVG = DSQRT(VG(1)**2 + VG(2)**2 + VG(3)**2)
      XAG = DSQRT(AG(1)**2 + AG(2)**2 + AG(3)**2)
      XXGI = DSQRT(XGI(1)**2 + XGI(2)**2 + XGI(3)**2)
      XVGI = DSQRT(VGI(1)**2 + VGI(2)**2 + VGI(3)**2)
      XAGI = DSQRT(AGI(1)**2 + AGI(2)**2 + AGI(3)**2)
      XADOTG = DSQRT(ADOTG(1)**2 + ADOTG(2)**2 + ADOTG(3)**2)
      XADF = DSQRT(ADF(1)**2 + ADF(2)**2 + ADF(3)**2)
      XADOTF = DSQRT(ADOTF(1)**2 + ADOTF(2)**2 + ADOTF(3)**2)
*
      NCEN = NCEN + 1
      if(rank.eq.0)then
      WRITE(77,177)TIME,TCEN0,STEPCC,XXG,XVG,XAG,XADOTG,XADF,XADOTF,
     *    NCEN
 177  FORMAT(1X,1P,9D13.3,I10)
      WRITE(78,177)TIME,(AG(K),K=1,3),(ADOTG(K),K=1,3)
      WRITE(79,177)TIME,(XG(K),K=1,3),(VG(K),K=1,3)
      CALL FLUSH(77)
      CALL FLUSH(78)
      CALL FLUSH(79)
      end if
*       Recompute AG,ADF,ADOTG,ADOTF with corrected values
*
      IFLAG = 1
*
      GO TO 555
*
      END IF
*
 100  CONTINUE
*
      RETURN
*
      END
