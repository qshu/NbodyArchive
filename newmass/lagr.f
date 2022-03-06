      SUBROUTINE LAGR
*
*
*       Lagrangian radii.
*       -----------------
*
      INCLUDE 'common6.h'
      INTEGER NSHELL(10)
      REAL*8  R2(NMAX)
      REAL*8  C(3),FLAGR(11),RLAGR(11),AVMASS(11)
*       Lagrangian radii at 1,2,5,7.5,10,15,20,50,75,85 % of total mass.
      DATA FLAGR  /1.E-2,2.E-2,5.E-2,7.5E-2,1.E-1,1.5E-1,2.E-1,
     &             5.E-1,7.5E-1,8.5E-1,1.0E0/
*
*
*       Set square radii of all single particles & c.m. bodies.
      NP = 0
      DO 10 I = 1,N
          NP = NP + 1
          R2(NP) = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
          JLIST(NP) = I
   10 CONTINUE
*
*       Sort square distances.
      CALL SORT1(NP,R2,JLIST)
*
*  Determine the Lagrangian radii for specified mass fractions.
*     RLAGR = Lagrangian radius
*     AVMASS = average mass of a spherical shell with radius R2(I)
*     NSHELL = particle counter within a shell
      ZM = 0.0D0
      ZMH = 0.5D0*ZMASS
      I = 0
*
      XMR1T = 0.D0
      XMR2T = 0.D0
      XMINT = 0.D0
      XN1T = 0.D0
      XN2T = 0.D0
      RI = 0.D0
      RHLAST = 0.D0
      R1LAST = 0.D0
      R2LAST = 0.D0
      RH1LST = 0.D0
      RH2LST = 0.D0
      II1 = 0
      II2 = 0
      TWOPI = 8.0D0*DATAN(1.D0)
      SX = 1.5*TWOPI/16.0
*
      DO 100 J = 1, 11
         AVMASS(J) = 0.D0
         NSHELL(J) = 0
 20      ILAST = I
         RLAST = RI
         IF(IM.GT.0)RHLAST=RHOPL(IM)
         I = I + 1
         IM = JLIST(I)
         ZM = ZM + BODY(IM)
         AVMASS(J) = AVMASS(J) + BODY(IM)
         RI = SQRT(X(1,IM)**2 + X(2,IM)**2 + X(3,IM)**2)
         RRR = RI/SX
         XMRPL(IM) = RRR**3/(1.D0+RRR*RRR)**1.5D0
         RHOPL(IM) = 3.D0/(2.D0*TWOPI*SX**3)/(1.D0+RRR*RRR)**2.5D0
         DVOL = 2.D0*TWOPI/3.D0*(RI**3-RLAST**3)
         XMINT = XMINT + XFUNC(IM)*(RHOPL(IM)+RHLAST)*DVOL/2.D0
*
         IF(IND(IM).EQ.2)THEN
         XMR2T = XMR2T + BODY(IM)
         DVOL = 2.D0*TWOPI/3.D0*(RI**3-R2LAST**3)
         XN2T = XN2T + RHOC(IM)*(RHOPL(IM)+RH2LST)*DVOL/2.D0
         II2 = II2 + 1
         R2LAST = RI
         RH2LST = RHOPL(IM)
         END IF
*
         IF(IND(IM).EQ.1)THEN
         XMR1T = XMR1T + BODY(IM)
         DVOL = 2.D0*TWOPI/3.D0*(RI**3-R1LAST**3)
         XN1T = XN1T + RHOC(IM)*(RHOPL(IM)+RH1LST)*DVOL/2.D0
         II1 = II1 + 1
         R1LAST = RI
         RH1LST = RHOPL(IM)
         END IF
*
         XMR1(IM) = XMR1T
         XMR2(IM) = XMR2T
         NSHELL(J) = NSHELL(J) + 1
         IF (ZM .LT. FLAGR(J)*ZMASS) GOTO 20
         RLAGR(J) = SQRT(R2(I))
         AVMASS(J) = AVMASS(J)/NSHELL(J)
 100  CONTINUE
*
         XN1T = XN1T*XM1TOT/XM1
         XN2T = XN2T*XM2TOT/XM2
*
      PRINT*,' LAGR: XMR1,2T,XMINT,XN1,2=',XMR1T,XMR2T,II1,II2,XMINT,
     *    XN1T,XN2T
*       Replace approximate half-mass radius by actual value.
      RSCALE = SQRT(R2(I))
*       Return for second pass to get correct total mass
*
      REWIND 13
*       Check output options (line printer or unit 7 or both).
          WRITE (6,30)  (FLAGR(K),K=1,11)
          WRITE (6,40)  (RLAGR(K),K=1,11)
          WRITE (6,50)  (AVMASS(K),K=1,11)
   30     FORMAT (/,3X,'FLGR:  ',1P,11(1X,D12.5))
   40     FORMAT (/,3X,'LAGR:  ',1P,11(1X,D12.5))
   50     FORMAT (/,3X,'AVMA:  ',1P,11(1X,D12.5))
*
          WRITE (13,30)  (FLAGR(K),K=1,11)
          WRITE (13,40)  (RLAGR(K),K=1,11)
          WRITE (13,50)  (AVMASS(K),K=1,11)
*
      RETURN
*
      END
