      SUBROUTINE KSREL(IPAIR,I1,XI,VI,FRELP,FRELD,VWHOLE,APPLY)
*
*
*       Relativistic perturbation on KS pair.
*       -------------------------------------
*
*       Coded by Pau Amaro-Seoane 2006.
*       -------------------------------
*
*       3PN,3.5PN and Spin terms added by Patrick Brem 2011.
*       ----------------------------------------------------
* 
      INCLUDE 'common6.h'
C      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
C     &                LISTC(LMAX)
 	  
       REAL*8 XI(6),VI(6),FRELP(3),FRELD(3),PI2
       REAL*8 RIJ, RIJ2,A(3),AD(3),KSAK,KSBK,ADK,BDK
       REAL*8 NX,NY,NZ,POWER1,POWER2,RCSCH
       REAL*8 KSV(3),KSX(3),KSA1(3),KSA2(3),KSV1(3),KSV2(3)
       REAL*8 CL, CL3, CL2, CL4, CL5, CL6, CL7
       REAL*8 A1,A2,A2_5,A3,A3_5,B1,B2,B2_5,B3,B3_5,ETA,M,RP,MOR,RPP
       REAL*8 A1D,A2D,A2_5D,A3D,A3_5D,B1D,B2D,B2_5D,B3D,B3_5D,VA
       REAL*8 KSIMPMOM2,KSECC,KSH,V1_V22,V12,V22,VWHOLE
       REAL*8 EKIN, EPOT,ENEWT
       REAL*8 XCOMPIMP,YCOMPIMP,ZCOMPIMP,FRP,GREL
* (c) patrick
       REAL*8 KSS(3),KSSIG(3),NVEC(3),AT(3)
       REAL*8 NCV(3),NCS(3),NCSIG(3),VCS(3),VCSIG(3)
       REAL*8 SDNCV,SIGDNCV,NDV,DM,S1(3),S2(3)
       REAL*8 C1_5(3),C2(3),C2_5(3),C1_5D(3),C2_5D(3),C2D(3)
       REAL*8 VDS,VDSIG,NDS,NDSIG,S1DLU,S2DLU
       REAL*8 SU1(3),SU2(3),SU(3),SV(3),SV1(3),SV2(3)
       REAL*8 XX1(3),XX2(3),L1(3),L2(3),L(3),LU(3) 
       REAL*8 VV1(3),VV2(3),LABS,SS1(3),SS2(3)
       REAL*8 NVDOT, NDOT(3), SNVDOT, SIGNVDOT
       REAL*8 NDOTCV(3),NCA(3),NDOTCS(3),NCSU(3)
       REAL*8 NDOTCSIG(3),NCSV(3),ACS(3),VCSU(3)
       REAL*8 ACSIG(3),VCSV(3),NSDOT,NSIGDOT,VSDOT,VSIGDOT
       REAL*8 XS(3),XA(3),NXA,NXS,XS2,XA2,XAD(3),XSD(3)
       REAL*8 NXADOT, NXSDOT, ESPIN, ANEWTON, SUM
       REAL*8 RATIO, ABSL, ABSA1, ABSA2, COSA, COSB, COSG
*  For energy correction
       REAL*8 AS(3),AP(3),POWER1S,POWER1P,POWER2S,POWER2P
       REAL*8 ADS(3),ADP(3),A25POWER(3),A25DPOWER(3)
       REAL*8 AT1(3),AT2(3),JVEC(3),JABS, SQM(3), CHI12,CHI22
*     variables such as NVEC and RIJ are recycled in the perturber sum at 101

       INTEGER I1,I2,K1,APPLY,IPAIR

C       WRITE (*,*) 'KSREL'
       
       PI2 = 9.869604401089359
        I2 = I1 + 1
*       calculate physical spin from 0<SPIN<1 in the loop
* Define relative coordinates and velocity components.
      DO 2 K = 1,3
        S1(K) = SPN(K,I1)*BODY(I1)*BODY(I1)/CLIGHT
        S2(K) = SPN(K,I2)*BODY(I2)*BODY(I2)/CLIGHT
	KSX(K) = XI(K)-XI(K+3)	
	KSV(K) = VI(K)-VI(K+3)
	KSS(K) = S1(K)+S2(K)
        KSSIG(K) = (BODY(I1)+BODY(I2))*(S2(K)/BODY(I2)-
     &                  S1(K)/BODY(I1))
        XS(K)  = 0.5*(SPN(K,I1)+SPN(K,I2))
        XA(K)  = 0.5*(SPN(K,I1)-SPN(K,I2))
    2 CONTINUE

*diagnostics
C      IF (I1.GT.N.OR.I2.GT.N) WRITE (*,*) 'ALARM!'

      RFLAG(IPAIR) = 1
	      
* Set speed of light CL and its powers.
        CL = CLIGHT
        CL2 = CL*CL
        CL3 = CL2*CL
        CL4 = CL2*CL2
        CL5 = CL4*CL
	CL6 = CL5*CL
	CL7 = CL6*CL
	  
        M = BODY(I1)+BODY(I2)
        RCSCH = 2.0*M/CL2
        ETA = BODY(I1)*BODY(I2)/(M*M)
        DM = BODY(I1)-BODY(I2)
	
        RIJ2 = KSX(1)**2+KSX(2)**2+KSX(3)**2
        RIJ = SQRT(RIJ2)
        MOR = M/RIJ
* lighter body over heavier body
        RATIO = BODY(I2)/BODY(I1)
        IF (BODY(I2).GT.BODY(I1)) THEN RATIO = 1/RATIO END IF
	  
        V1_V22 = KSV(1)*KSV(1)+KSV(2)*KSV(2)+KSV(3)*KSV(3)
        VWHOLE = SQRT(V1_V22)
	  
        NX = KSX(1)/RIJ
        NY = KSX(2)/RIJ
        NZ = KSX(3)/RIJ

        NVEC(1) = NX
        NVEC(2) = NY
        NVEC(3) = NZ

        CHI12 = SPN(1,I1)**2+SPN(2,I1)**2+SPN(3,I1)**2
        CHI22 = SPN(1,I2)**2+SPN(2,I2)**2+SPN(3,I2)**2

* Dot and Crossproducts needed for Spin-Orbit contribution
        NDV = NX*KSV(1) + NY*KSV(2) + NZ*KSV(3)
        VDS = KSV(1)*KSS(1)+KSV(2)*KSS(2)+KSV(3)*KSS(3)
        VDSIG = KSV(1)*KSSIG(1)+KSV(2)*KSSIG(2)+
     &       KSV(3)*KSSIG(3)
        NDS = NX*KSS(1)+NY*KSS(2)+NZ*KSS(3)
        NDSIG = NX*KSSIG(1)+NY*KSSIG(2)+NZ*KSSIG(3)
        CALL CROSSP(NVEC,KSV,NCV)
        CALL CROSSP(NVEC,KSS,NCS)
        CALL CROSSP(NVEC,KSSIG,NCSIG)
        CALL CROSSP(KSV,KSS,VCS)
        CALL CROSSP(KSV,KSSIG,VCSIG)
        SDNCV = KSS(1)*NCV(1)+KSS(2)*NCV(2)+KSS(3)*NCV(3)
        SIGDNCV = KSSIG(1)*NCV(1)+KSSIG(2)*NCV(2)+KSSIG(3)*NCV(3)
        CALL CROSSP(KSX,KSV,L)
        L(1) = L(1)*M*ETA
        L(2) = L(2)*M*ETA
        L(3) = L(3)*M*ETA
        LABS = SQRT(L(1)*L(1)+L(2)*L(2)+L(3)*L(3))
        LU(1) = L(1)/LABS
        LU(2) = L(2)/LABS
        LU(3) = L(3)/LABS 
        S1DLU = S1(1)*LU(1)+S1(2)*LU(2)+S1(3)*LU(3)
        S2DLU = S2(1)*LU(1)+S2(2)*LU(2)+S2(3)*LU(3)
        NXA = NVEC(1)*XA(1)+NVEC(2)*XA(2)+NVEC(3)*XA(3)
        NXS = NVEC(1)*XS(1)+NVEC(2)*XS(2)+NVEC(3)*XS(3)
        XA2 = XA(1)*XA(1)+XA(2)*XA(2)+XA(3)*XA(3)
        XS2 = XS(1)*XS(1)+XS(2)*XS(2)+XS(3)*XS(3)

   
* RP = RDOT
        RP = (KSX(1)*KSV(1)+KSX(2)*KSV(2)+KSX(3)*KSV(3))/RIJ

        IF (APPLY.EQ.0.OR.APPLY.EQ.1) THEN
************1PN*********
        A1 = 2.0*(2.0+ETA)*MOR-(1.0+3.0*ETA)*V1_V22 +1.5*ETA*RP*RP
        B1 = 2.0*(2.0-ETA)*RP	
************1PN*********
************2PN*********
	A2 = -0.75*(12.0+29.0*ETA)*MOR*MOR-
     &       ETA*(3.0-4.0*ETA)*V1_V22*V1_V22-
     &       1.875*ETA*(1.0-3.0*ETA)*RP*RP*RP*RP+
     &       0.5*ETA*(13.0-4.0*ETA)*MOR*V1_V22+
     &       (2.0+25.0*ETA+2.0*ETA*ETA)*MOR*RP*RP+
     &       1.5*ETA*(3.0-4.0*ETA)*V1_V22*RP*RP
        B2 = -0.5*RP*((4.0+41.0*ETA+8.0*ETA*ETA)*MOR-
     &       ETA*(15.0+4.0*ETA)*V1_V22+3.0*ETA*(3.0+2.0*ETA)*RP*RP)
************2PN*********
************2.5PN********
	 A2_5 = 1.6*ETA*MOR*RP*(17.0*MOR/3.0+3.0*V1_V22)
	 B2_5 = -1.6*ETA*MOR*(3.0*MOR+V1_V22)
************2.5PN********
************3PN********
	 A3 = MOR*MOR*MOR*(16.0+(1399.0/12.0-41.0*PI2/16.0)*ETA+
     &71.0*ETA*ETA/2.0)+ETA*(20827.0/840.0+123.0*PI2/64.0-ETA*ETA)
     &*MOR*MOR*V1_V22-(1.0+(22717.0/168.0+615.0*PI2/64.0)*ETA+
     &11.0*ETA*ETA/8.0-7.0*ETA*ETA*ETA)*MOR*MOR*RP*RP-
     &0.25*ETA*(11.0-49.0*ETA+52.0*ETA*ETA)*V1_V22*V1_V22*V1_V22+
     &35.0*ETA*(1.0-5.0*ETA+5.0*ETA*ETA)*RP*RP*RP*RP*RP*RP/16.0-
     &0.25*ETA*(75.0+32.0*ETA-40.0*ETA*ETA)*MOR*V1_V22*V1_V22-
     &0.5*ETA*(158.0-69.0*ETA-60.0*ETA*ETA)*MOR*RP*RP*RP*RP+
     &ETA*(121.0-16.0*ETA-20.0*ETA*ETA)*MOR*V1_V22*RP*RP+
     &3.0*ETA*(20.0-79.0*ETA+60.0*ETA*ETA)*V1_V22*V1_V22*RP*RP/8.0-
     &15.0*ETA*(4.0-18.0*ETA+17.0*ETA*ETA)*V1_V22*RP*RP*RP*RP/8.0

         B3 = RP*((4.0+(5849.0/840.0+123.0*PI2/32.0)*ETA-25.0*ETA*ETA-
     &8.0*ETA*ETA*ETA)*MOR*MOR+ETA*(65.0-152.0*ETA-48.0*ETA*ETA)*
     &V1_V22*V1_V22/8.0+15.0*ETA*(3.0-8.0*ETA-2.0*ETA*ETA)*RP*RP*RP*
     &RP/8.0+
     &ETA*(15.0+27.0*ETA+10.0*ETA*ETA)*MOR*V1_V22-ETA*(329.0+177.0*ETA+
     &108.0*ETA*ETA)*MOR*RP*RP/6.0-
     &3.0*ETA*(16.0-37.0*ETA-16.0*ETA*ETA)*V1_V22*RP*RP/4.0)
************3PN********
************3.5PN********
        A3_5 = MOR*ETA*(V1_V22*V1_V22*(-366./35.-12.*ETA)+
     &       V1_V22*RP*RP*(114.+12.*ETA)-112.*RP*RP*RP*RP+
     &       MOR*(V1_V22*(-692./35.+724.*ETA/15.)+RP*RP*
     &       (-294./5.-376.*ETA/5.)+MOR*(-3956./35.-184.*ETA/5.)))
 
         B3_5 = 8.0*ETA*MOR*((1325.0+546.0*ETA)*MOR*MOR/42.0+
     &(313.0+42.0*ETA)*V1_V22*V1_V22/28.0+75.0*RP*RP*RP*RP
     &-(205.0+777.0*ETA)*MOR*V1_V22/42.0+(205.0+
     &424.0*ETA)*MOR*RP*RP/12.0-3.0*(113.0+2.0*ETA)*V1_V22*RP*RP/4.0)/5.	
************3.5PN********
*******1.5PN SPIN-ORBIT****
      DO 1 K = 1,3
        C1_5(K) = (NVEC(K)*(12.*SDNCV+6.*DM*SIGDNCV/M)+
     &       9.*NDV*NCS(K)+3.*DM*NDV*NCSIG(K)/M -
     &       7.*VCS(K)-3.*DM*VCSIG(K)/M)/(RIJ2*RIJ)
*******1.5PN SPIN-ORBIT****
*******2 PN SPIN-SPIN******
        C2(K) = -MOR*MOR*MOR/RIJ*3.*ETA*(NVEC(K)*(XS2-XA2-
     &       5.*NXS*NXS+5.*NXA*NXA)+2.*(XS(K)*NXS-XA(K)*NXA))
*******2 PN SPIN-SPIN******
*******2.5PN SPIN-ORBIT****
        C2_5(K) = (NVEC(K)*(SDNCV*(-30.*ETA*NDV*NDV+
     &       24.*ETA*V1_V22-MOR*(38.+25.*ETA))+DM/M*SIGDNCV*
     &       (-15.*ETA*NDV*NDV+12.*ETA*V1_V22-MOR*(18.+14.5*
     &       ETA)))+NDV*KSV(K)*(SDNCV*(-9.+9.*ETA)+DM/M*
     &       SIGDNCV*(-3.+6.*ETA))+NCV(K)*(NDV*VDS*(-3.+3.*ETA)-
     &       8.*MOR*ETA*NDS-DM/M*(4.*MOR*ETA*NDSIG+3.*NDV*
     &       VDSIG))+NDV*NCS(K)*(-22.5*ETA*NDV*NDV+21.*ETA*
     &       V1_V22-MOR*(25.+15.*ETA))+DM/M*NDV*NCSIG(K)*(-15.*
     &       ETA*NDV*NDV+12.*ETA*V1_V22-MOR*(9.+8.5*ETA))+
     &       VCS(K)*(16.5*ETA*NDV*NDV+MOR*(21.+9.*ETA)-14.*ETA*
     &       V1_V22)+DM/M*VCSIG(K)*(9.*ETA*NDV*NDV-7.*ETA*V1_V22+
     &       MOR*(9.+4.5*ETA)))/(RIJ2*RIJ)
*******2.5PN SPIN-ORBIT****
*******2 PN SPIN QUADRUPOLE
        IF (CHI12.GT.0.AND.CHI22.GT.0.) THEN
        SQM(K) = -3.*MOR*MOR*MOR/RIJ*ETA/2.*(CHI12/RATIO*((1.-
     & 5.*(NVEC(1)*SPN(1,I1)+NVEC(2)*SPN(2,I1)+NVEC(3)*
     & SPN(3,I1))**2/CHI12)*NVEC(K)+2.*(NVEC(1)*SPN(1,I1)+
     & NVEC(2)*SPN(2,I1)+NVEC(3)*
     & SPN(3,I1))/CHI12*SPN(K,I1))+
     & RATIO*CHI22*((1.-
     & 5.*(NVEC(1)*SPN(1,I2)+NVEC(2)*SPN(2,I2)+NVEC(3)*
     & SPN(3,I2))**2/CHI22)*NVEC(K)+2.*(NVEC(1)*SPN(1,I2)+
     & NVEC(2)*SPN(2,I2)+NVEC(3)*
     & SPN(3,I2))/CHI12*SPN(K,I2)))
        ELSE
           SQM(K) = 0.0
        END IF
*******2 PN SPIN QUADRUPOLE
    1 CONTINUE

        KSAK = A1/CL2+A2/CL4+A2_5/CL5+A3/CL6+A3_5/CL7
	KSBK = B1/CL2+B2/CL4+B2_5/CL5+B3/CL6+B3_5/CL7
*     GAMREL for switch-on criterion
	GAMREL(1,IPAIR) = SQRT((A1*NX + B1*KSV(1))**2+(A1*NY+B1
     &  *KSV(2))**2+(A1*NZ+B1*KSV(3))**2)/CL2
	GAMREL(2,IPAIR) = SQRT((A2*NX + B2*KSV(1))**2+(A2*NY
     &  +B2*KSV(2))**2+(A2*NZ+B2*KSV(3))**2)/CL4
	GAMREL(3,IPAIR) = SQRT((A2_5*NX + B2_5*KSV(1))**2+(A2_5*NY
     &  +B2_5*KSV(2))**2+(A2_5*NZ+B2_5*KSV(3))**2)/CL5
	GAMREL(4,IPAIR) = SQRT((A3*NX + B3*KSV(1))**2+(A3*NY+B3*KSV(2)
     &  )**2+(A3*NZ+B3*KSV(3))**2)/CL6
	GAMREL(5,IPAIR) = SQRT((A3_5*NX + B3_5*KSV(1))**2+(A3_5*NY
     &  +B3_5*KSV(2))**2+(A3_5*NZ+B3_5*KSV(3))**2)/CL7
	GAMREL(6,IPAIR) = SQRT(C1_5(1)*C1_5(1)+C1_5(2)*C1_5(2)+
     &  C1_5(3)*C1_5(3))/CL2/MOR*RIJ
	GAMREL(7,IPAIR) = SQRT(C2(1)*C2(1)+C2(2)*C2(2)+C2(3)*C2(3)
     &  )/CL4/MOR*RIJ
	GAMREL(8,IPAIR) = SQRT(C2_5(1)*C2_5(1)+C2_5(2)*C2_5(2)+
     &  C2_5(3)*C2_5(3))/CL4/MOR*RIJ

        A(1) = MOR*(KSAK*NX + KSBK*KSV(1))/RIJ
        A(1) = A(1)+ C1_5(1)/CL2 + (C2(1)+C2_5(1)+SQM(1))/CL4

        A(2) = MOR*(KSAK*NY + KSBK*KSV(2))/RIJ
        A(2) = A(2)+ C1_5(2)/CL2 + (C2(2)+C2_5(2)+SQM(2))/CL4

        A(3) = MOR*(KSAK*NZ + KSBK*KSV(3))/RIJ
        A(3) = A(3)+ C1_5(3)/CL2 + (C2(3)+C2_5(3)+SQM(3))/CL4

        DO 345 K = 1,3
           FRELP(K) = A(K)
 345    CONTINUE

        END IF
        IF (APPLY.EQ.1.OR.APPLY.EQ.2) THEN

*     Total acceleration (+ newtonian)
        AT(1) = FRELP(1) - MOR*NX/RIJ
        AT(2) = FRELP(2) - MOR*NY/RIJ
        AT(3) = FRELP(3) - MOR*NZ/RIJ
        
        RPP = V1_V22/RIJ + AT(1)*NX+AT(2)*NY + AT(3)*NZ - RP*RP/RIJ
*** VA = KSV*VDOT needed when differentiating KSV**2
        VA = AT(1)*KSV(1) + AT(2)*KSV(2) + AT(3)*KSV(3)

**For Spin
      DO 4 K = 1,3
        NDOT(K) = (KSV(K)-NVEC(K)*RP)/RIJ
    4 CONTINUE
        NVDOT = NDOT(1)*KSV(1)+NDOT(2)*KSV(2)+NDOT(3)*KSV(3)+
     &       NVEC(1)*AT(1)+NVEC(2)*AT(2)+NVEC(3)*AT(3)

        CALL CROSSP(NDOT,KSV,NDOTCV)
        CALL CROSSP(NVEC,AT,NCA)


** FDOT'S

************1PN*********
	A1D = -2.0*(2.0+ETA)*MOR*RP/RIJ - 2.0*(1.0+3.0*ETA)*VA +
     &	       3.0*ETA*RP*RPP
        B1D = 2.0*(2.0-ETA)*RPP 
************1PN*********
************2PN*********
        A2D = 1.5*(12.0+29.0*ETA)*MOR*MOR*RP/RIJ - 
     &       ETA*(3.0-4.0*ETA)*4.0*V1_V22*VA - 7.5*ETA*(1.0-3.0*ETA)*
     &       RPP -
     &       0.5*ETA*(13.0-4.0*ETA)*MOR*RP*V1_V22/RIJ+
     &       ETA*(13.0-4.0*ETA)*MOR*VA -
     &       (2.0+25.0*ETA+2.0*ETA*ETA)*MOR*RP*RP*RP/RIJ+
     &       2.0*(2.0+25.0*ETA+2.0*ETA*ETA)*MOR*RP*RPP + 
     &       3.0*ETA*(3.0-4.0*ETA)*VA*RP*RP +
     &       3.0*ETA*(3.0-4.0*ETA)*V1_V22*RP*RPP

        B2D = -0.5*RPP*((4.0+41.0*ETA+8.0*ETA*ETA)*MOR -
     &       ETA*(15.0+4.0*ETA)*V1_V22+3.0*ETA*(3.0+2.0*ETA)*RP*RP) -
     &       0.5*RP*(-(4.0+41.0*ETA+8.0*ETA*ETA)*MOR*RP/RIJ - 
     &       2.0*ETA*(15.0+4.0*ETA)*VA + 6.0*ETA*(3.0+2.0*ETA)*RP*RPP)
************2PN*********
************2.5PN********
        A2_5D = -1.6*ETA*MOR*RP*RP*(17.0/3.0*MOR+3.0*V1_V22)/RIJ +
     &       1.6*ETA*MOR*RPP*(17.0/3.0*MOR+3.0*V1_V22)+
     &       1.6*ETA*MOR*RP*(-17.0*MOR*RP/3.0/RIJ+6.0*VA)
	 
	B2_5D = 1.6*ETA*MOR*RP*(3.0*MOR+V1_V22)/RIJ - 
     &       1.6*ETA*MOR*(-3.0*MOR*RP/RIJ+2.0*VA) 
************2.5PN********
************3PN********
        A3D = 6.*ETA*RP*RP*RP*RP*RP*RPP*(35.-175.*ETA+
     &       175.*ETA*ETA)/16. + ETA*(4.*RP*RP*RP*RPP*V1_V22 + 
     &       2.*RP*RP*RP*RP*VA)*(-15.+135.*ETA/2.-255.*ETA*ETA/4.)/2. +
     &       ETA*(2.*RP*RPP*V1_V22*V1_V22+4.*RP*RP*V1_V22*VA)/2.*
     &       (15.-237.*ETA/2.+45.*ETA*ETA) + 6.*V1_V22*V1_V22*VA*ETA*
     &       (-11./4.-49.*ETA/4.-13.*ETA*ETA) + MOR*(4.*RP*RP*RP*RPP*
     &       ETA*(-79.+69./2.*ETA+30.*ETA*ETA) + ETA*(2.*RP*RPP*V1_V22+
     &       2.*RP*RP*VA)*(121.-16.*ETA-20.*ETA*ETA)+4.*V1_V22*VA*ETA*
     &       (-75./4.-8.*ETA+10.*ETA*ETA)) - MOR*RP*((-79.+69.*ETA/
     &       2.+30.*ETA*ETA)*RP*RP*RP*RP*ETA+ETA*RP*RP*V1_V22*(121.-
     &       16.*ETA-20.*ETA*ETA)+ETA*V1_V22*V1_V22*(-75./4.-8.*ETA+10.*
     &       ETA*ETA))/RIJ - 2.*MOR*MOR*RP*(RP*RP*((-1.-615.*PI2*ETA/
     &       64.)-22717.*ETA/168.-11.*ETA*ETA/8.+7.*ETA*ETA*ETA)+ETA*
     &       V1_V22*((20827./840.+123.*PI2/64.)-ETA*ETA))/RIJ + MOR*MOR*
     &       (2.*RP*RPP*((-1.-615*PI2*ETA/64.)-22717.*ETA/168.-11.*
     &       ETA*ETA/8.+7*ETA*ETA*ETA)+2.*ETA*VA*((20827./840. +
     &       123.*PI2/64.)-ETA*ETA)) - 3.*MOR*MOR*MOR*RP*(16.+
     &       (1399./12.-41.*PI2/16.)*ETA+71.*ETA*ETA/2.)/RIJ

        B3D = 75.*RP*RP*RP*RP*RPP*ETA*(3./8.-ETA-.25*ETA*ETA)+
     &       ETA*(3.*RP*RP*RPP*V1_V22+2.*RP*RP*RP*VA)*(-12.+111.*
     &       ETA/4.+12.*ETA*ETA)+ETA*(RPP*V1_V22*V1_V22+4.*RP*V1_V22*
     &       VA)*(65./8.-19.*ETA-6.*ETA*ETA)-MOR*RP*(RP*RP*RP*
     &       ETA*(-329./6.-59.*ETA/2.-18.*ETA*ETA)+RP*V1_V22*ETA*
     &       (15.+27.*ETA+10.*ETA*ETA))/RIJ+MOR*(3.*RP*RP*RPP*ETA*
     &       (-329./6.-59.*ETA/2.-18.*ETA*ETA)+ETA*(RPP*V1_V22+
     &       2.*RP*VA)*(15.+27.*ETA+10.*ETA*ETA))-2.*MOR*MOR*RP*
     &       (RP*((4.+123.*PI2*ETA/32.)+5849.*ETA/840.-25.*ETA*
     &       ETA-8.*ETA*ETA*ETA))/RIJ+MOR*MOR*(RPP*((4.+123.*PI2*
     &       ETA/32.)+5849./840.*ETA-25.*ETA*ETA-8.*ETA*ETA*ETA))
************3PN********
************3.5PN********
        A3_5D = MOR*ETA*(-RP*(V1_V22*V1_V22*(-366./35.-12.*ETA)+
     &       V1_V22*RP*RP*(114.+12.*ETA)+RP*RP*RP*RP*(-112.))/RIJ+
     &       4*V1_V22*VA*(-366./35.-12.*ETA)+2.*(VA*RP*RP+RP*RPP*
     &       V1_V22)*(114.+12.*ETA)+4.*RP*RP*RP*RPP*(-112.)+MOR*
     &       (2*VA*(-692./35.+724.*ETA/15.)+2.*RP*RPP*(-294./5.-
     &       376.*ETA/5.)-2.*RP*(V1_V22*(-692./35.+724.*ETA/15.)+
     &       RP*RP*(-294./5.-376.*ETA/5.))/RIJ-3.*MOR*RP*(-3956./
     &       35.-184.*ETA/5.)/RIJ))

        B3_5D = MOR*ETA*(4.*V1_V22*VA*(626./35.+12.*ETA/5.)+
     &       2.*(VA*RP*RP+V1_V22*RP*RPP)*(-678./5.-12.*ETA/5.)+
     &       4.*RP*RP*RP*RPP*120.-RP*(V1_V22*V1_V22*(626./35.+
     &       12.*ETA/5.)+V1_V22*RP*RP*(-678./5.-12.*ETA/5.)+
     &       120.*RP*RP*RP*RP)/RIJ+MOR*(2.*VA*(-164./21.-148.*ETA/
     &       5.)+2*RP*RPP*(82./3.+848.*ETA/15.)-2.*RP*(V1_V22*
     &       (-164./21-148.*ETA/5.)+RP*RP*(82./3.+848.*ETA/15.))/
     &       RIJ-3.*MOR*RP*(1060./21.+104.*ETA/5.)/RIJ))

************3.5PN********

***********************

*       Simple Spin Euler evolution
      DO 3 K = 1,3
**********1PN**********
          SU1(K) = MOR*ETA*(NVEC(K)*(-4.*VDS-2.*DM/M*VDSIG)+
     &         KSV(K)*(3.*NDS+DM/M*NDSIG)+NDV*(2.*KSS(K)+DM/M*KSSIG(K)))
     &         /RIJ
          
          SV1(K) = MOR*(NVEC(K)*(VDSIG*(-2.+4.*ETA)-2.*DM/M*VDS)+
     &         KSV(K)*(NDSIG*(1.-ETA)+DM/M*NDS)+NDV*(KSSIG(K)*(1.-
     &         2.*ETA)+
     &         DM/M*KSS(K)))/RIJ
**********1PN**********
****1.5PN SPIN SPIN****
        SS1(K) = 0.5*(L(K)*(4.+3.*(BODY(I2)/BODY(I1)))+
     &       (S2(K)-3.*S2DLU*LU(K)))/(RIJ2*RIJ)

        SS2(K) = 0.5*(L(K)*(4.+3.*(BODY(I1)/BODY(I2)))+
     &       (S1(K)-3.*S1DLU*LU(K)))/(RIJ2*RIJ)

****1.5PN SPIN SPIN****
**********2PN**********

          SU2(K) = MOR*ETA/RIJ*(NVEC(K)*(VDS*(-2.*V1_V22+3.*NDV*NDV-
     &         6.*ETA*NDV*NDV+7.*MOR-8.*ETA*MOR)-14.*MOR*NDS*NDV+
     &         DM/M*VDSIG*ETA*(-3.*NDV*NDV-4.*MOR)+DM/M*MOR*NDSIG*NDV*
     &         (2.-ETA/2.))+KSV(K)*(NDS*(2.*V1_V22-4.*ETA*V1_V22-3.*NDV*
     &         NDV+7.5*ETA*NDV*NDV+4.*MOR-6.*ETA*MOR)+VDS*NDV*(2.-
     &         6.*ETA)+
     &         DM/M*NDSIG*(-1.5*ETA*V1_V22+3.*ETA*NDV*NDV-MOR-3.5*ETA*
     &         MOR)-3.*DM/M*VDSIG*NDV*ETA)+KSS(K)*NDV*(V1_V22-2.*ETA*
     &         V1_V22-1.5*NDV*NDV+3.*ETA*NDV*NDV-MOR+2.*ETA*MOR)+
     &         DM/M*KSSIG(K)*NDV*(-ETA*V1_V22+1.5*ETA*NDV*NDV+
     &         (ETA-1.)*MOR))

        SV2(K) = MOR/RIJ*(NVEC(K)*(VDSIG*ETA*(-2.*V1_V22+6.*ETA*NDV*
     &       NDV+(3.+8.*ETA)*MOR)+MOR*NDSIG*NDV*(2.-22.5*ETA+2.*
     &       ETA*ETA)+
     &       DM/M*VDS*ETA*(-3.*NDV*NDV-4.*MOR)+DM/M*MOR*NDS*NDV*(2.-
     &       0.5*ETA))+KSV(K)*(NDSIG*(0.5*ETA*V1_V22+2.*ETA*ETA*V1_V22-
     &       4.5*ETA*ETA*NDV*NDV+(4.5*ETA-1.+8.*ETA*ETA)*MOR)+VDSIG*NDV*
     &       ETA*(6.*ETA-1.)-3.*DM/M*VDS*NDV*ETA+DM/M*NDS*(-1.5*
     &       ETA*V1_V22+
     &       3.*ETA*NDV*NDV-(1.+3.5*ETA)*MOR))+KSSIG(K)*NDV*(2.*ETA*ETA*
     &       V1_V22-3.*ETA*ETA*NDV*NDV+(-1.+4.*ETA-2.*ETA*ETA)*MOR)+
     &       DM/M*KSS(K)*NDV*(-ETA*V1_V22+1.5*ETA*NDV*NDV+(-1.+ETA)*
     &       MOR))

**********2PN**********
    3 CONTINUE
**belongs to 1.5 SPIN SPIN**
        CALL CROSSP(SS1,S1,SS1)
        CALL CROSSP(SS2,S2,SS2)
**************
      DO 7 K = 1,3
          SU(K) = SU1(K)/CL2 + SU2(K)/CL4 + (SS1(K) + SS2(K))/CL2
          SV(K) = SV1(K)/CL2 + SV2(K)/CL4+M*(SS2(K)/BODY(I2)-SS1(K)/
     &       BODY(I1))/CL2
          KSS(K) = KSS(K) + SU(K)*STEP(I1)
          KSSIG(K) = KSSIG(K) + SV(K)*STEP(I1)
C* calculate actual spins from the total S and Sigma
          SPN(K,I1) = BODY(I1)*(M*KSS(K)-BODY(I2)*KSSIG(K))/M/M/
     &         BODY(I1)/BODY(I1)*CLIGHT

          SPN(K,I2) = BODY(I2)*(M*KSS(K)+BODY(I1)*KSSIG(K))/M/M/
     &         BODY(I2)/BODY(I2)*CLIGHT
C derivatives of tagoshi spin variables
          XAD(K) = 0.5/(M*M*BODY(I1)*BODY(I2))*(-SU(K)*M*DM-SV(K)*
     &          (BODY(I1)*BODY(I1)+BODY(I2)*BODY(I2)))
          XSD(K) = 0.5/(M*M*BODY(I1)*BODY(I2))*(SU(K)*M*M+SV(K)*
     &          (BODY(I1)*BODY(I1)-BODY(I2)*BODY(I2)))
    7 CONTINUE

**DERIVATIVES OF SPIN TERMS**
**Auxiliary terms**
        CALL CROSSP(NDOT,KSS,NDOTCS)
        CALL CROSSP(NVEC,SU,NCSU)
        CALL CROSSP(NDOT,KSSIG,NDOTCSIG)
        CALL CROSSP(NVEC,SV,NCSV)
        CALL CROSSP(AT,KSS,ACS)
        CALL CROSSP(KSV,SU,VCSU)
        CALL CROSSP(AT,KSSIG,ACSIG)
        CALL CROSSP(KSV,SV,VCSV)

        SNVDOT = SU(1)*NCV(1)+SU(2)*NCV(2)+SU(3)*NCV(3)+
     &       KSS(1)*NDOTCV(1)+KSS(2)*NDOTCV(2)+KSS(3)*NDOTCV(3)+
     &       KSS(1)*NCA(1)+KSS(2)*NCA(2)+KSS(3)*NCA(3)

        SIGNVDOT = SV(1)*NCV(1)+SV(2)*NCV(2)+SV(3)*NCV(3)+
     &       KSSIG(1)*NDOTCV(1)+KSSIG(2)*NDOTCV(2)+KSSIG(3)*NDOTCV(3)+
     &       KSSIG(1)*NCA(1)+KSSIG(2)*NCA(2)+KSSIG(3)*NCA(3)

        NSDOT = NDOT(1)*KSS(1)+NDOT(2)*KSS(2)+NDOT(3)*KSS(3)+
     &       NVEC(1)*SU(1)+NVEC(2)*SU(2)+NVEC(3)*SU(3)
        NSIGDOT = NDOT(1)*KSSIG(1)+NDOT(2)*KSSIG(2)+NDOT(3)*KSSIG(3)+
     &       NVEC(1)*SV(1)+NVEC(2)*SV(2)+NVEC(3)*SV(3)
        VSDOT = AT(1)*KSS(1)+AT(2)*KSS(2)+AT(3)*KSS(3)+
     &       KSV(1)*SU(1)+KSV(2)*SU(2)+KSV(3)*SU(3)
        VSIGDOT = AT(1)*KSSIG(1)+AT(2)*KSSIG(2)+AT(3)*KSSIG(3)+
     &       KSV(1)*SV(1)+KSV(2)*SV(2)+KSV(3)*SV(3)

        NXSDOT = NDOT(1)*XS(1)+NDOT(2)*XS(2)+NDOT(3)*XS(3)+
     &       NVEC(1)*XSD(1)+NVEC(2)*XSD(2)+NVEC(3)*XSD(3)
        NXADOT = NDOT(1)*XA(1)+NDOT(2)*XA(2)+NDOT(3)*XA(3)+
     &       NVEC(1)*XAD(1)+NVEC(2)*XAD(2)+NVEC(3)*XAD(3)
******end of auxiliary terms
********1.5 PN********
      DO 6 K = 1,3
        C1_5D(K) = -3.*RP/RIJ*C1_5(K)+(NDOT(K)*(12.*SDNCV+6.*DM/M*
     &       SIGDNCV)+NVEC(K)*(12.*SNVDOT+6.*DM/M*SIGNVDOT)+9.*NVDOT*
     &       NCS(K)+9.*NDV*(NDOTCS(K)+NCSU(K))+3.*DM/M*(NVDOT*NCSIG(K)+
     &       NDV*(NDOTCSIG(K)+NCSV(K)))-7.*(ACS(K)+VCSU(K))-3.*DM/M*
     &       (ACSIG(K)+VCSV(K)))/(RIJ2*RIJ)
*******1.5 PN*********
******2PN SS**********
        C2D(K) = -4.*RP/RIJ*C2(K)-MOR*MOR*MOR*3.*ETA/RIJ*(NDOT(K)*
     &       (XS2-XA2-5.*NXS*NXS+5.*NXA*NXA)+NVEC(K)*(2.*(XS(1)*XSD(1)+
     &       XS(2)*XSD(2)+XS(3)*XSD(3)-XA(1)*XAD(1)-XA(2)*XAD(2)-
     &       XA(3)*XAD(3))-10.*NXS*NXSDOT+10.*NXA*NXADOT)+2.*(XSD(K)*
     &       NXS+XS(K)*NXSDOT-XAD(K)*NXA-XA(K)*NXADOT))
******2PN SS**********

*******2.5 PN*********
        C2_5D(K) = -3.*RP/RIJ*C2_5(K)+(NDOT(K)*(SDNCV*(-30.*ETA*
     &       NDV*NDV+24.*ETA*V1_V22-MOR*(38.+25.*ETA))+DM/M*SIGDNCV*
     &       (-15.*ETA*NDV*NDV+12.*ETA*V1_V22-MOR*(18.+14.5*ETA)))+
     &       NVEC(K)*(SNVDOT*(-30.*ETA*NDV*NDV+24.*ETA*V1_V22-MOR*
     &       (38.+25.*ETA))+SDNCV*(-60.*ETA*NDV*NVDOT+48.*ETA*VA+
     &       MOR*RP/RIJ*(38.+25.*ETA))+DM/M*SIGNVDOT*(-15.*ETA*NDV*
     &       NDV+12.*ETA*V1_V22-MOR*(18.+14.5*ETA))+DM/M*SIGDNCV*
     &       (-30.*ETA*NDV*NVDOT+24.*ETA*VA+MOR*RP/RIJ*(18.+14.5*ETA)))+
     &       (NVDOT*KSV(K)+NDV*AT(K))*(SDNCV*(-9.+9.*ETA)+DM/M*SIGDNCV*
     &       (-3.+6.*ETA))+NDV*KSV(K)*(SNVDOT*(-9.+9.*ETA)+DM/M*
     &       SIGNVDOT*(-3.+6.*ETA))+(NDOTCV(K)+NCA(K))*(NDV*VDS*(-3.+
     &       3.*ETA)-8.*MOR*ETA*NDS-DM/M*(4.*MOR*ETA*NDSIG+3.*NDV*VDSIG)
     &       )+NCV(K)*((NVDOT*VDS+NDV*VSDOT)*(-3.+3.*ETA)-8.*ETA*MOR*
     &       (NSDOT-RP/RIJ*NDS)-DM/M*(4.*ETA*MOR*(NSIGDOT-RP/RIJ*NDSIG)+
     &       3.*(NVDOT*VDSIG+NDV*VSIGDOT)))+(NVDOT*NCS(K)+NDV*
     &       (NDOTCS(K)+NCSU(K)))*(-22.5*ETA*NDV*NDV+21.*ETA*V1_V22-
     &       MOR*(25.+15.*ETA))+NDV*NCS(K)*(-45.*ETA*NDV*NVDOT+42.*ETA*
     &       VA+MOR*RP/RIJ*(25.+15.*ETA))+DM/M*(NVDOT*NCSIG(K)+NDV*
     &       (NDOTCSIG(K)+NCSV(K)))*(-15.*ETA*NDV*NDV+12.*ETA*V1_V22-
     &       MOR*(9.+8.5*ETA))+DM/M*NDV*NCSIG(K)*(-30.*ETA*NDV*NVDOT+
     &       24.*ETA*VA+MOR*RP/RIJ*(9.+8.5*ETA))+(ACS(K)+VCSU(K))*
     &       (16.5*ETA*NDV*NDV+MOR*(21.+9.*ETA)-14.*ETA*V1_V22)+
     &       VCS(K)*(33.*ETA*NDV*NVDOT-MOR*RP/RIJ*(21.+9.*ETA)-
     &       28.*ETA*VA)+DM/M*(ACSIG(K)+VCSV(K))*(9.*ETA*NDV*NDV-
     &       7.*ETA*V1_V22+MOR*(9.+4.5*ETA))+DM/M*VCSIG(K)*(18.*
     &       ETA*NDV*NVDOT-14.*ETA*VA-MOR*RP/RIJ*(9.+4.5*ETA)))/
     &       (RIJ2*RIJ)
*******2.5 PN*********
    6 CONTINUE
        DEBUGC2D = SQRT(C2D(1)**2+C2D(2)**2+C2D(3)**2)

* Sum up the Adots
        ADK = A1D/CL2+A2D/CL4+A2_5D/CL5+A3D/CL6+A3_5D/CL7
	BDK = B1D/CL2+B2D/CL4+B2_5D/CL5+B3D/CL6+B3_5D/CL7

        AD(1) = -2.0*MOR*RP*(KSAK*NX+KSBK*KSV(1))/RIJ2 + 
     &       MOR*(ADK*NX+BDK*KSV(1))/RIJ + MOR*(KSAK*(KSV(1)-NX*RP)/RIJ+
     &         KSBK*AT(1))/RIJ + C1_5D(1)/CL2 + C2D(1)/CL4 +
     &       C2_5D(1)/CL4

        AD(2) = -2.0*MOR*RP*(KSAK*NY+KSBK*KSV(2))/RIJ2 + 
     &       MOR*(ADK*NY+BDK*KSV(2))/RIJ + MOR*(KSAK*(KSV(2)-NY*RP)/RIJ+
     &         KSBK*AT(2))/RIJ + C1_5D(2)/CL2 + C2D(2)/CL4 +
     &       C2_5D(2)/CL4

        AD(3) = -2.0*MOR*RP*(KSAK*NZ+KSBK*KSV(3))/RIJ2 + 
     &       MOR*(ADK*NZ+BDK*KSV(3))/RIJ + MOR*(KSAK*(KSV(3)-NZ*RP)/RIJ+
     &         KSBK*AT(3))/RIJ + C1_5D(3)/CL2 + C2D(3)/CL4 +
     &       C2_5D(3)/CL4
***********************	  
        DO 346 K = 1,3
           FRELD(K) = AD(K)
 346    CONTINUE

* ------------------------------------------------------------
* Caluculations of the relativistic energy (used in adjust.F).
* ------------------------------------------------------------
*
* This calculation approximates the relativistic energy associated
* with the binary using a semi-classical radiative approach:
* E = int(F dot v)dt ~ F dot v *timestep
*                      + (dF/dt dot v + F dot dv/dt)*timestep^2
* and should provide a self-consistent relativistic energy for
* high values of c
*
* First order correction
* ----------------------
	 A2_5 = 1.6*ETA*MOR*RP*(17.0*MOR/3.0+3.0*V1_V22)/CL5
	 B2_5 = -1.6*ETA*MOR*(3.0*MOR+V1_V22)/CL5
      DO 173 MM = 1,3
         A25POWER(MM) = MOR*(A2_5*NVEC(MM) + B2_5*KSV(MM))/RIJ
C         AP(MM) = BODY(I2)*MOR*(A2_5*NVEC(MM) + B2_5*KSV(MM))/RIJ/M/CL5
         AP(MM) = BODY(I2)*FRELP(MM)/M
         AS(MM) = -BODY(I1)*FRELP(MM)/M
C         AS(MM) = -BODY(I1)/BODY(I2)*AP(MM)
 173  CONTINUE
*
* Calculte the power radiated at each order.
      POWER1P = (AP(1)*VI(1)+AP(2)*VI(2)+
     & AP(3)*VI(3))*BODY(I1)
      POWER1S = (AS(1)*VI(4)+AS(2)*VI(5)+
     & AS(3)*VI(6))*BODY(I2)
      DIFFS(1) = POWER1P
      DIFFS(2) = POWER1S
      DIFFS(3) = (POWER1P+POWER1S)
*
* RELTSTEP is defined in ksint.F
C      RELENERGY(IPAIR) = RELENERGY(IPAIR)
C     &     + (POWER1P + POWER1S)*STEP(I1)
C      PREVRELPOTENERGY = PREVRELPOTENERGY + 
C     &   (POWER1P + POWER1S)*STEP(I1)
C      write (*,*) "RELTSTEP AND STEP",RELTSTEP(IPAIR),STEP(IPAIR)
*
* Second order correction
* -----------------------
        A2_5D = -(1.6*ETA*MOR*RP*RP*(17.0/3.0*MOR+3.0*V1_V22)/RIJ+
     &       1.6*ETA*MOR*RPP*(17.0/3.0*MOR+3.0*V1_V22)+
     &       1.6*ETA*MOR*RP*(-17.0*MOR*RP/3.0/RIJ+6.0*VA))/CL5
	 
	B2_5D = (1.6*ETA*MOR*RP*(3.0*MOR+V1_V22)/RIJ - 
     &       1.6*ETA*MOR*(-3.0*MOR*RP/RIJ+2.0*VA))/CL5


      DO 174 NN = 1,3
        A25DPOWER(NN) = -2.0*MOR*RP*(A2_5*NX+B2_5*KSV(1))/RIJ2 + 
     &       MOR*(A2_5D*NX+B2_5D*KSV(1))/RIJ + 
     &        MOR*(A2_5*(KSV(1)-NX*RP)/RIJ+
     &         B2_5D*AT(1))/RIJ
        AT1(NN) = BODY(I2)*AT(NN)/M
        AT2(NN) = -BODY(I1)*AT(NN)/M
         ADP(NN) = BODY(I2)*FRELD(NN)/M
         ADS(NN) = -BODY(I1)*FRELD(NN)/M
 174  CONTINUE
*
      POWER2P = 0.5*(AP(1)*AT1(1)+AP(2)*AT1(2)+AP(3)*AT1(3)
     &     + ADP(1)*VI(1)+ADP(2)*VI(2)+ADP(3)*VI(3))*BODY(I1)
      POWER2S = 0.5*(AS(1)*AT2(1)+AS(2)*AT2(2)+AS(3)*AT2(3)
     &     + ADS(1)*VI(4)+ADS(2)*VI(5)+ADS(3)*VI(6))*BODY(I2)
      DIFFS(4) = POWER2P + POWER2S
C      RELENERGY(IPAIR) = RELENERGY(IPAIR)+
C     &     (POWER2P+POWER2S)*STEP(I1)*STEP(I1)

*     =================================================================
*     Calculate additional influence of perturbers on energy correction
*     =================================================================
*
*     Sum over perturberlist
*     Acc of pair particles given as AP(3) and AS(3) above
      NNB2 = LIST(1,I1) + 1
      DO 101 J = 2, NNB2
         K = LIST(J,I1)
         A1 = X(1,K) - XI(1)
         A2 = X(2,K) - XI(2)
         A3 = X(3,K) - XI(3)
         RIJ2 = A1*A1 + A2*A2 + A3*A3
         RIJ = SQRT(RIJ2)
         RELPOTENERGY(IPAIR) = RELPOTENERGY(IPAIR)+
     &    BODY(I1)*BODY(K)*(AP(1)*A1+AP(2)*A2+AP(3)*A3)*
     &    STEP(I1)*STEP(I1)/(2.*RIJ2*RIJ)
*     And for 2nd body in pair
         A1 = X(1,K) - XI(4)
         A2 = X(2,K) - XI(5)
         A3 = X(3,K) - XI(6)
         RIJ2 = A1*A1 + A2*A2 + A3*A3
         RIJ = SQRT(RIJ2)
         RELPOTENERGY(IPAIR) = RELPOTENERGY(IPAIR)+
     &    BODY(I2)*BODY(K)*(AS(1)*A1+AS(2)*A2+AS(3)*A3)*
     &    STEP(I1)*STEP(I1)/(2.*RIJ2*RIJ)
 101  CONTINUE
*     ===============================================================
*
*
* RELTSTEP is defined in ksint.F
C      WRITE (*,*) 'debugg', (FRELP(1)*KSV(1)+FRELP(2)*KSV(2)+
C     &     FRELP(3)*KSV(3))*RELTSTEP(IPAIR)*M*ETA
C      RELENERGY(IPAIR) = RELENERGY(IPAIR)
C     &    + (FRELP(1)*KSV(1)+FRELP(2)*KSV(2)+
C     &     FRELP(3)*KSV(3))*STEP(I1)*M*ETA

C     &     + (POWER2P + POWER2S)*RELTSTEP(IPAIR)*RELTSTEP(IPAIR)

C      RELSUM = RELSUM + (POWER1P+POWER1S)*RELTSTEP(IPAIR)
C      RELSUM = RELSUM + (A25POWER(1)*KSV(1)+A25POWER(2)*KSV(2)+
C     &     A25POWER(3)*KSV(3))*RELTSTEP(IPAIR)*M
C      RELSUM = RELSUM + (FRELP(1)*KSV(1)+FRELP(2)*KSV(2)+
C     &     FRELP(3)*KSV(3))*RELTSTEP(I1)*M*ETA
C      RELENERGY(IPAIR) = RELENERGY(IPAIR) + 
C     & 0.5*(FRELD(1)*KSV(1)+FRELD(2)*KSV(2)+
C     & FRELD(3)*KSV(3)+
C     & FRELP(1)*AT(1)+FRELP(2)*AT(2)+FRELP(3)*AT(3))*
C     & STEP(I1)*STEP(I1)*M*ETA
*      WRITE (*,*) 'debug ', POWER1P, POWER1S, POWER2P, POWER2S
*      write (*,*) "2nd order ",RELENERGY2_5(IPAIR)
*
*      CALL RENG (IPAIR,I1,KSX,KSV)
*
* -------------------------------------------
* End calculation of the relativistic energy.
* -------------------------------------------
*

      END IF

C      WRITE (*,*) 'Relativistic'
*checking derivative and other tests***
*        DIFFS(1) = (C2_5(1)-DIFFS(2))/STEP(1)
*        DIFFS(2) = C2_5(1)
** CHECK DERIVATIVES **
*        DIFFS(1) = (DIFFS(1)-C2D(3))/C2D(3)
C         DIFFS(2) = SQRT(L(1)*L(1)+L(2)*L(2)+L(3)*L(3))*CL/(BODY(1)+
C     &   BODY(2))**2
C         DIFFS(3) = SQRT(SPN(1,1)**2+SPN(2,1)**2+SPN(3,1)**2)
C         DIFFS(4) = SQRT(SPN(1,2)**2+SPN(2,2)**2+SPN(3,2)**2)
C         DIFFS(1) = SQRT((L(1)+SPN(1,1)*BODY(1)*BODY(1)/CL+SPN(1,2)*
C     &   BODY(2)*BODY(2)/CL)**2+(L(2)+SPN(2,1)*BODY(1)*BODY(1)/CL+
C     &   SPN(2,2)*BODY(2)*BODY(2)/CL)**2+(L(3)+SPN(3,1)*BODY(1)*
C     &   BODY(1)/CL+SPN(3,2)*BODY(2)*BODY(2)/CL)**2)*CL/(BODY(1)+
C     &   BODY(2))**2
C        IF ((A3D.GT.1.0D30).OR.(A3D.LT.-1.0D30)) THEN
*          DIFFS(1) = (DIFFS(1)-C2_5D(1))/C2_5D(1)
C        ELSE
C          DIFFS(1) = 0.
C        END IF
*/checking derivative**
C        DIFFS(1) = SQRT((SU(1)*SU(1)+SU(2)*SU(2)+SU(3)*SU(3))/
C     &  (KSS(1)*KSS(1)+KSS(2)*KSS(2)+KSS(3)*KSS(3)))*STEP(1)

** PN strength survey
C       ANEWTON = MOR/RIJ
C       DIFFS(1) = SQRT((A1*NX+B1*KSV(1))**2+(A1*NY+B1*KSV(2))**2+
C     & (A1*NZ+B1*KSV(3))**2)/CL2
C       DIFFS(2) = SQRT((A2*NX+B2*KSV(1))**2+(A2*NY+B2*KSV(2))**2+
C     & (A2*NZ+B2*KSV(3))**2)/CL4
C      DIFFS(3) = SQRT((A2_5*NX+B2_5*KSV(1))**2+(A2_5*NX+B2_5*KSV(1))**2+
C     & (A2_5*NZ+B2_5*KSV(3))**2)/CL5
C       DIFFS(4) = SQRT((A3*NX+B3*KSV(1))**2+(A3*NY+B3*KSV(2))**2+
C     & (A3*NZ+B3*KSV(3))**2)/CL6
C      DIFFS(5) = SQRT((A3_5*NX+B3_5*KSV(1))**2+(A3_5*NY+B3_5*KSV(2))**2+
C     & (A3_5*NZ+B3_5*KSV(3))**2)/CL7
C       DIFFS(6) = SQRT(C1_5(1)**2+C1_5(2)**2+C1_5(3)**2)/ANEWTON/CL2
C       DIFFS(7) = SQRT(C2(1)**2+C2(2)**2+C2(3)**2)/ANEWTON/CL4
C       DIFFS(8) = SQRT(C2_5(1)**2+C2_5(2)**2+C2_5(3)**2)/ANEWTON/CL4
C       DIFFS(9) = SQRT(SQM(1)**2+SQM(2)**2+SQM(3)**2)/ANEWTON/CL4
C       DIFFS(5) = (SPN(1,1)*LU(1)+SPN(2,1)*LU(2)+SPN(3,1)*LU(3))/
C     & SQRT(SPN(1,1)**2+SPN(2,1)**2+SPN(3,1)**2)
C       DIFFS(6) = (SPN(1,2)*LU(1)+SPN(2,2)*LU(2)+SPN(3,2)*LU(3))/
C     & SQRT(SPN(1,2)**2+SPN(2,2)**2+SPN(3,2)**2)
C       DIFFS(7) = (SPN(1,1)*SPN(1,2)+SPN(2,1)*SPN(2,2)+
C     &  SPN(3,1)*SPN(3,2))/DIFFS(3)/DIFFS(4)
C       DIFFS(5) = 180.*ACOS(DIFFS(5))/3.14
C       DIFFS(6) = 180.*ACOS(DIFFS(6))/3.14
C       DIFFS(7) = 180.*ACOS(DIFFS(7))/3.14

**/PN bla

** final spin determination after rezzolla
       IF (APPLY.EQ.3) THEN
         ABSA1 = SQRT(SPN(1,1)**2+SPN(2,1)**2+SPN(3,1)**2)
         ABSA2 = SQRT(SPN(1,2)**2+SPN(2,2)**2+SPN(3,2)**2)
         IF (ABSA1.LE.1.0E-10) THEN
           COSA = 0.
           COSB = 0.
         ELSE
         COSA = (SPN(1,1)*SPN(1,2)+SPN(2,1)*SPN(2,2)+
     &     SPN(3,1)*SPN(3,2))/ABSA1/ABSA2
         COSB = (SPN(1,1)*L(1)+SPN(2,1)*L(1)+SPN(3,1)*L(3))/
     &     ABSA1/LABS
         END IF
         IF (ABSA2.LE.1.0E-10) THEN
           COSG = 0.
         ELSE
         COSG = (SPN(1,2)*L(1)+SPN(2,2)*L(1)+SPN(3,2)*L(3))/
     &     ABSA2/LABS
         END IF
         ABSL = -0.129/(1+RATIO*RATIO)**2*(ABSA1**2+ABSA2**2*RATIO**4+
     &     2.*ABSA1*ABSA2*RATIO**2*COSA)+(-3.84*ETA-2.686+2.)/
     &     (1+RATIO**2)*(ABSA1*COSB+ABSA2*RATIO**2*COSG)+
     &     3.4641-3.454*ETA+2.353*ETA*ETA
         AFINABS = 1/((1+RATIO)**2)*SQRT(ABSA1**2+ABSA2**2*RATIO**4+2.*
     &     ABSA2*ABSA1*RATIO**2*COSA+2.*(ABSA1*COSB+ABSA2*RATIO**2*
     &     COSG)*ABSL*RATIO+ABSL*ABSL*RATIO*RATIO)
         JVEC(1) = KSS(1) + L(1)
         JVEC(2) = KSS(2) + L(2)
         JVEC(3) = KSS(3) + L(3)
         JABS = SQRT(JVEC(1)**2+JVEC(2)**2+JVEC(3)**2)
         AFIN(1) = AFINABS*JVEC(1)/JABS
         AFIN(2) = AFINABS*JVEC(2)/JABS
         AFIN(3) = AFINABS*JVEC(3)/JABS
         WRITE (*,*) 'FINALSPIN', AFIN(1), AFIN(2), AFIN(3)

      END IF
*     If the actual total spin is lower than the prediction,
*     it is definitely time to merge!
C         IF (DIFFS(1).LT.AFIN) THEN
C            WRITE (*,*) 'SPIN CRITERION FULFILLED, Merge in next cycle!'
*     Set some flag here to merge in next iteration...
C         END IF
            
C         DIFFS(8) = AFIN

**plot complete L,S1,S2 information
C         DIFFS(1) = L(1)+SPN(1,1)*BODY(1)**2/CL+SPN(1,2)*BODY(2)**2/CL
C         DIFFS(2) = L(2)+SPN(2,1)*BODY(1)**2/CL+SPN(2,2)*BODY(2)**2/CL
C         DIFFS(3) = L(3)+SPN(3,1)*BODY(1)**2/CL+SPN(3,2)*BODY(2)**2/CL
C         LABS = DIFFS(1)*DIFFS(1)+DIFFS(2)*DIFFS(2)+DIFFS(3)*DIFFS(3)
C         DIFFS(1) = DIFFS(1)/SQRT(LABS)
C         DIFFS(2) = DIFFS(2)/SQRT(LABS)
C         DIFFS(3) = DIFFS(3)/SQRT(LABS)
C         DIFFS(4) = SPN(1,1)*BODY(1)**2/CL
C         DIFFS(5) = SPN(2,1)*BODY(1)**2/CL
C         DIFFS(6) = SPN(3,1)*BODY(1)**2/CL
C         DIFFS(7) = SPN(1,2)*BODY(2)**2/CL
C         DIFFS(8) = SPN(2,2)*BODY(2)**2/CL
C         DIFFS(9) = SPN(3,2)*BODY(2)**2/CL
C         DIFFS(4) = LU(1)
C         DIFFS(5) = LU(2)
C         DIFFS(6) = LU(3)
**********************
	 V12 = VI(1)*VI(1)+VI(2)*VI(2)+VI(3)*VI(3)
         V22 = VI(4)*VI(4)+VI(5)*VI(5)+VI(6)*VI(6)
	  
         VWHOLE = MAX(SQRT(V12),SQRT(V22))

	 EKIN = BODY(I1)*V12/2.0 + BODY(I2)*V22/2.0
         EPOTEN = -BODY(I1)*BODY(I2)/RIJ
*** add spin energy
C         ESPIN = BODY(I1)*BODY(I1)*(SQRT(SPN(1,I1)*SPN(1,I1)+
C     &       SPN(2,I1)*SPN(2,I1)+SPN(3,I1)*SPN(3,I1)))+BODY(I2)*
C     &       BODY(I2)*(SQRT(SPN(1,I2)*SPN(1,I2)+
C     &       SPN(2,I2)*SPN(2,I2)+SPN(3,I2)*SPN(3,I2)))

	 ENEWT = EKIN + EPOTEN
C         WRITE (*,*) 'ENERGYCHECK',TIME+TOFF, RELENERGY(IPAIR), ENEWT
C         IF (I1.EQ.1) WRITE (*,*) TIME+TOFF, DIFFS(1), DIFFS(6),
C     &   DIFFS(7), DIFFS(8)
         KSA1(1) =  BODY(I2)*A(1)/M
         KSA1(2) =  BODY(I2)*A(2)/M
         KSA1(3) =  BODY(I2)*A(3)/M
	
         KSA2(1) =  -BODY(I1)*A(1)/M
	 KSA2(2) =  -BODY(I1)*A(2)/M
         KSA2(3) =  -BODY(I1)*A(3)/M
	          
	 POWER1 = (KSA1(1)*VI(1)+KSA1(2)*VI(2)+KSA1(3)*VI(3))
     &                 *BODY(I1)
         POWER2 = (KSA2(1)*VI(4)+KSA2(2)*VI(5)+KSA2(3)*VI(6))
     &	          *BODY(I2)

          RETURN	  
	  END
