      SUBROUTINE KSREL(IPAIR,I1,XI,VI,FRELP1,FRELP2,FRELP2_5,FRELD1,
     &                 FRELD2,FRELD2_5,VWHOLE)
*
* Relativistic perturbation on KS pair.
* -------------------------------------
*
* --------------------------------------------------------
* Created by Jonathan M.B. Downing (02.2007) based on work
* by Gabor Kupi and Pau Amaro-Seoane in NBODY4.
* --------------------------------------------------------
*
      INCLUDE 'common6.h'
* --------------------------------------------------------
      COMMON/CHAINC/ XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &               LISTC(LMAX)
*
* Names and Counters
      INTEGER NAME1,NAME2,I,J,K,L,II,JJ,KK,LL
*
* Input
      REAL*8 XI(6),VI(6)
*
* Constants
      REAL*8 PI2,CL,CL2,CL4,CL5,CL6,CL7
*
* Calculated from input
      REAL*8 RIJ,RIJ2,ETA,M_BOTH,RP,MOR,RPP
      REAL*8 KSV(3),KSX(3),KSV1(3),KSV2(3),N12(3)
      REAL*8 KSADOT1(3),DSADOT2(3)
*
* Parameters for acceleration
      REAL*8 A1,A2,A2_5,B1,B2,B2_5
      REAL*8 KSAK1,KSBK1,KSAK2,KSBK2,KSAK2_5,KSBK2_5,KSAK,KSBK
      REAL*8 AC1(3),AC2(3),AC2_5(3),A(3),AT(3)
      REAL*8 FRELP1(3),FRELP2(3),FRELP2_5(3)
*
* Paramaters for derivative of acceleration
      REAL*8 A1D,A2D,A2_5D,B1D,B2D,B2_5D,VA
      REAL*8 ADK1,BDK1,ADK2,BDK2,ADK2_5,BDK2_5
      REAL*8 ACD1(3),ACD2(3),ACD2_5(3),AD(3),AD1(3),AD2(3)
      REAL*8 FRELD1(3),FRELD2(3),FRELD2_5(3)
*
* Diagnostics
      REAL*8 FRP1_2,FRP2_2,FRP2_5_2,FRP1,FRP2,FRP2_5
      REAL*8 AP1(3),AS1(3),AP2(3),AS2(3),AP2_5(3),AS2_5(3)
      REAL*8 ADP1(3),ADS1(3),ADP2(3),ADS2(3),ADP2_5(3),ADS2_5(3)
      REAL*8 POWER1P(3),POWER1S(3),POWER2P(3),POWER2S(3)
      REAL*8 KSIMPMOM2,KSECC,KSH
      SAVE IIIK
* ---------------------------------------------------------------
*
* Write the headers for the relativistic datafiles (once per simulation).
      IF (FIRSTREL.EQ.0) THEN
*
         WRITE (35,7641) 'TIME','JPAIR','I','J','MASS(I)','MASS(J)',
     &                   'KSTAR(I)','KSTAR(J)','EREL','SEMI','ZN',
     &                   'ORBSEP','RI','ECC','GAMMA','GAMMA_1PN',
     &                   'GAMMA_2PN','GAMMA_2.5PN','X1','Y1','Z1',
     &        'X2','Y2','Z2','Vx1','Vy1','Vz1','Vx2','Vy2','Vz2'
 7641    FORMAT (A12,3A6,2A14,2A10,22A14)
*
         FIRSTREL = 1
*
      END IF
*
* Report a new relativistic encounter in the output file.
      I2 = I1+1
      IF (KSTHEFIRST.EQ.0) THEN
*
         WRITE (*,*) ' '
         WRITE (*,7643) TIME+TOFF,VOVERC,NAME(I1),NAME(I2)
 7643    FORMAT ('Relativistic Event -->  Time = ',F20.8,
     &           ' VOVERC = ',F20.8,' I1 = ',I8,' I2 = ',I8)
*
         KSTHEFIRST = 1
*
      END IF
*
* Uncomment this GOTO to get additional (e,a,etc.) output about the
* KS perturbers of particles 1 and/or 2 without PN terms.
*
*      GOTO 1362
*
* Set a flag showing that the binary is relativistic.
      RFLAG(IPAIR) = 1
*
* Initialization of constants and KS quantities for PN calculation.
      PI2 = 9.869604401089359
*
* Speed of light, CL, and its powers
      CL = CLIGHT
      CL2 = CL*CL
      CL4 = CL2*CL2
      CL5 = CL4*CL
      CL6 = CL5*CL
      CL7 = CL6*CL
*
* Temporary variables
      FRP1_2 = 0.0
      FRP2_2 = 0.0
      FRP2_5_2 = 0.0
*
* Center of mass properties -- mass
      M_BOTH = BODY(I1)+BODY(I2)
      MU = BODY(I1)*BODY(I2)/M_BOTH
      RCSCH = 2.0*M/CL2
      ETA = BODY(I1)*BODY(I2)/(M_BOTH*M_BOTH)
*
* Initialization of relativistic paramaters and input conversion.
      DO 111 L = 1,3
*
         FRELP1(L) = 0.0
         FRELP1(L+3) = 0.0
         FRELP2(L) = 0.0
         FRELP2(L+3) = 0.0
         FRELP2_5(L) = 0.0
         FRELP2_5(L+3) = 0.0
*
         FRELD1(L) = 0.0
         FRELD1(L+3) = 0.0
         FRELD2(L) = 0.0
         FRELD2(L+3) = 0.0
         FRELD2_5(L) = 0.0
         FRELD2_5(L+3) = 0.0
*
         KSX(L) = XI(L)-XI(L+3)
*
         KSV(L) = VI(L)-VI(L+3)
*
 111  CONTINUE
*
* Center of mass properties -- radius and velocity
      RIJ2 = KSX(1)*KSX(1)+KSX(2)*KSX(2)+KSX(3)*KSX(3)
      RIJ = SQRT(RIJ2)
      MOR = M_BOTH/RIJ
      V1_V22 = KSV(1)*KSV(1)+KSV(2)*KSV(2)+KSV(3)*KSV(3)
      VWHOLE = SQRT(V1_V22)
*
* Unit vectors.
      N12(1) = KSX(1)/RIJ
      N12(2) = KSX(2)/RIJ
      N12(3) = KSX(3)/RIJ
*
      RP = (KSX(1)*KSV(1)+KSX(2)*KSV(2)+KSX(3)*KSV(3))/RIJ
*
* ------------------------------------
* Calculation of the PN accelerations.
* ------------------------------------
*
* --------------
* 1PN correction
*
      IF (KZ(41).EQ.0.OR.KZ(41).EQ.1.OR.KZ(41).EQ.4) THEN !See define.f
*
* Calculation of the 1PN coefficients.
         A1 = 2.0*(2.0+ETA)*MOR-(1.0+3.0*ETA)*V1_V22+1.5*ETA*RP*RP
         B1 = 2.0*(2.0-ETA)*RP
*
* Calculation of the 1PN acceleration.
         KSAK1 = A1/CL2
         KSBK1 = B1/CL2
*
         DO 121 I = 1,3
            AC1(I) = MOR*(KSAK1*N12(I) + KSBK1*KSV(I))/RIJ
            FRELP1(I) = AC1(I)
            FRP1_2 = FRP1_2 + FRELP1(I)*FRELP1(I)
 121     CONTINUE
*
         FRP1 = SQRT(FRP1_2)
         GAMMAREL1(IPAIR) = FRP1*RIJ*RIJ/M_BOTH
      END IF
*
* End 1PN correction
* ------------------
*
* --------------
* 2PN correction
*
      IF (KZ(41).EQ.0.OR.KZ(41).EQ.2.OR.KZ(41).EQ.4) THEN    !See define.f
*
* Calculation of the 2PN coefficients
         A2 = -0.75*(12.0+29.0*ETA)*MOR*MOR 
     &        -ETA*(3.0-4.0*ETA)*V1_V22*V1_V22 
     &        -1.875*ETA*(1.0-3.0*ETA)*RP*RP*RP*RP
     &        +0.5*ETA*(13.0-4.0*ETA)*MOR*V1_V22
     &        +(2.0+25.0*ETA+2.0*ETA*ETA)*MOR*RP*RP
     &        +1.5*ETA*(3.0-4.0*ETA)*V1_V22*RP*RP
*
         B2 = -0.5*RP*((4.0+41.0*ETA+8.0*ETA*ETA)*MOR
     &        -ETA*(15.0+4.0*ETA)*V1_V22+3.0*ETA*(3.0+2.0*ETA)*RP*RP)
*
* Calculation of the 2PN acceleration
         KSAK2 = A2/CL4
         KSBK2 = B2/CL4
*
         DO 131 J = 1,3
            AC2(J) = MOR*(KSAK2*N12(J) + KSBK2*KSV(J))/RIJ
            FRELP2(J) = AC2(J)
            FRP2_2 = FRP2_2 + FRELP2(J)*FRELP2(J)
 131     CONTINUE
*
         FRP2 = SQRT(FRP2_2)
         GAMMAREL2(IPAIR) = FRP2*RIJ*RIJ/M_BOTH
      END IF
*
* End 2PN correction
* ------------------
*
* ----------------
* 2.5PN correction
*
      IF (KZ(41).EQ.0.OR.3) THEN !See define.f
*
* Calculation of the 2.5PN coefficients
         A2_5 = 1.6*ETA*MOR*RP*(17.0*MOR/3.0+3.0*V1_V22)
         B2_5 = -1.6*ETA*MOR*(3.0*MOR+V1_V22)
*
* Calculation of the 2.5PN acceleration
         KSAK2_5 = A2_5/CL5
         KABK2_5 = B2_5/CL5
*
         DO 141 K = 1,3
            AC2_5(K) = MOR*(KSAK2_5*N12(K) + KSBK2_5*KSV(K))/RIJ
            FRELP2_5(K) = AC2_5(K)
            FRP2_5_2 = FRP2_5_2 + FRELP2_5(K)*FRELP2_5(K)
 141     CONTINUE
*
         FRP2_5 = SQRT(FRP2_5_2)
         GAMMAREL2_5(IPAIR) = FRP2_5_2*RIJ*RIJ/M_BOTH
      END IF
*
* End 2.5PN correction
* --------------------
*
*       Original by Pau has 3PN and 3.5PN expressions here but they are
*       commented out.  The input for this code only goes to 2.5PN so
*       I have not included the 3PN and 3.5PN terms.  They may be included
*       in a more advanced version of this code at a later date. JMBD
*
* ---------------------------------------
* End calculation of the PN acceleration.
* ---------------------------------------
*
* Parameters for calculation of PN acceleration derivatives
      KSAK = KSAK1 + KSAK2 + KSAK2_5
      KSBK = KSBK1 + KSBK2 + KSBK2_5
*
* PN acceleration
      A(1) = MOR*(KSAK*N12(1) + KSBK*KSV(1))/RIJ
      A(2) = MOR*(KSAK*N12(2) + KSBK*KSV(2))/RIJ
      A(3) = MOR*(KSAK*N12(3) + KSBK*KSV(3))/RIJ
*
* PN acceleration + newtonian acceleration for derivatives
      AT(1) = -MOR/RIJ + A(1)
      AT(2) = -MOR/RIJ + A(2)
      AT(3) = -MOR/RIJ + A(3)
*
* More PN derivative parameters
      RPP = V1_V22/RIJ + AT(1)*N12(1)+AT(2)*N12(2)+AT(3)*N12(3)
     &      - RP*RP/RIJ
      VA = AT(1)*KSV(1)+AT(2)*KSV(2)+AT(3)*KSV(3)
*
* ----------------------------------------------------------
* Calculation of the 1st derivative of the PN accelerations.
* ----------------------------------------------------------
*
* ----------------------------
* 1PN acceleration derivative.
*
      IF (KZ(41).EQ.0.OR.KZ(41).EQ.1.OR.KZ(41).EQ.4) THEN  !See define.f
*
* Calculation of the 1PN coefficients.
         A1D = -2.0*(2.0+ETA)*MOR*RP/RIJ - 2.0*(1.0+3.0*ETA)*VA
     &         +3.0*ETA*RP*RPP
         B1D = 2.0*(2.0-ETA)*RPP
*
* Calculation of the 1PN acceleration derivative.
         ADK1 = A1D/CL2
         BDK1 = B1D/CL2
*
         DO 151 II = 1,3
            ACD1(II) = -2.0*MOR*RP*(KSAK1*N12(II)+KSBK1*KSV(II))/RIJ2
     &            +MOR*(ADK1*N12(II)+BDK1*KSV(II))/RIJ
     &            +MOR*(KSAK1*(KSV(II)-N12(II)*RP)/RIJ+KSBK1*AT(II))/RIJ
*
            FRELD1(II) = ACD1(II)
 151     CONTINUE
*
      END IF
*
* End 1PN acceleration derivative.
* --------------------------------
*
* ----------------------------
* 2PN acceleration derivative.
*
      IF (KZ(41).EQ.0.OR.KZ(41).EQ.2.OR.KZ(41).EQ.4) THEN     !See define.f
*
* Calculation of the 2PN coefficients.
         A2D = 1.5*(12.0+29.0*ETA)*MOR*MOR*RP/RIJ
     &         -ETA*(3.0-4.0*ETA)*4.0*V1_V22*VA
     &         -7.5*ETA*(1.0-3.0*ETA)*RPP
     &         -0.5*ETA*(13.0-4.0*ETA)*MOR*RP*V1_V22/RIJ
     &         +ETA*(13.0-4.0*ETA)*MOR*VA
     &         -(2.0+25.0*ETA+2.0*ETA*ETA)*MOR*RP*RP*RP/RIJ
     &         +2.0*(2.0+25.0*ETA+2.0*ETA*ETA)*MOR*RP*RPP
     &         +3.0*ETA*(3.0-4.0*ETA)*VA*RP*RP
     &         +3.0*ETA*(3.0-4.0*ETA)*V1_V22*RP*RPP
*
         B2D = -0.5*RPP*((4.0+41.0*ETA+8.0*ETA*ETA)*MOR
     &         -ETA*(15.0+4.0*ETA)*V1_V22+3.0*ETA*(3.0+2.0*ETA)*RP*RP)
     &         -0.5*RP*(-(4.0+41.0*ETA+8.0*ETA*ETA)*MOR*RP/RIJ
     &         -2.0*ETA*(15.0+4.0*ETA)*VA
     &         +6.0*ETA*(3.0+2.0*ETA)*RP*RPP)
*
* Calculation of the 2PN acceleration derivative.
         ADK2 = A2D/CL4
         BDK2 = B2D/CL4
*
         DO 161 JJ = 1,3
            ACD2(JJ) = -2.0*MOR*RP*(KSAK2*N12(JJ)+KSBK2*KSV(JJ))/RIJ2
     &            +MOR*(ADK2*N12(JJ)+BDK2*KSV(JJ))/RIJ
     &            +MOR*(KSAK2*(KSV(JJ)-N12(JJ)*RP)/RIJ+KSBK2*AT(JJ))/RIJ
*
            FRELD2(JJ) = ACD2(JJ)
 161     CONTINUE
*
      END IF
*
* End 2PN acceleration derivative.
* --------------------------------
*
* ------------------------------
* 2.5PN acceleration derivative.
*
      IF (KZ(41).EQ.0.OR.3) THEN               !See define.f
*
* Calculation of the 2.5PN coefficients
         A2_5D = -1.6*ETA*MOR*RP*RP*(17.0/3.0*MOR+3.0*V1_V22)/RIJ
     &           +1.6*ETA*MOR*RPP*(17.0/3.0*MOR+3.0*V1_V22)
     &           +1.6*ETA*MOR*RP*(-17.0*MOR*RP/3.0/RIJ+6.0*VA)
*
         B2_5D = 1.6*ETA*MOR*RP*(3.0*MOR+V1_V22)/RIJ
     &           -1.6*ETA*MOR*(-3.0*MOR*RP/RIJ+2.0*VA)
*
* Calculation of the 2.5PN acceleration derivative.
         ADK2_5 = A2_5D/CL5
         BDK2_5 = B2_5D/CL5
*
         DO 171 KK = 1,3
         ACD2_5(KK) = -2.0*MOR*RP*(KSAK2_5*N12(KK)+KSBK2_5*KSV(KK))/RIJ2
     &        +MOR*(ADK2_5*N12(KK)+BDK2_5*KSV(KK))/RIJ
     &        +MOR*(KSAK2_5*(KSV(KK)-N12(KK)*RP)/RIJ+KSBK2_5*AT(KK))/RIJ
*
            FRELD2_5(KK) = ACD2_5(KK)
 171     CONTINUE
*
      END IF
*
* End 2.5PN acceleration derivative.
* ----------------------------------
*
*       Original by Pau has 3PN and 3.5PN expressions here but they are
*       commented out.  The input for this code only goes to 2.5PN so
*       I have not included the 3PN and 3.5PN terms.  They may be included
*       in a more advanced version of this code at a later date. JMBD
*
* ---------------------------------------------------------
* End calculation of the 1st derivative of PN acceleration.
* ---------------------------------------------------------
*
* Total derivative of PN acceleration
      DO 172 LL = 1,3
         AD(LL) = ACD1(LL) + ACD2(LL) + ACD2_5(LL)
 172  CONTINUE
*
* Uncomment the CONTINUE to get additional (e,a,etc.) output about
* the KS perturbers of particles 1 and/or 2 without PN terms.
*
* 1362 CONTINUE
*
* ------------------------------------------------------------
* Caluculations of the relativistic energy (used in adjust.F).
* ------------------------------------------------------------
*
* This calculation approximates the relativistic energy associated
* with the binary using a semi-classical radiative approach:
* E = int(F dot v)dt ~ F dot v *timestep
*                      + (dF/dt dot v + F dot dv/dt)*timestep
* and should provide a self-consistent relativistic energy for
* high values of c
*
* First order correction
* ----------------------
      DO 173 MM = 1,3
         AP1(MM) = BODY(I2)*AC1(MM)/M_BOTH
         AS1(MM) = -BODY(I1)*AC1(MM)/M_BOTH
         AP2(MM) = BODY(I2)*AC2(MM)/M_BOTH
         AS2(MM) = -BODY(I1)*AC2(MM)/M_BOTH
         AP2_5(MM) = BODY(I2)*AC2_5(MM)/M_BOTH
         AS2_5(MM) = - BODY(I1)*AC2_5(MM)/M_BOTH
 173  CONTINUE
*
* Calculte the power radiated at each order.
      POWER1P(1) = (AP1(1)*VI(1)+AP1(2)*VI(2)+AP1(3)*VI(3))*BODY(I1)
      POWER1P(2) = (AP2(1)*VI(1)+AP2(2)*VI(2)+AP2(3)*VI(3))*BODY(I1)
      POWER1P(3) = (AP2_5(1)*VI(1)+AP2_5(2)*VI(2)+AP2_5(3)*VI(3))
     &     *BODY(I1)
      POWER1S(1) = (AS1(1)*VI(4)+AS1(2)*VI(5)+AS1(3)*VI(6))*BODY(I2)
      POWER1S(2) = (AS2(1)*VI(4)+AS2(2)*VI(5)+AS2(3)*VI(6))*BODY(I2)
      POWER1S(3) = (AS2_5(1)*VI(4)+AS2_5(2)*VI(5)+AS2_5(3)*VI(6))
     &     *BODY(I2)
*
* RELTSTEP is defined in ksint.F
      RELENERGY1(IPAIR) = RELENERGY1(IPAIR)
     &     + (POWER1P(1) + POWER1S(1))*RELTSTEP(IPAIR)
      RELENERGY2(IPAIR) = RELENERGY2(IPAIR)
     &     + (POWER1P(2) + POWER1S(2))*RELTSTEP(IPAIR)
      RELENERGY2_5(IPAIR) = RELENERGY2_5(IPAIR)
     &     + (POWER1P(3) + POWER1S(3))*RELTSTEP(IPAIR)
*      write (*,*) "1st order ",RELENERGY2_5(IPAIR)
*
* Second order correction
* -----------------------
      DO 174 NN = 1,3
         ADP1(NN) = BODY(I2)*ACD1(NN)/M_BOTH
         ADS1(NN) = -BODY(I1)*ACD1(NN)/M_BOTH
         ADP2(NN) = BODY(I2)*ACD2(NN)/M_BOTH
      ADS2(NN) = -BODY(I1)*ACD2(NN)/M_BOTH
         ADP2_5(NN) = BODY(I2)*ACD2_5(NN)/M_BOTH
         ADS2_5(NN) = - BODY(I1)*ACD2_5(NN)/M_BOTH
 174  CONTINUE
*
      POWER2P(1) = 0.5*(AP1(1)*AP1(1)+AP1(2)*AP1(2)+AP1(3)*AP1(3)
     &     + ADP1(1)*VI(1)+ADP1(2)*VI(2)+ADP1(3)*VI(3))*BODY(I1)
      POWER2P(2) = 0.5*(AP2(1)*AP2(1)+AP2(2)*AP2(2)+AP2(3)*AP2(3)
     &     + ADP2(1)*VI(1)+ADP2(2)*VI(2)+ADP2(3)*VI(3))*BODY(I1)
      POWER2P(3) = 0.5*(AP2_5(1)*AP2_5(1)+AP2_5(2)*AP2_5(2)
     &     + AP2_5(3)*AP2_5(3)+ADP2_5(1)*VI(1)+ADP2_5(2)*VI(2)
     &     + ADP2_5(3)*VI(3))*BODY(I1)
      POWER2S(1) = 0.5*(AS1(1)*AS1(1)+AS1(2)*AS1(2)+AS1(3)*AS1(3)
     &     + ADS1(1)*VI(4)+ADS1(2)*VI(5)+ADS1(3)*VI(6))*BODY(I2)
      POWER2S(2) = 0.5*(AS2(1)*AS2(1)+AS2(2)*AS2(2)+AS2(3)*AS2(3)
     &     + ADS2(1)*VI(4)+ADS2(2)*VI(5)+ADS2(3)*VI(6))*BODY(I2)
      POWER2S(3) = 0.5*(AS2_5(1)*AS2_5(1)+AS2_5(2)*AS2_5(2)
     &     + AS2_5(3)*AS2_5(3)+ADS2_5(1)*VI(4)+ADS2_5(2)*VI(5)
     &     + ADS2_5(3)*VI(6))*BODY(I2)
*
* RELTSTEP is defined in ksint.F
      RELENERGY1(IPAIR) = RELENERGY1(IPAIR)
     &     + (POWER2P(1) + POWER2S(1))*RELTSTEP(IPAIR)*RELTSTEP(IPAIR)
      RELENERGY2(IPAIR) = RELENERGY2(IPAIR)
     &     + (POWER2P(2) + POWER2S(2))*RELTSTEP(IPAIR)*RELTSTEP(IPAIR)
      RELENERGY2_5(IPAIR) = RELENERGY2_5(IPAIR)
     &     + (POWER2P(3) + POWER2S(3))*RELTSTEP(IPAIR)*RELTSTEP(IPAIR)
*      write (*,*) "2nd order ",RELENERGY2_5(IPAIR)
*
*      CALL RENG (IPAIR,I1,KSX,KSV)
*
* -------------------------------------------
* End calculation of the relativistic energy.
* -------------------------------------------
*
*
* Calculation of Newtonian properties of the KS pair.
      V12 = VI(1)*VI(1) + VI(2)*VI(2) + VI(3)*VI(3)
      V22 = VI(4)*VI(4) + VI(5)*VI(5) + VI(6)*VI(6)
*
* Another diagnostic.
      VWHOLE = MAX(SQRT(V12),SQRT(V22))
*
      IF (ILL(IPAIR).EQ.0) THEN
*
         KSV1(1) = BODY(I2)*KSV(1)/M_BOTH
         KSV1(2) = BODY(I2)*KSV(2)/M_BOTH
         KSV1(3) = BODY(I2)*KSV(3)/M_BOTH
*
         KSV2(1) = -BODY(I1)*KSV(1)/M_BOTH
         KSV2(2) = -BODY(I1)*KSV(2)/M_BOTH
         KSV2(3) = -BODY(I1)*KSV(3)/M_BOTH
*
         RELENERGY0(IPAIR) = + RELENERGY2_5(IPAIR)
*
         ILL(IPAIR) = 1
*
      END IF
*
* Not sure what the maximum length for an integer is, hence:
*
      IIIK = IIIK + 1
      IF (IIIK.GE.10000000) IIIK = 0
*
* Set in (IIIK/NUMBER)*NUMBER.EQ.IIIK NUMBER to a lower value to
* get more output - this is because (IIIK/NUMBER) is an integer,
* NUMBER must always be expressed as an interger.
*
      IF ((IIIK/1000000)*1000000.EQ.IIIK) THEN
*
         XCOMPIMP = KSX(2)*KSV(3) - KSX(3)*KSV(2)
         YCOMPIMP = KSX(3)*KSV(1) - KSX(1)*KSV(3)
         ZCOMPIMP = KSX(1)*KSV(2) - KSX(2)*KSV(1)
*
         KSIMPMOM2 = XCOMPIMP*XCOMPIMP + YCOMPIMP*YCOMPIMP
     &               + ZCOMPIMP*ZCOMPIMP
*
         KSH = 0.5*V1_V22 - (BODY(I1)+BODY(I2))/RIJ
*
         IF (1.0+2.0*KSH*KSIMPMOM2/(BODY(I1)+BODY(I2))
     &      /(BODY(I1)+BODY(I2)).GE.0.0) THEN
*
            KSECC = SQRT(1.0+2.0*KSH*KSIMPMOM2/(BODY(I1)+BODY(I2))
     &              /(BODY(I1)+BODY(I2)))
*
         END IF
*
         KSV1(1) = BODY(I2)*KSV(1)/M_BOTH
         KSV1(2) = BODY(I2)*KSV(2)/M_BOTH
         KSV1(3) = BODY(I2)*KSV(3)/M_BOTH
*
         KSV2(1) = -BODY(I1)*KSV(1)/M_BOTH
         KSV2(2) = -BODY(I1)*KSV(2)/M_BOTH
         KSV2(3) = -BODY(I1)*KSV(3)/M_BOTH
*
      END IF
*
      RETURN
*
      END
