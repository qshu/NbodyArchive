      SUBROUTINE KSREL(IPAIR,I1,XI,VI,FRELP,FRELD,VWHOLE)
*
*
*       Relativistic perturbation on KS pair.
*       ------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      REAL*8 XI(6),VI(6),FRELP(6),FRELD(6),TF(3),TD(3),RELACC(6)
      REAL*8 CL, CL2, CL4, CL5, NX,NY,NZ,NV1,NV2,NA1,NA2
      REAL*8 V12,V22,V1V2, V1A1,V1A2,V2A1,V2A2, RIJ, RIJ2
      REAL*8 V1_4,V2_5,V3_6,V4_1,V5_2,V6_3,V1_V22
      REAL*8 A1_4,A2_5,A3_6,A4_1,A5_2,A6_3
      REAL*8 A2R2N,A2R2V,A2R3N,A2DR2N,A2DR2V,A2DR2A
      REAL*8 A2DR3N,A2DR3V,A2DR4N,A2DR4V,A2DR5V
      REAL*8 A4R2N,A4R2V,A4R3N,A4R3V,A4R4N, A4DR2N,A4DR2V
      REAL*8 A4DR2A,A4DR3N,A4DR3V,A4DR3A,A4DR4N,A4DR4V,A4DR5N,A4DR5V
      REAL*8 A5R3N,A5R3V,A5R4N,A5R4V,A5DR3N,A5DR3V,A5DR3A,A5DR4N
      REAL*8 A5DR4A,A5DR5N,A5DR5V,AX5,AY5,AZ5
      REAL*8 CMA(3),CMVI(3), FRELP2_5(6)
      REAL*8 POWER1,POWER2,RELENERGY,EKIN,EPOTEN,ENEWT
      INTEGER IIIk
      REAL*8 NORMV(3),XNOW(3),N3_N1,N2_N1,COSA,COSB,SINA,SINB
      REAL*8 NETA1,NKSI1,NETA2,NKSI2,SEMIA,SEMIB,KSECC,KSH,KSIMPMOM2

      I2 = I1 + 1

         FRELP(1) = 0.0
         FRELP(2) = 0.0
         FRELP(3) = 0.0

	 FRELD(1) = 0.0
         FRELD(2) = 0.0
         FRELD(3) = 0.0

          RIJ2 = (XI(1)-XI(4))*(XI(1)-XI(4))+(XI(2)-XI(5))*(XI(2)-XI(5))
     &     +(XI(3)-XI(6))*(XI(3)-XI(6))
          RIJ = SQRT(RIJ2)

          NX = (XI(1)-XI(4))/RIJ
          NY = (XI(2)-XI(5))/RIJ
          NZ = (XI(3)-XI(6))/RIJ

         CMVI(1) =(BODY(I1)*VI(1)+BODY(I2)*VI(4))/(BODY(I1)+BODY(I2))
         CMVI(2) =(BODY(I1)*VI(2)+BODY(I2)*VI(5))/(BODY(I1)+BODY(I2))
         CMVI(3) =(BODY(I1)*VI(3)+BODY(I2)*VI(6))/(BODY(I1)+BODY(I2))

         VI(1) = VI(1) - CMVI(1)
         VI(2) = VI(2) - CMVI(2)
         VI(3) = VI(3) - CMVI(3)
         VI(4) = VI(4) - CMVI(1)
         VI(5) = VI(5) - CMVI(2)
         VI(6) = VI(6) - CMVI(3)

         RELACC(1) = FRELP(1) - NX*BODY(I2)/RIJ2
         RELACC(2) = FRELP(2) - NY*BODY(I2)/RIJ2
         RELACC(3) = FRELP(3) - NZ*BODY(I2)/RIJ2
         RELACC(4) = FRELP(4) + NX*BODY(I1)/RIJ2
         RELACC(5) = FRELP(5) + NY*BODY(I1)/RIJ2
         RELACC(6) = FRELP(6) + NZ*BODY(I1)/RIJ2


         CMA(1) =(BODY(I1)*RELACC(1)+BODY(I2)*
&	  RELACC(4))/(BODY(I1)+BODY(I2))
         CMA(2) =(BODY(I1)*RELACC(2)+BODY(I2)*
&	  RELACC(5))/(BODY(I1)+BODY(I2))
         CMA(3) =(BODY(I1)*RELACC(3)+BODY(I2)*
&	  RELACC(6))/(BODY(I1)+BODY(I2))




         FRELP(1) =  - CMA(1)
         FRELP(2) =  - CMA(2)
         FRELP(3) =  - CMA(3)
         FRELP(4) =  - CMA(1)
         FRELP(5) =  - CMA(2)
         FRELP(6) =  - CMA(3)

* Speed of light CL and its powers
      CL = CLIGHT
      CL2 = CL*CL
      CL4 = CL2*CL2
      CL5 = CL4*CL


*******************************************************************************************
*Notation: AiRj(N/V), ith term in summ of acceleration, j r^{-j}, N/V term is multiplied by n or v1-v2 vector.


*******************************************************************************************
******************************  Perturbation on I1 body starts  ***************************
*******************************************************************************************
          V12 =  VI(1)*VI(1)+VI(2)*VI(2)+VI(3)*VI(3)
          V22 =  VI(4)*VI(4)+VI(5)*VI(5)+VI(6)*VI(6)
          V1V2 = VI(1)*VI(4)+VI(2)*VI(5)+VI(3)*VI(6)


	  NV1 = NX*VI(1)+NY*VI(2)+NZ*VI(3)
          NV2 = NX*VI(4)+NY*VI(5)+NZ*VI(6)
          NA1 =  NX*RELACC(1)+NY*RELACC(2)+NZ*RELACC(3)
          NA2 =  NX*RELACC(4)+NY*RELACC(5)+NZ*RELACC(6)
          V1A1 = VI(1)*RELACC(1)+VI(2)*RELACC(2)+VI(3)*RELACC(3)
          V1A2 = VI(1)*RELACC(4)+VI(2)*RELACC(5)+VI(3)*RELACC(6)
          V2A1 = VI(4)*RELACC(1)+VI(5)*RELACC(2)+VI(6)*RELACC(3)
          V2A2 = VI(4)*RELACC(4)+VI(5)*RELACC(5)+VI(6)*RELACC(6)



          V1_4 = VI(1)-VI(4)
          V2_5 = VI(2)-VI(5)
          V3_6 = VI(3)-VI(6)
          A1_4 = RELACC(1)-RELACC(4)
          A2_5 = RELACC(2)-RELACC(5)
          A3_6 = RELACC(3)-RELACC(6)

          V1_V22 = V1_4*V1_4+V2_5*V2_5+V3_6*V3_6

          VWHOLE = MAX(SQRT(V12),SQRT(V22))
*********************************   A2, A2DOT  START  *************************************

          A2R2N = BODY(I2)*(-V12-2.0*V22+4.0*V1V2+1.5*NV2*NV2)
          A2R2V = BODY(I2)*(4.0*NV1-3.0*NV2)
          A2R3N = BODY(I2)*(5.0*BODY(I1)+4.0*BODY(I2))

          A2DR2N = BODY(I2)*(-2.0*V1A1+4.0*V1A2+4.0*V2A1-4.0*V2A2
     &     +3.0*NV2*NA2)
          A2DR2V = BODY(I2)*(4.0*NA1-3.0*NA2)
          A2DR2A = BODY(I2)*(4.0*NV1-3.0*NV2)
          A2DR3N = BODY(I2)*(-3.0*(NV2-NV1)*(V12+2.0*V22-4.0*V1V2)+
     &     3.0*NV2*(V1V2-V22)+7.5*NV2*NV2*(NV2-NV1))
          A2DR3V = BODY(I2)*(3.0*V12+V22-3.0*V1V2
     &     +(NV2-NV1)*(12.0*NV1-9.0*NV2)+1.5*NV2*NV2)
          A2DR4N = BODY(I2)*(4.0*(5.0*BODY(I1)+4.0*BODY(I2))*(NV2-NV1))
          A2DR4V = BODY(I2)*(5.0*BODY(I1)+4.0*BODY(I2))


          FRELP(1) = FRELP(1) + (A2R2N*NX + A2R2V*V1_4
     &     + A2R3N*NX/RIJ)/RIJ2/CL2
          FRELP(2) = FRELP(2) + (A2R2N*NY + A2R2V*V2_5
     &     + A2R3N*NY/RIJ)/RIJ2/CL2
          FRELP(3) = FRELP(3) + (A2R2N*NZ + A2R2V*V3_6
     &     + A2R3N*NZ/RIJ)/RIJ2/CL2

          FRELD(1) = (A2DR2N*NX+A2DR2V*V1_4
     &    +A2DR2A*A1_4+(A2DR3N*NX+A2DR3V*V1_4+
     &    (A2DR4N*NX+A2DR4V*V1_4)/RIJ)/RIJ )/RIJ2/CL2

          FRELD(2) = (A2DR2N*NY+A2DR2V*V2_5
     &    +A2DR2A*A2_5+(A2DR3N*NY+A2DR3V*V2_5+
     &    (A2DR4N*NY+A2DR4V*V2_5)/RIJ)/RIJ )/RIJ2/CL2

          FRELD(3) =  (A2DR2N*NZ+A2DR2V*V3_6
     &    +A2DR2A*A3_6+(A2DR3N*NZ+A2DR3V*V3_6+
     &    (A2DR4N*NZ+A2DR4V*V3_6)/RIJ)/RIJ )/RIJ2/CL2


**********************************   A2, A2DOT  END    *************************************
**********************************   A4, A4DOT  START  *************************************

          A4R2N = BODY(I2)*(-2.0*V22*V22+4.0*V22*V1V2-2.0*V1V2*V1V2
     &     +NV2*NV2*(1.5*V12+4.5*V22-6.0*V1V2-1.875*NV2*NV2))
          A4R2V = BODY(I2)*(V12*NV2+4.0*V22*NV1-5.0*V22*NV2
     &     +4.0*V1V2*(NV2-NV1)-6.0*NV1*NV2*NV2+4.5*NV2*NV2*NV2)
          A4R3N = BODY(I2)*(BODY(I1)*(-3.75*V12+1.25*V22-2.5*V1V2
     &     +19.5*NV1*NV1-39.0*NV1*NV2+8.5*NV2*NV2)
     &     +BODY(I2)*(4.0*V22-8.0*V1V2+2.0*NV1*NV1-4.0*NV1*NV2
     &     -6.0*NV2*NV2))
          A4R3V = BODY(I2)*(BODY(I1)*(-15.75*NV1+13.75*NV2)
     &     -2.0*BODY(I2)*(NV1+NV2))
          A4R4N = BODY(I2)*(-14.25*BODY(I1)*BODY(I1)
     &     -9.0*BODY(I2)*BODY(I2)-34.5*BODY(I1)*BODY(I2))

          A4DR2N = BODY(I2)*(8.0*V2A2*(V1V2-V22)+4.0*V22*(V2A1+V1A2)
     &     -4.0*V1V2*(V2A1+V1A2)+3.0*V1A1*NV2*NV2+3.0*V12*NV2*NA2
     &     +9.0*V2A2*NV2*NV2+9.0*V22*NV2*NA2-6.0*NV2*NV2*(V2A1+V1A2)
     &     -12.0*V1V2*NV2*NA2-7.5*NV2*NV2*NV2*NA2)

          A4DR2V = BODY(I2)*(2.0*V1A1*NV2+V12*NA2+8.0*V2A2*NV1
     &    +4.0*V22*NA1-10.0*V2A2*NV2-5.0*V22*NA2
     &    +4.0*(V2A1+V1A2)*(NV2-NV1)+4.0*V1V2*(NA2-NA1)
     &    -6.0*NV2*NV2*NA1-12.0*NV1*NV2*NA2+13.5*NV2*NV2*NA2)

           A4DR2A = BODY(I2)*(V12*NV2+4.0*V22*NV1-5.0*V22*NV2
     &      +4.0*V1V2*(NV2-NV1)-6.0*NV1*NV2*NV2+4.5*NV2*NV2*NV2)

           A4DR3N = BODY(I2)*(-6.0*V22*V22*(NV2-NV1)
     &      +12.0*V22*V1V2*(NV2-NV1)-6.0*V1V2*V1V2*(NV2-NV1)
     &      +3.0*V12*NV2*(V1V2-V22)+7.5*V12*NV2*NV2*(NV2-NV1)
     &      +9.0*V22*NV2*(V1V2-V22)+22.5*V22*NV2*NV2*(NV2-NV1)
     &      -12.0*V1V2*NV2*(V1V2-V22)-30.0*V1V2*NV2*NV2*(NV2-NV1)
     &      -7.5*NV2*NV2*NV2*(V1V2-V22)-13.125*NV2*NV2*NV2*NV2*(NV2-NV1)
     &      +BODY(I1)*(-7.5*V1A1+2.5*V2A2-2.5*V2A1-2.5*V1A2
     &      +39.0*NV1*(NA1-NA2)-39.0*NV2*NA1
     &      +17.0*NV2*NA2)+BODY(I2)*(8.0*V2A2-8.0*V2A1
     &      -8.0*V1A2+4.0*NV1*(NA1-NA2)-4.0*NV2*NA1-12.0*NV2*NA2))

           A4DR3V = BODY(I2)*(-2.0*V22*V22+4.0*V22*V1V2-2.0*V1V2*V1V2
     &      +1.5*V12*NV2*NV2+4.5*V22*NV2*NV2-6.0*V1V2*NV2*NV2
     &      -1.875*NV2*NV2*NV2*NV2+V12*(V1V2-V22)+3.0*V12*NV2*(NV2-NV1)
     &      +4.0*V22*(V12-V1V2)+12.0*V22*NV1*(NV2-NV1)
     &      -5.0*V22*(V1V2-V22)
     &      -15.0*V22*NV2*(NV2-NV1)-4.0*V1V2*(V12-V1V2)
     &      -12.0*V1V2*V1V2*NV1*(NV2-NV1)+4.0*V1V2*(V1V2-V22)
     &      +12.0*V1V2*V1V2*NV2*(NV2-NV1)-6.0*NV2*NV2*(V12-V1V2)
     &      -12.0*NV1*NV2*(V1V2-V22)-30.0*NV1*NV2*NV2*(NV2-NV1)
     &      +13.5*NV2*NV2*(V1V2-V22)+22.5*NV2*NV2*NV2*(NV2-NV1)
     &      +NA1*(-15.75*BODY(I1)-2.0*BODY(I2))
     &      +NA2*(13.75*BODY(I1)-2.0*BODY(I2)))

           A4DR3A = BODY(I2)*(NV1*(-15.75*BODY(I1)-2.0*BODY(I2))
     &      +NV2*(13.75*BODY(I1)-2.0*BODY(I2)))

           A4DR4N = BODY(I2)*(BODY(I1)*(-15.0*V12*(NV2-NV1)
     &      +5.0*V22*(NV2-NV1)-10.0*V1V2*(NV2-NV1)
     &      +39.0*NV1*(V12-V1V2)+117.0*NV1*NV1*(NV2-NV1)
     &      -39.0*NV1*(V1V2-V22)-39.0*NV2*(V12-V1V2)
     &      -234.0*NV1*NV2*(NV2-NV1)+17.0*NV2*(V1V2-V22)
     &      +51.0*NV2*NV2*(NV2-NV1))+BODY(I2)*(+16.0*V22*(NV2-NV1)
     &      -32.0*V1V2*(NV2-NV1)+4.0*NV1*(V12-V1V2)
     &      +12.0*NV1*NV1*(NV2-NV1)
     &      -4.0*NV1*(V1V2-V22)-4.0*NV2*(V12-V1V2)
     &      -24.0*NV1*NV2*(NV2-NV1)
     &      -12.0*NV2*(V1V2-V22)-36.0*NV2*NV2*(NV2-NV1)))

           A4DR4V = BODY(I2)*(BODY(I1)*(-3.75*V12+1.25*V22-2.5*V1V2
     &      +19.5*NV1*NV1-39.0*NV1*NV2+8.5*NV2*NV2)+BODY(I2)*(4.0*V22
     &      -8.0*V1V2+2.0*NV1*NV1-4.0*NV1*NV2-6.0*NV2*NV2)
     &      +(-15.75*BODY(I1)-2.0*BODY(I2))*(V12-V1V2
     &      -4.0*NV1*(NV2-NV1))+(13.75*BODY(I1)-2.0*BODY(I2))
     &      *(V1V2-V22-4.0*NV2*(NV2-NV1)))

           A4DR5N = BODY(I2)*(-71.25*BODY(I1)*BODY(I1)
     &      -45.0*BODY(I2)*BODY(I2)-172.5*BODY(I1)*BODY(I2))*(NV2-NV1)

           A4DR5V = BODY(I2)*(-14.25*BODY(I1)*BODY(I1)
     &      -9.0*BODY(I2)*BODY(I2)-34.5*BODY(I1)*BODY(I2))


          FRELP(1) = FRELP(1) + (A4R2N*NX+A4R2V*V1_4+
     &     (A4R3N*NX+A4R3V*V1_4+A4R4N*NX/RIJ)/RIJ)/RIJ2/CL4

          FRELP(2) = FRELP(2) + (A4R2N*NY+A4R2V*V2_5+
     &     (A4R3N*NY+A4R3V*V2_5+A4R4N*NY/RIJ)/RIJ)/RIJ2/CL4

          FRELP(3) = FRELP(3) + (A4R2N*NZ+A4R2V*V3_6+
     &     (A4R3N*NZ+A4R3V*V3_6+A4R4N*NZ/RIJ)/RIJ)/RIJ2/CL4


          FRELD(1) = FRELD(1) + (A4DR2N*NX+A4DR2V*V1_4
     &     +A4DR2A*A1_4+(A4DR3N*NX+A4DR3V*V1_4
     &     +A4DR3A*A1_4+(A4DR4N*NX+A4DR4V*V1_4
     &     +(A4DR5N*NX+A4DR5V*V1_4)/RIJ)/RIJ)/RIJ)/RIJ2/CL4

          FRELD(2) = FRELD(2) + (A4DR2N*NY+A4DR2V*V2_5
     &     +A4DR2A*A2_5+(A4DR3N*NY+A4DR3V*V2_5
     &     +A4DR3A*A2_5+(A4DR4N*NY+A4DR4V*V2_5
     &     +(A4DR5N*NY+A4DR5V*V2_5)/RIJ)/RIJ)/RIJ)/RIJ2/CL4

          FRELD(3) = FRELD(3) + (A4DR2N*NZ+A4DR2V*V3_6
     &     +A4DR2A*A3_6+(A4DR3N*NZ+A4DR3V*V3_6
     &     +A4DR3A*A3_6+(A4DR4N*NZ+A4DR4V*V3_6
     &     +(A4DR5N*NZ+A4DR5V*V3_6)/RIJ)/RIJ)/RIJ)/RIJ2/CL4

**********************************   A4, A4DOT  END    *************************************
**********************************   A5, A5DOT  START  *************************************

          A5R3N = 0.8*BODY(I2)*BODY(I1)*(3.0*(NV1-NV2)*V1_V22)
          A5R3V = 0.8*BODY(I2)*BODY(I1)*(-V1_V22)
          A5R4N = 0.8*BODY(I2)*BODY(I1)*((NV1-NV2)*(-6.0*BODY(I1)
     &      +17.33333333333*BODY(I2)))
          A5R4V = 0.8*BODY(I2)*BODY(I1)*(2.0*BODY(I1)-8.0*BODY(I2))

               A5DR3N = 0.8*BODY(I2)*BODY(I1)*(3.0*V1_V22*(NA1-NA2)
     &          +6.0*(NV1-NV2)*(V1A1+V2A2-V2A1-V1A2))
               A5DR3V = 0.8*BODY(I2)*BODY(I1)*
     &          (-2.0*(V1A1+V2A2-V2A1-V1A2))
               A5DR3A = 0.8*BODY(I2)*BODY(I1)*(-V1_V22)
               A5DR4N = 0.8*BODY(I2)*BODY(I1)*(3.0*V1_V22*
     &          (V12+V22-2.0*V1V2)
     &          -15.0*V1_V22*(NV1-NV2)*(NV1-NV2)+(17.333333333*BODY(I2)
     &          -6.0*BODY(I1))*(NA1-NA2))
               A5DR4A = 0.8*BODY(I2)*BODY(I1)
     &          *(2.0*BODY(I1)-8.0*BODY(I2))
               A5DR5N = 0.8*BODY(I2)*BODY(I1)*((17.33333333333333
     &          *BODY(I2)-6.0*BODY(I1))*(V12-2.0*V1V2+V22
     &          -6.0*(NV2-NV1)*(NV2-NV1)))
               A5DR5V = 0.8*BODY(I2)*BODY(I1)*((-14.0*BODY(I1)
     &          +49.3333333333333*BODY(I2))*(NV1-NV2))
               A5DR4V = 0.8*BODY(I2)*BODY(I1)*(6.0*V1_V22*(NV1-NV2))


          FRELP2_5(1) = (A5R3N*NX + A5R3V*V1_4
     &     + (A5R4N*NX+A5R4V*V1_4)/RIJ)/RIJ2/RIJ/CL5
          FRELP2_5(2) = (A5R3N*NY + A5R3V*V2_5
     &     + (A5R4N*NY+A5R4V*V2_5)/RIJ)/RIJ2/RIJ/CL5
          FRELP2_5(3) = (A5R3N*NZ + A5R3V*V3_6
     &     + (A5R4N*NZ+A5R4V*V3_6)/RIJ)/RIJ2/RIJ/CL5

          FRELP(1) = FRELP(1) + FRELP2_5(1) 
          FRELP(2) = FRELP(2) + FRELP2_5(2) 
          FRELP(3) = FRELP(3) + FRELP2_5(3) 



          FRELD(1) = FRELD(1) + (A5DR3N*NX+A5DR3V*V1_4
     &    +A5DR3A*A1_4+(A5DR4N*NX
     &    +A5DR4A*A1_4+(A5DR5N*NX+A5DR5V*V1_4)
     &    /RIJ)/RIJ)/RIJ/RIJ2/CL5
     &    + A5DR4V*V1_4/RIJ2/RIJ2/CL5
          FRELD(2) = FRELD(2) + (A5DR3N*NY+A5DR3V*V2_5
     &    +A5DR3A*A2_5+(A5DR4N*NY
     &    +A5DR4A*A2_5+(A5DR5N*NY+A5DR5V*V2_5)
     &    /RIJ)/RIJ)/RIJ/RIJ2/CL5
     &    + A5DR4V*V2_5/RIJ2/RIJ2/CL5
          FRELD(3) = FRELD(3) + (A5DR3N*NZ+A5DR3V*V3_6
     &    +A5DR3A*A3_6+(A5DR4N*NZ
     &    +A5DR4A*A3_6+(A5DR5N*NZ+A5DR5V*V3_6)
     &    /RIJ)/RIJ)/RIJ/RIJ2/CL5
     &    + A5DR4V*V3_6/RIJ2/RIJ2/CL5


*********************************   A5, A5DOT  END    *************************************

*******************************************************************************************
******************************  Perturbation on I1 body ends    ***************************
*******************************************************************************************
*_________________________________________________________________________________________*
*******************************************************************************************
******************************  Perturbation on I2 body starts  ***************************
*******************************************************************************************
* Exchanging 1<->4, 2<->5, 3<->6

          V22 = VI(4)*VI(4)+VI(5)*VI(5)+VI(6)*VI(6)
          V12 = VI(1)*VI(1)+VI(2)*VI(2)+VI(3)*VI(3)

          NX = -NX
          NY = -NY
          NZ = -NZ
          NV2 = NX*VI(4)+NY*VI(5)+NZ*VI(6)
          NV1 = NX*VI(1)+NY*VI(2)+NZ*VI(3)
          NA2 =  NX*RELACC(4)+NY*RELACC(5)+NZ*RELACC(6)
          NA1 =  NX*RELACC(1)+NY*RELACC(2)+NZ*RELACC(3)
          V2A2 = VI(4)*RELACC(4)+VI(5)*RELACC(5)+VI(6)*RELACC(6)
          V2A1 = VI(4)*RELACC(1)+VI(5)*RELACC(2)+VI(6)*RELACC(3)
          V1A2 = VI(1)*RELACC(4)+VI(2)*RELACC(5)+VI(3)*RELACC(6)
          V1A1 = VI(1)*RELACC(1)+VI(2)*RELACC(2)+VI(3)*RELACC(3)

          V4_1 = VI(4)-VI(1)
          V5_2 = VI(5)-VI(2)
          V6_3 = VI(6)-VI(3)
          A4_1 = RELACC(4)-RELACC(1)
          A5_2 = RELACC(5)-RELACC(2)
          A6_3 = RELACC(6)-RELACC(3)

          V2_V12 = V4_1*V4_1+V5_2*V5_2+V6_3*V6_3

**********************************   A2, A2DOT  START  *************************************
          A2R2N = BODY(I1)*(-V22-2.0*V12+4.0*V1V2+1.5*NV1*NV1)
          A2R2V = BODY(I1)*(4.0*NV2-3.0*NV1)

          A2R3N = BODY(I1)*(5.0*BODY(I2)+4.0*BODY(I1))

          A2DR2N = BODY(I1)*(-2.0*V2A2+4.0*V2A1+4.0*V1A2-4.0*V1A1
     &     +3.0*NV1*NA1)
          A2DR2V = BODY(I1)*(4.0*NA2-3.0*NA1)
          A2DR2A = BODY(I1)*(4.0*NV2-3.0*NV1)
          A2DR3N = BODY(I1)*(-3.0*(NV1-NV2)*(V22+2.0*V12-4.0*V1V2)+
     &     3.0*NV1*(V1V2-V12)+7.5*NV1*NV1*(NV1-NV2))
          A2DR3V = BODY(I1)*(3.0*V22+V12-3.0*V1V2
     &     +(NV1-NV2)*(12.0*NV2-9.0*NV1)+1.5*NV1*NV1)
          A2DR4N = BODY(I1)*(4.0*(5.0*BODY(I2)+4.0*BODY(I1))*(NV1-NV2))
          A2DR4V = BODY(I1)*(5.0*BODY(I2)+4.0*BODY(I1))


          FRELP(4) = FRELP(4) +  (A2R2N*NX + A2R2V*V4_1
     &     + A2R3N*NX/RIJ)/RIJ2/CL2
          FRELP(5) = FRELP(5) +  (A2R2N*NY + A2R2V*V5_2
     &     + A2R3N*NY/RIJ)/RIJ2/CL2
          FRELP(6) = FRELP(6) +  (A2R2N*NZ + A2R2V*V6_3
     &     + A2R3N*NZ/RIJ)/RIJ2/CL2


          FRELD(4) =  (A2DR2N*NX+A2DR2V*V4_1
     &    +A2DR2A*A4_1+(A2DR3N*NX+A2DR3V*V4_1+
     &    (A2DR4N*NX+A2DR4V*V4_1)/RIJ)/RIJ )/RIJ2/CL2

          FRELD(5) =  (A2DR2N*NY+A2DR2V*V5_2
     &    +A2DR2A*A5_2+(A2DR3N*NY+A2DR3V*V5_2+
     &    (A2DR4N*NY+A2DR4V*V5_2)/RIJ)/RIJ )/RIJ2/CL2

          FRELD(6) =  (A2DR2N*NZ+A2DR2V*V6_3
     &    +A2DR2A*A6_3+(A2DR3N*NZ+A2DR3V*V6_3+
     &    (A2DR4N*NZ+A2DR4V*V6_3)/RIJ)/RIJ )/RIJ2/CL2


*********************************   A2, A2DOT  END    *************************************
*********************************   A4, A4DOT  START  *************************************

          A4R2N = BODY(I1)*(-2.0*V12*V12+4.0*V12*V1V2-2.0*V1V2*V1V2
     &     +NV1*NV1*(1.5*V22+4.5*V12-6.0*V1V2-1.875*NV1*NV1))
          A4R2V = BODY(I1)*(V22*NV1+4.0*V12*NV2-5.0*V12*NV1
     &     +4.0*V1V2*(NV1-NV2)-6.0*NV2*NV1*NV1+4.5*NV1*NV1*NV1)
          A4R3N = BODY(I1)*(BODY(I2)*(-3.75*V22+1.25*V12-2.5*V1V2
     &     +19.5*NV2*NV2-39.0*NV2*NV1+8.5*NV1*NV1)
     &     +BODY(I1)*(4.0*V12-8.0*V1V2+2.0*NV2*NV2-4.0*NV2*NV1
     &     -6.0*NV1*NV1))
          A4R3V = BODY(I1)*(BODY(I2)*(-15.75*NV2+13.75*NV1)
     &     -2.0*BODY(I1)*(NV2+NV1))
          A4R4N = BODY(I1)*(-14.25*BODY(I2)*BODY(I2)
     &     -9.0*BODY(I1)*BODY(I1)-34.5*BODY(I2)*BODY(I1))

          A4DR2N = BODY(I1)*(8.0*V1A1*(V1V2-V12)+4.0*V12*(V1A2+V2A1)
     &     -4.0*V1V2*(V1A2+V2A1)+3.0*V2A2*NV1*NV1+3.0*V22*NV1*NA1
     &     +9.0*V1A1*NV1*NV1+9.0*V12*NV1*NA1-6.0*NV1*NV1*(V1A2+V2A1)
     &     -12.0*V1V2*NV1*NA1-7.5*NV1*NV1*NV1*NA1)

          A4DR2V = BODY(I1)*(2.0*V2A2*NV1+V22*NA1+8.0*V1A1*NV2
     &    +4.0*V12*NA2-10.0*V1A1*NV1-5.0*V12*NA1
     &    +4.0*(V1A2+V2A1)*(NV1-NV2)+4.0*V1V2*(NA1-NA2)
     &    -6.0*NV1*NV1*NA2-12.0*NV2*NV1*NA1+13.5*NV1*NV1*NA1)

           A4DR2A = BODY(I1)*(V22*NV1+4.0*V12*NV2-5.0*V12*NV1
     &      +4.0*V1V2*(NV1-NV2)-6.0*NV2*NV1*NV1+4.5*NV1*NV1*NV1)

           A4DR3N = BODY(I1)*(-6.0*V12*V12*(NV1-NV2)
     &      +12.0*V12*V1V2*(NV1-NV2)-6.0*V1V2*V1V2*(NV1-NV2)
     &      +3.0*V22*NV1*(V1V2-V12)+7.5*V22*NV1*NV1*(NV1-NV2)
     &      +9.0*V12*NV1*(V1V2-V12)+22.5*V12*NV1*NV1*(NV1-NV2)
     &      -12.0*V1V2*NV1*(V1V2-V12)-30.0*V1V2*NV1*NV1*(NV1-NV2)
     &      -7.5*NV1*NV1*NV1*(V1V2-V12)
     &      -13.125*NV1*NV1*NV1*NV1*(NV1-NV2)
     &      +BODY(I2)*(-7.5*V2A2+2.5*V1A1-2.5*V1A2-2.5*V2A1
     &      +39.0*NV2*(NA2-NA1)-39.0*NV1*NA2
     &      +17.0*NV1*NA1)+BODY(I1)*(8.0*V1A1-8.0*V1A2
     &      -8.0*V2A1+4.0*NV2*(NA2-NA1)-4.0*NV1*NA2-12.0*NV1*NA1))

           A4DR3V = BODY(I1)*(-2.0*V12*V12+4.0*V12*V1V2-2.0*V1V2*V1V2
     &      +1.5*V22*NV1*NV1+4.5*V12*NV1*NV1-6.0*V1V2*NV1*NV1
     &      -1.875*NV1*NV1*NV1*NV1+V22*(V1V2-V12)+3.0*V22*NV1*(NV1-NV2)
     &      +4.0*V12*(V22-V1V2)+12.0*V12*NV2*(NV1-NV2)
     &      -5.0*V12*(V1V2-V12)
     &      -15.0*V12*NV1*(NV1-NV2)-4.0*V1V2*(V22-V1V2)
     &      -12.0*V1V2*V1V2*NV2*(NV1-NV2)+4.0*V1V2*(V1V2-V12)
     &      +12.0*V1V2*V1V2*NV1*(NV1-NV2)-6.0*NV1*NV1*(V22-V1V2)
     &      -12.0*NV2*NV1*(V1V2-V12)-30.0*NV2*NV1*NV1*(NV1-NV2)
     &      +13.5*NV1*NV1*(V1V2-V12)+22.5*NV1*NV1*NV1*(NV1-NV2)
     &      +NA2*(-15.75*BODY(I2)-2.0*BODY(I1))
     &      +NA1*(13.75*BODY(I2)-2.0*BODY(I1)))

           A4DR3A = BODY(I1)*(NV2*(-15.75*BODY(I2)-2.0*BODY(I1))
     &      +NV1*(13.75*BODY(I2)-2.0*BODY(I1)))

           A4DR4N = BODY(I1)*(BODY(I2)*(-15.0*V22*(NV1-NV2)
     &      +5.0*V12*(NV1-NV2)-10.0*V1V2*(NV1-NV2)
     &      +39.0*NV2*(V22-V1V2)+117.0*NV2*NV2*(NV1-NV2)
     &      -39.0*NV2*(V1V2-V12)-39.0*NV1*(V22-V1V2)
     &      -234.0*NV2*NV1*(NV1-NV2)+17.0*NV1*(V1V2-V12)
     &      +51.0*NV1*NV1*(NV1-NV2))+BODY(I1)*(+16.0*V12*(NV1-NV2)
     &      -32.0*V1V2*(NV1-NV2)+4.0*NV2*(V22-V1V2)
     &      +12.0*NV2*NV2*(NV1-NV2)-24.0*NV2*NV1*(NV1-NV2)
     &      -4.0*NV2*(V1V2-V12)-4.0*NV1*(V22-V1V2)
     &      -12.0*NV1*(V1V2-V12)-36.0*NV1*NV1*(NV1-NV2)))

           A4DR4V = BODY(I1)*(BODY(I2)*(-3.75*V22+1.25*V12-2.5*V1V2
     &      +19.5*NV2*NV2-39.0*NV2*NV1+8.5*NV1*NV1)+BODY(I1)*(4.0*V12
     &      -8.0*V1V2+2.0*NV2*NV2-4.0*NV2*NV1-6.0*NV1*NV1)
     &      +(-15.75*BODY(I2)-2.0*BODY(I1))*(V22-V1V2
     &      -4.0*NV2*(NV1-NV2))+(13.75*BODY(I2)-2.0*BODY(I1))
     &      *(V1V2-V12-4.0*NV1*(NV1-NV2)))

           A4DR5N = BODY(I1)*(-71.25*BODY(I2)*BODY(I2)
     &      -45.0*BODY(I1)*BODY(I1)-172.5*BODY(I2)*BODY(I1))*(NV1-NV2)
           A4DR5V = BODY(I1)*(-14.25*BODY(I2)*BODY(I2)
     &      -9.0*BODY(I1)*BODY(I1)-34.5*BODY(I2)*BODY(I1))


          FRELP(4) = FRELP(4) + (A4R2N*NX+A4R2V*V4_1+
     &     (A4R3N*NX+A4R3V*V4_1+A4R4N*NX/RIJ)/RIJ)/RIJ2/CL4

          FRELP(5) = FRELP(5) + (A4R2N*NY+A4R2V*V5_2+
     &     (A4R3N*NY+A4R3V*V5_2+A4R4N*NY/RIJ)/RIJ)/RIJ2/CL4

          FRELP(6) = FRELP(6) + (A4R2N*NZ+A4R2V*V6_3+
     &     (A4R3N*NZ+A4R3V*V6_3+A4R4N*NZ/RIJ)/RIJ)/RIJ2 /CL4


          FRELD(4) = FRELD(4) + (A4DR2N*NX+A4DR2V*V4_1
     &     +A4DR2A*A4_1+(A4DR3N*NX+A4DR3V*V4_1
     &     +A4DR3A*A4_1+(A4DR4N*NX+A4DR4V*V4_1
     &     +(A4DR5N*NX+A4DR5V*V4_1)/RIJ)/RIJ)/RIJ)/RIJ2/CL4

          FRELD(5) = FRELD(5) + (A4DR2N*NY+A4DR2V*V5_2
     &     +A4DR2A*A5_2+(A4DR3N*NY+A4DR3V*V5_2
     &     +A4DR3A*A5_2+(A4DR4N*NY+A4DR4V*V5_2
     &     +(A4DR5N*NY+A4DR5V*V5_2)/RIJ)/RIJ)/RIJ)/RIJ2/CL4

          FRELD(6) = FRELD(6) + (A4DR2N*NZ+A4DR2V*V6_3
     &     +A4DR2A*A6_3+(A4DR3N*NZ+A4DR3V*V6_3
     &     +A4DR3A*A6_3+(A4DR4N*NZ+A4DR4V*V6_3
     &     +(A4DR5N*NZ+A4DR5V*V6_3)/RIJ)/RIJ)/RIJ)/RIJ2/CL4

**********************************   A4, A4DOT  END    *************************************
**********************************   A5, A5DOT  START  *************************************

          A5R3N = 0.8*BODY(I1)*BODY(I2)*(3.0*(NV2-NV1)*V2_V12)
          A5R3V = 0.8*BODY(I1)*BODY(I2)*(-V2_V12)
          A5R4N = 0.8*BODY(I1)*BODY(I2)*((NV2-NV1)*(-6.0*BODY(I2)
     &      +17.33333333333*BODY(I1)))
          A5R4V = 0.8*BODY(I1)*BODY(I2)*(2.0*BODY(I2)-8.0*BODY(I1))

               A5DR3N = 0.8*BODY(I1)*BODY(I2)*(3.0*V2_V12*(NA2-NA1)
     &          +6.0*(NV2-NV1)*(V2A2+V1A1-V1A2-V2A1))
               A5DR3V = 0.8*BODY(I1)*BODY(I2)*
     &          (-2.0*(V2A2+V1A1-V1A2-V2A1))
               A5DR3A = 0.8*BODY(I1)*BODY(I2)*(-V2_V12)
               A5DR4N = 0.8*BODY(I1)*BODY(I2)*(3.0*V2_V12*
     &          (V22+V12-2.0*V1V2)
     &          -15.0*V2_V12*(NV2-NV1)*(NV2-NV1)+(17.333333333*BODY(I1)
     &          -6.0*BODY(I2))*(NA2-NA1))
               A5DR4A = 0.8*BODY(I1)*BODY(I2)
     &          *(2.0*BODY(I2)-8.0*BODY(I1))
               A5DR5N = 0.8*BODY(I1)*BODY(I2)*((17.33333333333333
     &          *BODY(I1)-6.0*BODY(I2))*(V22-2.0*V1V2+V12
     &          -6.0*(NV1-NV2)*(NV1-NV2)))
               A5DR5V = 0.8*BODY(I1)*BODY(I2)*((-14.0*BODY(I2)
     &          +49.3333333333333*BODY(I1))*(NV2-NV1))
               A5DR4V =0.8*BODY(I2)*BODY(I1)*(6.0*V1_V22*(NV2-NV1))


          FRELP2_5(4) = (A5R3N*NX + A5R3V*V4_1
     &     + (A5R4N*NX+A5R4V*V4_1)/RIJ)/RIJ2/RIJ/CL5
          FRELP2_5(5) = (A5R3N*NY + A5R3V*V5_2
     &     + (A5R4N*NY+A5R4V*V5_2)/RIJ)/RIJ2/RIJ/CL5
          FRELP2_5(6) = (A5R3N*NZ + A5R3V*V6_3
     &     + (A5R4N*NZ+A5R4V*V6_3)/RIJ)/RIJ2/RIJ/CL5


          FRELP(4) = FRELP(4) + FRELP2_5(4)
          FRELP(5) = FRELP(5) + FRELP2_5(5)
          FRELP(6) = FRELP(6) + FRELP2_5(6)



          FRELD(4) = FRELD(4) + (A5DR3N*NX+A5DR3V*V4_1
     &    +A5DR3A*A4_1+(A5DR4N*NX
     &    +A5DR4A*A4_1+(A5DR5N*NX+A5DR5V*V4_1)
     &    /RIJ)/RIJ)/RIJ/RIJ2/CL5
     &    + A5DR4V*V4_1/RIJ2/RIJ2/CL5
          FRELD(5) = FRELD(5) + (A5DR3N*NY+A5DR3V*V5_2
     &    +A5DR3A*A5_2+(A5DR4N*NY
     &    +A5DR4A*A5_2+(A5DR5N*NY+A5DR5V*V5_2)
     &    /RIJ)/RIJ)/RIJ/RIJ2/CL5
     &    + A5DR4V*V5_2/RIJ2/RIJ2/CL5
          FRELD(6) = FRELD(6) + (A5DR3N*NZ+A5DR3V*V6_3
     &    +A5DR3A*A6_3+(A5DR4N*NZ
     &    +A5DR4A*A6_3+(A5DR5N*NZ+A5DR5V*V6_3)
     &    /RIJ)/RIJ)/RIJ/RIJ2/CL5
     &    + A5DR4V*V6_3/RIJ2/RIJ2/CL5

**********************************   A5, A5DOT  END    *************************************

*******************************************************************************************
******************************  Perturbation on I2 body ends    ***************************
*******************************************************************************************

	 EKIN = BODY(I1)*V12/2.0 + BODY(I2)*V22/2.0
         EPOTEN = -BODY(I1)*BODY(I2)/RIJ
	 ENEWT(IPAIR) = EKIN + EPOTEN


         POWER1 = (FRELP(1)*VI(1) + FRELP(2)*VI(2) + FRELP(3)*VI(3))
&	 *BODY(I1)
         POWER2 = (FRELP(4)*VI(4) + FRELP(5)*VI(5) + FRELP(6)*VI(6))
&	 *BODY(I2)

         POWER1_2_5 = (FRELP2_5(1)*VI(1) + FRELP2_5(2)*VI(2) +
& FRELP2_5(3)*VI(3))*BODY(I1)
         POWER2_2_5 = (FRELP2_5(4)*VI(4) + FRELP2_5(5)*VI(5) + 
& FRELP2_5(6)*VI(6))*BODY(I2)

      	  RELENERGY2_5(IPAIR) = RELENERGY2_5(IPAIR) + 
&(POWER1_2_5 + POWER2_2_5)*RELTSTEP(IPAIR)


      	  RELENERGY(IPAIR) = RELENERGY(IPAIR) + (POWER1 + POWER2)
&	 *RELTSTEP(IPAIR)

	 IF(ILL(IPAIR).EQ.0) THEN
	  WRITE(*,112)'ENTER KSREL',TIME,NAME(I1),NAME(I2)
&,(POWER1 + POWER2)*RELTSTEP(IPAIR),ENEWT(IPAIR),SQRT(V12),SQRT(V22)
&,RIJ,(XI(KKK),KKK=1,6),(VI(KKK),KKK=1,6),IPAIR
           CALL FLUSH(6)
	 RELENERGY0(IPAIR) = RELENERGY(IPAIR)
	 RELENERGY02_5(IPAIR) = RELENERGY2_5(IPAIR)
         ENEWT0(IPAIR) = ENEWT(IPAIR)

         ILL(IPAIR) = 1
	 END IF
 112      FORMAT(A11,F35.25,2I6,17F35.25,I6)
 



          IIIk = IIIk + 1
          IF(IIIk.EQ.1000000) IIIk = 0

	  IF((IIIk/100000)*100000.EQ.IIIk) THEN
               
         RX = XI(1) - XI(4)
         RY = XI(2) - XI(5)
         RZ = XI(3) - XI(6)  
          
	 XCOMPIMP = RY*V3_6 - RZ*V2_5
	 YCOMPIMP = RZ*V1_4 - RX*V3_6
	 ZCOMPIMP = RX*V2_5 - RY*V1_4

          KSIMPMOM2 = XCOMPIMP*XCOMPIMP + YCOMPIMP*YCOMPIMP 
&+ ZCOMPIMP*ZCOMPIMP 

         KSH = 0.5*V1_V22 - (BODY(I1)+BODY(I2))/RIJ

          IF(1.0+2.0*KSH*KSIMPMOM2/(BODY(I1)+BODY(I2))
&/(BODY(I1)+BODY(I2)).GE.0.0) THEN
           KSECC = SQRT(1.0+2.0*KSH*KSIMPMOM2/(BODY(I1)+BODY(I2))
&/(BODY(I1)+BODY(I2)))
           END IF 

           SEMI = -0.5*0.002/H(IPAIR)
           ECC = 1.0 - R(IPAIR)/SEMI
           FRP=SQRT((FRELP(1)-FRELP(4))**2 + 
&(FRELP(2)-FRELP(5))**2 + (FRELP(3)-FRELP(6))**2)

           GREL = FRP*RIJ**2/(BODY(I1)+BODY(I2))


           WRITE(*,23)'KS RUNS',TIME,RIJ/
&	   (2.0*(BODY(I1)+BODY(I2))/CL2),NAME(I1),NAME(I2),
&        (XI(KKK),KKK=1,6),(VI(KKK),KKK=1,6),ENEWT(IPAIR)
&        ,RELENERGY(IPAIR),GREL,KSECC,RELENERGY2_5(IPAIR)

           CALL FLUSH(6)
	  END IF
  23      FORMAT(A7,2F25.12,2I6,17F35.22)




         VI(1) = VI(1) + CMVI(1)
         VI(2) = VI(2) + CMVI(2)
         VI(3) = VI(3) + CMVI(3)
         VI(4) = VI(4) + CMVI(1)
         VI(5) = VI(5) + CMVI(2)
         VI(6) = VI(6) + CMVI(3)



         FRELP(1) = FRELP(1) + CMA(1)
         FRELP(2) = FRELP(2) + CMA(2)
         FRELP(3) = FRELP(3) + CMA(3)
         FRELP(4) = FRELP(4) + CMA(1)
         FRELP(5) = FRELP(5) + CMA(2)
         FRELP(6) = FRELP(6) + CMA(3)



*       Set the relative perturbing force and first derivative.
   50 DO 55 K = 1,3
          FRELP(K) = FRELP(K) - FRELP(K+3)
          FRELD(K) = FRELD(K) - FRELD(K+3)
          TF(K) = 0.0D0
          TD(K) = 0.0D0
   55 CONTINUE
*

      RETURN
*
      END
