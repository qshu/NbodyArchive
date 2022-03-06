      subroutine ksrel_output(time,xi,vi,n1,n2,m1,m2,dt,np,hi,gi_pn,
     &     gi_nopn,spni1,ipair,TIME2,TNEW2)
*
*     Output for PN events
*     ------------------------------------
*
      include 'params.h'
      include 'postnewton.h'
      REAL*8 time,XI(6),VI(6),M1,M2,DT,HI,GI_PN,GI_NOPN
      INTEGER N1,N2,NP,IPAIR
      REAL*8 XX12(3),VV12(3),V1_2,V2_2,V1_22,MU,M12,ETA,RIJ,ETEMP
      REAL*8 LANG(3),L2,L2PN,ARAD,POTMIN,POTMINPN
      REAL*8 SEMI,SEMIPN,ECC2,ECC,ECCPN,ECC2PN
*      REAL*8 ECC2_FULL,ECC2_FULL_PN
      REAL*8 GE,TZ,PI,A_EIN
      REAL*8 EFULL,ECC2_FULL,LFULL(3),DL2_FULL,HH
      REAL*8 DL_FULL 
      REAL*8 VWHOLE, RP, LPN1,LPN2,LPN3 
      REAL*8 L_PN1(3),L_PN2(3),L_PN3(3),L_PN2_5(3), L_PN(3)
      REAL*8 CL2,CL4,CL5,CL6,ETA2,ETA3,RIJ2, EDOT,ADOT,TGR
      INTEGER SPNI1,SPNI2, TIME2
      REAL*8 STEP(NMAX)
**    
*      SAVE ICOUNT
*      DATA ICOUNT/0/

      SPNI2 = SPNI1 + 1
      
      CL2= CLIGHTPN*CLIGHTPN
      CL4 = CL2*CL2
      CL6 = CL4*CL2
      CL5 = CL4*CLIGHTPN
      PI = 3.141592654
      PI2 = PI*PI

      XX12(1) = XI(1)-XI(4)
      XX12(2) = XI(2)-XI(5)
      XX12(3) = XI(3)-XI(6)
      VV12(1) = VI(1)-VI(4)
      VV12(2) = VI(2)-VI(5)
      VV12(3) = VI(3)-VI(6)         
      V1_2 = VI(1)**2+VI(2)**2+VI(3)**2
      V2_2 = VI(4)**2+VI(5)**2+VI(6)**2
      V1_22 = VV12(1)**2+VV12(2)**2+VV12(3)**2
      VWHOLE = SQRT(V1_22)
      M12 = M1+M2
      MU = M1*M2/M12
      ETA = MU/M12
      ETA2 = ETA*ETA
      ETA3 = ETA2*ETA
      RIJ = SQRT(XX12(1)**2+XX12(2)**2+XX12(3)**2)
      RIJ2 = RIJ*RIJ
      RP = (XX12(1)*VV12(1)+XX12(2)*VV12(2)+XX12(3)*VV12(3))/RIJ
      MOR = M12/RIJ
      ETEMP = (0.5*V1_22 - M12/RIJ)
      EFULL = ETEMP + PREVRELENERGY 
C      ETEMP = HI
      LANG(1) = (XX12(2)*VV12(3)-XX12(3)*VV12(2))
      LANG(2) = (XX12(3)*VV12(1)-XX12(1)*VV12(3))
      LANG(3) = (XX12(1)*VV12(2)-XX12(2)*VV12(1)) 
      L2 = (LANG(1)**2+LANG(2)**2+LANG(3)**2) !/(MU)**2
         
      LPN1 = (1./CL2)*(0.5*(1.0-3.0*ETA)*V1_22+(3.0+ETA)*MOR)
      LPN2 = (1./CL4)*(3.* V1_22*V1_22/8.
     &      - 21.*ETA*V1_22*V1_22/8. + 39.*ETA2*V1_22*V1_22/8.
     &      + MOR*(-RP*RP*ETA-5*RP*RP*ETA2/2. + 7.*V1_22/2
     &      -5.*ETA*V1_22 -9.*ETA2*V1_22/2.)+(M12*MOR)*
     &      (7./2. -41.*ETA/4. + ETA2))
      LPN3 = (1./CL6)*(1./16.*( 5. - 59.*ETA +238.*ETA2 -
     &        323.*ETA3)*(VWHOLE)**6.
     &        +1./8.*( 33. - 142.*ETA +106.*ETA2 + 195.*ETA3)
     &        *(VWHOLE**4.)*MOR
     &        - 1./4.*ETA*( 12. -7.*ETA - 75.*ETA2)*(VWHOLE**2.)
     &        *RP*RP*MOR
     &        + 3./8.*ETA*(2.-2.*ETA- 11.*ETA2)*(RP)**4.*M12/MOR
     &        + 1./12.*(135.- 322.*ETA+315.*ETA2-108.*ETA3)*
     &        (VWHOLE**2.)*(MOR)**2.
     &        + 1./24.*( 12.-287.*ETA-951.*ETA2-324.*ETA3)*
     &        RP*RP*(MOR)**2.
     &        + (5./2.- 1./1120.*(20796.-1435.*PI2)*ETA -7.*ETA2+ETA3))

      L_PN1(1)   = LPN1*LANG(1)
      L_PN1(2)   = LPN1*LANG(2)
      L_PN1(3)   = LPN1*LANG(3)
      L_PN2(1)   = LPN2*LANG(1)
      L_PN2(2)   = LPN2*LANG(2)
      L_PN2(3)   = LPN2*LANG(3)
      L_PN3(1)   = LPN3*LANG(1)
      L_PN3(2)   = LPN3*LANG(2)
      L_PN3(3)   = LPN3*LANG(3)

      L_PN2_5(1) = 8./5.*1./CL5/MU *(M12**3.)*RP*ETA2/RIJ2*LANG(1)
      L_PN2_5(2) = 8./5.*1./CL5/MU *(M12**3.)*RP*ETA2/RIJ2*LANG(2)
      L_PN2_5(3) = 8./5.*1./CL5/MU *(M12**3.)*RP*ETA2/RIJ2*LANG(3)

      L_PN(1) = L_PN1(1)+L_PN2(1)+L_PN3(1)
      L_PN(2) = L_PN1(2)+L_PN2(2)+L_PN3(2)
      L_PN(3) = L_PN1(3)+L_PN2(3)+L_PN3(3)

      LFULL(1) = LANG(1)+L_PN(1)
      LFULL(2) = LANG(2)+L_PN(2)
      LFULL(3) = LANG(3)+L_PN(3) 
      DL2_FULL = LFULL(1)*LFULL(1)+LFULL(2)*LFULL(2)+LFULL(3)*LFULL(3)
      DL_FULL = SQRT(DL2_FULL)    
      HH = DL_FULL/M12  
      L2PN = L2*(1.0+LPN1/CL2)**2
*
      
      ARAD = -0.50*M12/(ETEMP)*(1.0 + (-2.0*ETEMP)/(4.0*CL2)*(-7.0+ETA)+
     &      ((-2.0*ETEMP)**2/(16.0*CL2*CL2))*(1.0+ETA*ETA+
     &      16.0/(-2.*ETEMP*HH*HH)*(-4.0+7.0*ETA))+
     &     (((-2.*ETEMP)**3)/6720.0*CL6)*(105.0-105.0*ETA+105.0*ETA**3+
     &      1.0/(-2.0*ETEMP*HH*HH)*(26880.0+4305.0*PI2*ETA-215408.0*ETA
     &      +47040.0*ETA*ETA)-4.0/((-2.0*ETEMP*HH*HH)**2)
     &      *(53760.0-176024.0
     &      *ETA+4305.0*PI2*ETA+15120.0*ETA*ETA))) 

*      IF ((ETEMP+PN1ENERGY(2,IPAIR)/MU).LT.0.) THEN
      POTMIN = -(M12*M12/(2.*L2))
      POTMINPN = -(M12*M12/(2.*L2PN))
      IF ((ETEMP.LT.0.).AND.((ETEMP.GE.POTMIN)
     &     .OR.(ETEMP.GE.POTMINPN))) THEN ! adding restriction over ETEMP after AND
         SEMI = -M12/(2.*ETEMP)
         SEMIPN = -M12/(2.*(ETEMP-PN1ENERGY(2,IPAIR)/MU))
         ECC2 = 1.+(2.*ETEMP*L2)/(M12*M12)
         ECC = SQRT(ECC2)
         ECC2PN = 1.+(2.*(ETEMP+PN1ENERGY(2,IPAIR)/MU)*L2PN)/
     &        (M12*M12)
         ECC2PN = MAX(ECC2PN,1D-08)
         ECCPN = SQRT(ECC2PN)
         ECC2_FULL = MAX(ECC2_FULL,1D-08)
         ECC2_FULL_PN = SQRT(ECC2_FULL)

*     Einstein shift
         A_EIN = 6.0*PI*M12/(SEMI*CL2*(1-ECC2))
c*     Merging timescale estimation
*         GE = (1.0 - ECC2)**3.5
*     &        /(1.0 + (73.0/24.0 + 37.0*ECC2/96.0)*ECC2)
c         TZ = 5.0/64.0*CLIGHTPN**5*GE*SEMI**4/(M1*M2*M12)
*       Form da/dt & de/dt according to Peters 1964.
c         ADOT =  64.0/5.0*M1*M2*(M1+M2)/
c     &   (CLIGHTPN**5*SEMI**3*(1.0 - ECC2)**3.5)
c         ADOT = ADOT*(1.0 + 73.0/24.0*ECC2 + 37.0/96.0*ECC2**2)
c         EDOT = 304.0/15.0*ECC*M1*M2*(M1 + M2)/(CLIGHTPN**5*SEMI**4)
c         EDOT = EDOT/(1.0 - ECC2)**2.5*(1.0 + 121.0/304.0*ECC2)
c         TGR = SEMI/ADOT
 
*         ICOUNT = ICOUNT+1
*         I1 = 2*IPAIR - 1
*         I = NMAX + IPAIR
*         STEP(I1) = DT
*        IF (A_EIN.GT.1.0D-4.AND.MOD(ICOUNT,100).EQ.0) THEN
*        IF(A_EIN.GT.1.0D-4.AND.STEP(I1).GT.10*STEP(I))THEN
*       COUNTBIN = COUNTBIN + 1
       TIME2 = TIME*10000
       IF(TNEW2.EQ.TIME2)GOTO 1111
       IF(A_EIN.GT.1.0D-4.AND.MOD(TIME2,50).EQ.0)THEN!.AND.(MOD(COUNTBIN,100)
*     &    .EQ.0))THEN
         TNEW2=TIME2
          
*         I1 = 2*IPAIR - 1
*         I = N + IPAIR
         

*        IF(TZ.LE.5.0E-04)THEN
         WRITE (50,*) !TIME2,TNEW2,
     & TIME,TZ,TGR,V1_22/CLIGHTPN,A_EIN, 
*,GI_PN,GI_NOPN,DT,
     &        NP, 
     &          N1,N2,M1,M2,
*SEMIPN,ECCPN
*     &           ,ARAD
*     &        ,ECC2_FULL_PN,
*     &        ,SEMI,ECC,
     &        2*(M1+M2)/CL2,
     &        XI(1:6),VI(1:6),XX12(1:3),VV12(1:3)
*     &        HI*MU,PREVRELENERGY,PN1ENERGY(2,IPAIR),
*     &        DIFFS(4),DIFFS(5)
*     &        SPN(1,SPNI1),SPN(2,SPNI1),
*     &        SPN(3,SPNI1), SPN(1,SPNI2),SPN(2,SPNI2),SPN(3,SPNI2),
*     &        SPNI1,SPNI2
 
        ENDIF
        TIME2 = TIME2+1
1111   END IF

      call flush(50)

      RETURN

      END
      
