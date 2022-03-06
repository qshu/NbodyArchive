      SUBROUTINE GIANT2(L,I,W,Q,WSCALE,QSCALE,ZN,QL)
*
*
*       Structure constants of giant star (chain version).
*       --------------------------------------------------
*
*       Theory of Rosemary Mardling, Ap. J. XX, YYY, 1995.
*       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               RP(NTMAX),ES(NTMAX),CM(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      REAL*8  WW(6),QQ(6),W(2),Q(2),WSCALE(2),QSCALE(2),SW(2)
      DATA WW  /2.119d0,3.113d0,8.175d0,3.742d0,4.953d0,9.413d0/
      DATA QQ  /0.4909d0,0.4219d0,0.2372d0,0.4677d0,0.3560d0,0.1519d0/
      DATA A0,A1,A2,A3 /0.944525d0,-0.392030d0,6.01655D-02,-3.34790D-03/
      DATA B0,B1,B2  /-0.3789d0,1.481d0,-0.1018d0/
      DATA C0,C1,C2,C3  /1.455775d0,3.691069d0,-10.42117d0,11.23818d0/
      DATA E0,E1,E2,E3  /1.934977d0,2.214222d0,-4.855796d0,4.025394d0/
*
*
*       Set typical core mass of 0.3/0.5 Msun for chain binary in next bin.
      IC = NCHAOS + 1
      CM(L,IC) = (0.3d0 + 0.1d0*(ISTAR(I) - 3.d0))/ZMBAR
*
*       Form ratio of core mass and mass.
      SIG = CM(L,IC)/M(I)
*
*       Include safety check on mass ratio just in case.
      IF (SIG.GT.0.9d0.OR.SIG.LT.0.d0) THEN
          if(rank.eq.0)
     &    WRITE (6,5)  IC, ISTAR(I), I, M(I)*ZMBAR, CM(L,IC)*ZMBAR
    5     FORMAT (' WARNING!    GIANT2    IC K* I M MC ',2I4,I6,2F7.2)
          SIG = 0.9d0
          CM(L,IC) = 0.9d0*M(I)
      END IF
*
*       Define mass, core mass, radius, envelope mass and luminosity in S.U.
      ZM = M(I)*ZMBAR
      ZCMC = CM(L,IC)*ZMBAR
      RSI = SIZE(I)*SU
      ZME = ZM - ZCMC
*       Obtain L^{1/3} from giant relation L = 1.98D+05*M_c^6.
      ZL3 = 58.3d0*ZCMC**2
*       Evaluate damping constant from Zahn's theory (R.M. 13/5/97).
*       FAC = (GM)^{1/2}*M_{env}^{1/3)/(R^{5/6}*L^{1/3}) = 8.48D+03 for S.U.
      QL = 8.48D+03*SQRT(ZM)*(ZME/ZL3)**0.33d0/RSI**0.833d0
*
*       Set effective frequencies, overlap integrals and structure constants.
      DO 10 K = 1,2
          K1 = 3*K - 2
          IF (K.EQ.1) THEN
              SW(K) = ((C3*SIG + C2)*SIG + C1)*SIG + C0
          ELSE
              SW(K) = ((E3*SIG + E2)*SIG + E1)*SIG + E0
           END IF
           W(K) = SW(K)**2
          Q(K) = ((A3*SW(K) + A2)*SW(K) + A1)*SW(K) + A0
          WSCALE(K) = SQRT(W(K)/WW(K1))
          QSCALE(K) = (Q(K)/QQ(K1)/WSCALE(K))**2
          QSCALE(K) = MAX(QSCALE(K),0.0001D0)
   10 CONTINUE
*
*       Evaluate new polytropic index.
      ZN = (B2*SW(1) + B1)*SW(1) + B0
*     if(rank.eq.0)
*    &WRITE (24,20)  IC, IPAIR, ISTAR(J), CM(L,IC)/M(I), ZN
*  20 FORMAT (' GIANT:    IC KS K* MC/M R/R0 n ',3I4,2F6.2)
*     CALL FLUSH(24)
*
*       Include warning if n > 5.
*     IF (rank.eq.0.and.ZN.GE.5.0) THEN
*         WRITE (6,30)  IC, ISTAR(J), CM(L,IC)/M(I),
*    &                  SIZE(I)/RIN, ZN
*  30     FORMAT (' GIANT:    WARNING!    IC K* MC/M R/R0 n ',
*    &                                    2I4,F6.2,F6.1,F6.2)
*     END IF

      RETURN
*
      END
