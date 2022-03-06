      SUBROUTINE DISRUPT(I,rt)
*
*     Disrupting the bound close star
*     Similar to cmbody.f
*
      INCLUDE 'common6.h'
      COMMON/BHole/   BHx(3), BHv(3),BHTIME
      COMMON/DISR/    IDISR 
      REAL*8 rdot(3), L2, rrel(3), r2, rdot2, mu, EN, BODYold,EN1,
     &       XII(NMAX)
*
      REAL*8 rt, rnorm(3),SMBHOLD
      INTEGER I
*
      IF(I.GT.N) then
         GO TO 20
      END IF
      BODYold=BODY(I)
*     For potential energy correction:
         NNB = LIST(1,NTOT+1)
         CALL XVPRED(NTOT+1,NNB)
*       Replace second component temporarily by arbitrary body.
      DO 12 L = 1,NNB
          JPERT(L) = LIST(L+1,NTOT+1)
          JLIST(L) = LIST(L+1,NTOT+1)
*     changed nbpot direcly!!!(23.07.07)
          if (JPERT(L).eq.I) then
*     PROBLEM when jpert gt n
             JPERT(L) = (NTOT+1) + 1
             if (JPERT(L).eq.I) JPERT(L) = (NTOT+1)+2
             IF (JPERT(L).GE.N) then 
                JPERT(L) = (NTOT+1)-1
                if (JPERT(L).eq.I) JPERT(L) = (NTOT+1)-2
             end if
*     END PROBLEM
          END IF
 12    CONTINUE
*
       CALL NBPOTMBH(NNB,NNB,POT1)
       PRINT*,'POT1======================',POT1,r2
*       print*, 'DISR gets pot1: ', pot1, name(I), name(CURJ),
*     &      nnb,(name(JPERT(m)),
*     &      m=1,NNB) 

*
*     Calculate the reduced mass:
      mu=SMBH*BODYold/(SMBH+BODYold)
*
      rdot(1)=XDOT(1,I)-BHV(1)
      rdot(2)=XDOT(2,I)-BHV(2)
      rdot(3)=XDOT(3,I)-BHV(3)
*
      rdot2=rdot(1)**2+rdot(2)**2+rdot(3)**2
*
      rrel(1)=X(1,I)-BHx(1)
      rrel(2)=X(2,I)-BHx(2)
      rrel(3)=X(3,I)-BHx(3)
*
      r2=rrel(1)**2+rrel(2)**2+rrel(3)**2
*
         DO 9 K=1,3
            rnorm(K)=rrel(K)/sqrt(r2)
 9       CONTINUE
*     Radial velocity:
         vr=rnorm(1)*rDOT(1)+rnorm(2)*rDOT(2)+rnorm(3)*rDOT(3)
*
*     Calculate the kinetik Energy:
*
      EN1=0.5D0*BODY(I)*rdot2
*
      PRINT*,'EN,EN1=============',EN,EN1
*
*     Energy-conservation: 
      E2pot = EN1
*
*     Accrete the star if bound:
*      IF (EN.lt.0.0)then 
*     set BH pos. and vel. to c.m.:
      T0(NTOT+1) = TIME
*
*     accrete the mass:
         SMBH=SMBH+BODYold
*       Create new body from c.m. and initialize zero mass ghost in #I2.
*     add the twobody-potential-energy to ECOLL:
*
*       Ensure that ghost will escape next output (far from fast escapers).
*
      BODY(I) = 0.0
      LIST(1,I) = 0
      T0(I) = TADJ + DTADJ 
*      STEP(I) = 1.0D+06
*      STEPR(I) = 1.0D+06
      RI = SQRT(X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
      VI = SQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
*
      DO 30 K = 1,3
*       Ensure that ghost will escape next output (far from fast escapers).
          X0(K,I) = 1000.0*RSCALE*X(K,I)/RI
          X(K,I) = X0(K,I)
          X0DOT(K,I) = SQRT(0.004*ZMASS/RSCALE)*XDOT(K,I)/VI
          XDOT(K,I) = X0DOT(K,I)
          F(K,I) = 0.0D0
          FDOT(K,I) = 0.0D0
          D0(K,I) = 0.0
          D1(K,I) = 0.0
          D2(K,I) = 0.0D0
          D3(K,I) = 0.0D0
          D0R(K,I) = 0.0
          D1R(K,I) = 0.0
          D2R(K,I) = 0.0D0
          D3R(K,I) = 0.0D0
   30 CONTINUE
*
*
       CALL NBPOTMBH(NNB,NNB,POT2)
       PRINT*,'POT2======================',POT2,E2pot
*       
       DP = POT2 - POT1
       ECOLL = ECOLL + DP
       ECOLL = ECOLL + E2POT
       E(10) = E(10) + E2POT
*        
       IF(NAME(I).EQ.2759) PRINT*,'RI,VI++++++++++++++',RI,VI
*
*       Remove the ghost particle from neighbour lists containing curj.
      JPERT(1) = I
*
      JLIST(1) = I
*     CALL NBREM(I,1,NNB)
*     see nbrem.f (30) for details.
*ChOM      CALL NBREM(curj,1,NNB)      
*
*       Remove ghost from list of curj (use NTOT as dummy here).
      JPERT(1) = NTOT+1
*ChOM      CALL NBREM(NTOT,1,1)

*       Set IPHASE < 0 to ensure updating of time-step sequence.
      IPHASE = -1
      IDISR = 1
         print*, 'DISRUPT disrupts', Name(I), E2POT, DP
        GO TO 100
*
 20     CONTINUE
*      Star to disrupt is a binary:
        PRINT*,' Star to disrupt is a binary:'
        KSPAIR=I-N
        IPHASE = 2
*     Check for rare case of merged binary component.
        IF (NAME(I).LT.0) IPHASE = 7
        PRINT*,'RSCALE==========',RSCALE
        print*, 'DISRUPT terminated a binary', TIME+TOFF, I-N

 100    RETURN
        END
