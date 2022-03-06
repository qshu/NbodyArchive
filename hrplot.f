      SUBROUTINE HRPLOT
*
*
*       HR diagnostics of evolving stars.
*       ---------------------------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/RCHE/ CMRCH(13,MMAX),NAMER(2,MMAX),KSTARR(3,MMAX)
      REAL*8 TSCLS(20),LUMS(10),GB(10),TM,TN,TDISP
      REAL*8 M0,M01,M1,M2,RM,AGE,LUM,MC,RCC,OSPIN,OSPIN2
      REAL*8 MENV,RENV,K2,K3
      PARAMETER(K3=0.21d0)
      REAL*8 RDD(NMAX),XDD(NMAX),MTOT
      REAL*4 XS(3),VS(3),R1,ZL1,R2,ZL2
      REAL*4 XJ(3,6),VJ(3,6),BODYJ(6)
      INTEGER ICHECK(NMAX),MLIST(NMAX),NCNT,NJ,INDEXJ(KMAX)
      LOGICAL EVOUT,FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
*
*
      EVOUT = .FALSE.
      TPHYS = (TPLOT+TOFF)*TSTAR
      IF(TPLOT.GT.0.0)THEN
         T6 = TPHYS - EPOCH0
         CALL MTURN(TURN,T6,ZPARS)
      ENDIF
      TDISP = TOFF*TSTAR
*
      WRITE(82,100)NPAIRS,TPHYS
      NS = N - 2*NPAIRS
      IMERGE = 0
      NJ = 0
      NCNT = 0
      WRITE(83,101)NS,TPHYS,TCR,TSCALE
      WRITE(83,102)NC,RC,RBAR,RTIDE,(RDENS(K),K=1,3)
      WRITE(83,103)ZMBAR,TURN,RSCALE
 100  FORMAT(I8,F9.1)
 101  FORMAT(I8,F9.1,2F6.2)
 102  FORMAT(I8,3F6.2,3F10.5)
 103  FORMAT(F10.4,2F6.2)
      IF(EVOUT)THEN
         WRITE(93,*)TPHYS
         WRITE(92,*)TPHYS
      ENDIF
*
      DO 15 I = 1,N
         ICHECK(I) = 0
 15   CONTINUE
*
*       Write formatted data bank on unit 99.
      IF(KZ(25).GE.2)THEN
         IF(FIRST)THEN
            OPEN(UNIT=98,STATUS='NEW',FORM='FORMATTED',FILE='BSS')
            FIRST = .FALSE.
         ENDIF
         IBS = 0
         II = 0
         MTOT = 0.D0
      ENDIF
*
      IF(NRSAVE.GT.0) WRITE(6,*)' HRPLOT NRSAVE ',NRSAVE
      DO 20 I = 1,N
         M1 = BODY(I)*ZMBAR
         M2 = 0.D0
         J = 0
         DO 22 K = 1,NRSAVE
            IF(NAMER(1,K).EQ.NAME(I)) GOTO 20
            IF(NAMER(2,K).EQ.NAME(I)) GOTO 20
 22      CONTINUE
*       Binary Stars.
         IF(I.LT.IFIRST)THEN
            JPAIR = KVEC(I)
            J = 2*JPAIR
            IF(I.EQ.J) GOTO 20
            ICM = N + JPAIR
*       Determine merger & ghost index for negative c.m. name.
            IF(NAME(ICM).LT.0)THEN
               CALL FINDJ(I,J,IMERGE)
               IF(NAME(J).GT.NZERO) GOTO 20
               M1 = CM(1,IMERGE)*ZMBAR
               M2 = CM(2,IMERGE)*ZMBAR
               HJ = HM(IMERGE)
               RJ = SQRT(XREL(1,IMERGE)**2 + XREL(2,IMERGE)**2
     &                                     + XREL(3,IMERGE)**2)
               KW = KSTARM(IMERGE)
               IF(NAME(2*JPAIR).GT.0.AND.NAME(2*JPAIR).LE.NZERO)THEN
                  NJ = NJ + 1
                  INDEXJ(NJ) = 2*JPAIR
               ENDIF
            ELSE
*       Check for ghost binary.
               IF(M1.EQ.0.0)THEN
                  IM = 0
*       Search merger table to identify corresponding index of c.m. name.
                  DO 24 K = 1,NMERGE
                     IF(NAMEM(K).EQ.NAME(ICM))THEN
                        IM = K
                     ENDIF
 24               CONTINUE
                  IF(IM.EQ.0) GOTO 20
*       Copy masses for first component.
                  M1 = CM(3,IM)*ZMBAR
                  M2 = CM(4,IM)*ZMBAR
               ELSE
                  M2 = BODY(J)*ZMBAR
               ENDIF
               RJ = R(JPAIR)
               HJ = H(JPAIR)
            ENDIF
*
            IF(M2.LT.0.0001D0)THEN
               WRITE(6,*)' HRPLOT SMALL M2 ',M2
               WRITE(6,*)' JP ICM J1 J2 ',JPAIR,ICM,I,J
               WRITE(6,*)' NM2 KW2 M0 AGE ',NAME(J),KW2,M0,AGE
               M2 = 0.D0
               J = 0
               IF(NAME(I).EQ.0) GOTO 20
            ELSE
               BODYI = (M1 + M2)/ZMBAR
               SEMI = -0.5D0*BODYI/HJ
               ECC2 = (1.D0 - RJ/SEMI)**2
               KW = KSTAR(ICM)
            ENDIF
*
            DO 26 , K = 1,3
               XS(K) = X(K,ICM)
               VS(K) = XDOT(K,ICM)
 26         CONTINUE
         ELSE
*       Replace ghost mass of single star with original value from merged KS.
            IF(M1.EQ.0.0)THEN
               IM = 0
*       Search merger table to identify corresponding index.
               DO 28 , K = 1,NMERGE
                  IF(NAMEM(K).EQ.NAME(I))THEN
                     IM = K
                  ENDIF
 28            CONTINUE
               IF(IM.EQ.0) GOTO 20
               M1 = CM(2,IM)*ZMBAR
            ENDIF
            IF(ICHECK(I).GT.0) GOTO 20
            IF(NAME(I).GT.2*NBIN0) GOTO 50
            DO 40 JJ = I+1,N
               IF(IABS(NAME(I) - NAME(JJ)).EQ.1)THEN
*       Original binaries that are not regularized.
                  J = JJ
                  RIJ2 = 0.0
                  VIJ2 = 0.0
                  RDOT = 0.0
                  DO 42 K = 1,3
                     RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
                     VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
                     RDOT = RDOT + 
     &                      (X(K,I) - X(K,J))*(XDOT(K,I) - XDOT(K,J))
 42               CONTINUE
                  BODYI = BODY(I) + BODY(J)
                  SEMI = 2.0/SQRT(RIJ2) - VIJ2/BODYI
                  IF(SEMI.LE.0.0) GOTO 50
                  IF(ICHECK(J).GT.0) GOTO 50
                  SEMI = 1.D0/SEMI
                  M2 = BODY(J)*ZMBAR
                  ECC2 = (1.D0 - SQRT(RIJ2)/SEMI)**2 + 
     &                    RDOT**2/(SEMI*BODYI)
                  DO 44 , K = 1,3
                     XS(K) = (BODY(I)*X(K,I) + BODY(J)*X(K,J))/BODYI
                     VS(K) = (BODY(I)*XDOT(K,I)+BODY(J)*XDOT(K,J))/BODYI
 44               CONTINUE
                  KW = 0
                  GOTO 60
               ENDIF
 40         CONTINUE
*       Single stars (skip rare chain subsystem).
 50         CONTINUE
            IF(NAME(I).EQ.0) GOTO 20
            DO 52 , K = 1,3
               XS(K) = X(K,I)
               VS(K) = XDOT(K,I)
 52         CONTINUE
         ENDIF
 60      CONTINUE
         IF(M1.LT.1.0D-12) GOTO 20
         M0 = BODY0(I)*ZMBAR
         MC = 0.D0
         KW1 = KSTAR(I)
         CALL star(KW1,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Obtain stellar parameters at current epoch.
         AGE = MAX(TPLOT,TEV0(I))*TSTAR - EPOCH(I)
         
         CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM,LUM,KW1,MC,RCC,MENV,RENV,K2)
         R1 = LOG10(RM)
         ZL1 = LOG10(LUM)
         IF(M1.LT.1.0D-12)THEN
            WRITE(6,*)' HRPLOT M1 BLOW ',i,name(i),kw1,m1
            TEV(I) = TIME + 0.1D0/TSTAR
            IF(I.LT.IFIRST)THEN
               JPAIR = KVEC(I)
               ICM = N + JPAIR
               TEV(ICM) = TEV(I) + 0.1D0/TSTAR
            ENDIF
            GOTO 20
         ENDIF
         OSPIN = MAX(SPIN(I)*SPNFAC,1.0D-10)/
     &           (K2*RM*RM*(M1-MC)+K3*RCC*RCC*MC)
         IF(KZ(25).GE.2)THEN
            IF(KW1.LE.1.AND.(M1.GT.1.02*TURN.OR.TM.LT.TPHYS))THEN
               IF(IBS.EQ.0)THEN
                  WRITE(98,111)TPHYS,N,NBS,RC,TURN
                  IBS = 1
               ENDIF
               WRITE(98,112)NAME(I),M0,M1,EPOCH(I)+TDISP,
     &                      (XS(K),K=1,3)
            ENDIF
         ENDIF
         ICHECK(I) = 1
         IF(M2.EQ.0.0)THEN
*       Single star output.
            WRITE(83,121)NAME(I),KW1,M1,ZL1,R1,(XS(K),K=1,3),
     &                   (VS(K),K=1,3),EPOCH(I)+TDISP,OSPIN
            IF(EVOUT)THEN
               WRITE(93,1121)KW1,M0,M1,(XS(K),K=1,3),
     &                       (VS(K),K=1,3),EPOCH(I)+TDISP,OSPIN
            ENDIF
            NCNT = NCNT + 1
         ELSE
*       Binary star output.
            IF(J.LE.0)THEN
               WRITE(6,*)' HRPLOT J <= 0 ',I,J,NAME(I),NAME(J),M1,M2
               GOTO 20
            ENDIF
            ICHECK(J) = 1
            M01 = M0
            M0 = BODY0(J)*ZMBAR
            MC = 0.D0
            KW2 = KSTAR(J)
            CALL star(KW2,M0,M2,TM,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Obtain stellar parameters at current epoch.
            AGE = MAX(TPLOT,TEV0(J))*TSTAR - EPOCH(J)
            CALL hrdiag(M0,AGE,M2,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                  RM,LUM,KW2,MC,RCC,MENV,RENV,K2)
            R2 = LOG10(RM)
            ZL2 = LOG10(LUM)
            OSPIN2 = MAX(SPIN(J)*SPNFAC,1.0D-10)/
     &               (K2*RM*RM*(M1-MC)+K3*RCC*RCC*MC)
            IF(KZ(25).GE.2)THEN
               IF(KW2.LE.1.AND.(M2.GT.1.02*TURN.OR.TM.LT.TPHYS))THEN
                  IF(IBS.EQ.0)THEN
                     WRITE(98,111)TPHYS,N,NBS,RC,TURN
                     IBS = 1
                  ENDIF
                  WRITE(98,112)NAME(J),M0,M2,EPOCH(J)+TDISP,
     &                         (XS(K),K=1,3)
               ENDIF
            ENDIF
            ECC = SQRT(ECC2)
            PB = DAYS*SEMI*SQRT(ABS(SEMI)/BODYI)
            PB = MIN(PB,999999.9D0)
            PB = LOG10(ABS(PB))
            SEMI = LOG10(ABS(SEMI*SU))
            WRITE(82,123)NAME(I),NAME(J),KW1,KW2,KW,
     &                   ECC,PB,SEMI,M1,M2,ZL1,ZL2,R1,R2,
     &                   (XS(K),K=1,3),(VS(K),K=1,3),
     &                   EPOCH(I)+TDISP,EPOCH(J)+TDISP,OSPIN,OSPIN2
            IF(EVOUT)THEN
               WRITE(92,1123)KW1,KW2,ECC,SEMI,M01,M0,M1,M2,
     &                       (XS(K),K=1,3),(VS(K),K=1,3),
     &                       EPOCH(I)+TDISP,EPOCH(J)+TDISP,OSPIN,OSPIN2
            ENDIF
            NCNT = NCNT + 2
         ENDIF
         IF(KZ(25).GE.3)THEN
            II = II + 1
            RI = XS(1)**2 + XS(2)**2 + XS(3)**2
            RDD(II) = SQRT(RI)
            XDD(II) = M1 + M2
            MLIST(II) = II
            MTOT = MTOT + XDD(II)
         ENDIF
 20   CONTINUE
*
* Outer ghosts. 
*
      DO 65 L = 1,NJ
         J = INDEXJ(L)
         M1 = BODY(J)*ZMBAR
         M0 = BODY0(J)*ZMBAR
         MC = 0.D0
         KW1 = KSTAR(J)
         CALL star(KW1,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Obtain stellar parameters at current epoch.
         AGE = MAX(TPLOT,TEV0(J))*TSTAR - EPOCH(J)
         CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM,LUM,KW1,MC,RCC,MENV,RENV,K2)
         R1 = LOG10(RM)
         ZL1 = LOG10(LUM)
         OSPIN = MAX(SPIN(J)*SPNFAC,1.0D-10)/
     &           (K2*RM*RM*(M1-MC)+K3*RCC*RCC*MC)
         DO 66 K = 1,3
            XS(K) = X(K,J)
            VS(K) = XDOT(K,J)
 66      CONTINUE
         WRITE(83,121)NAME(J),KW1,M1,ZL1,R1,(XS(K),K=1,3),(VS(K),K=1,3),
     &                EPOCH(J)+TDISP,OSPIN
         IF(EVOUT)THEN
            WRITE(93,1121)KW1,M0,M1,(XS(K),K=1,3),(VS(K),K=1,3),
     &                    EPOCH(J)+TDISP,OSPIN
         ENDIF
         NCNT = NCNT + 1
 65   CONTINUE
*
* Chain regularization (case NAME(ICM) = 0).
*
      IF(NCH.GT.0)THEN
         CALL CHDATA(XJ,VJ,BODYJ)
         DO 67 L = 1,NCH
*       Copy global address from common JLIST (set in CHDATA).
            J = JLIST(L)
            M1 = BODYJ(L)*ZMBAR
            M0 = BODY0(J)*ZMBAR
            MC = 0.D0
            KW1 = KSTAR(J)
            CALL star(KW1,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Obtain stellar parameters at current epoch.
            AGE = MAX(TPLOT,TEV0(J))*TSTAR - EPOCH(J)
            CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                  RM,LUM,KW1,MC,RCC,MENV,RENV,K2)
            R1 = LOG10(RM)
            ZL1 = LOG10(LUM)
            OSPIN = MAX(SPIN(J)*SPNFAC,1.0D-10)/
     &              (K2*RM*RM*(M1-MC)+K3*RCC*RCC*MC)
            DO 68 K = 1,3
               XS(K) = XJ(K,L)
               VS(K) = VJ(K,L)
 68         CONTINUE
            WRITE(83,121)NAME(J),KW1,M1,ZL1,R1,(XS(K),K=1,3),
     &                   (VS(K),K=1,3),EPOCH(J)+TDISP,OSPIN
            IF(EVOUT)THEN
               WRITE(93,1121)KW1,M0,M1,(XS(K),K=1,3),(VS(K),K=1,3),
     &                       EPOCH(J)+TDISP,OSPIN
            ENDIF
            NCNT = NCNT + 1
 67      CONTINUE
      ENDIF
*
* Saved Roche Binaries.
*
      DO 70 , I = 1,NRSAVE
         M0 = CMRCH(8,I)*ZMBAR
         M1 = CMRCH(9,I)*ZMBAR
         MC = 0.D0
         KW1 = KSTARR(2,I)
         CALL star(KW1,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
         AGE = TPLOT*TSTAR - CMRCH(10,I)
         CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM,LUM,KW1,MC,RCC,MENV,RENV,K2)
         R1 = LOG10(RM)
         ZL1 = LOG10(LUM)
         DO 72 , K = 1,3
            XS(K) = CMRCH(K+1,I)
            VS(K) = CMRCH(K+4,I)
 72      CONTINUE
         IF(KZ(25).GE.2)THEN
            IF(KW1.LE.1.AND.(M1.GT.1.02*TURN.OR.TM.LT.TPHYS))THEN
               IF(IBS.EQ.0)THEN
                  WRITE(98,111)TPHYS,N,NBS,RC,TURN
                  IBS = 1
               ENDIF
               WRITE(98,112)NAMER(1,I),M0,M1,CMRCH(10,I),(XS(K),K=1,3)
            ENDIF
         ENDIF
         M01 = M0
         M0 = CMRCH(11,I)*ZMBAR
         M2 = CMRCH(12,I)*ZMBAR
         MC = 0.D0
         KW2 = KSTARR(3,I)
         CALL star(KW2,M0,M2,TM,TN,TSCLS,LUMS,GB,ZPARS)
         AGE = TPLOT*TSTAR - CMRCH(13,I)
         CALL hrdiag(M0,AGE,M2,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM,LUM,KW2,MC,RCC,MENV,RENV,K2)
         R2 = LOG10(RM)
         ZL2 = LOG10(LUM)
         IF(KZ(25).GE.2)THEN
            IF(KW2.LE.1.AND.(M2.GT.1.02*TURN.OR.TM.LT.TPHYS))THEN
               IF(IBS.EQ.0)THEN
                  WRITE(98,111)TPHYS,N,NBS,RC,TURN
                  IBS = 1
               ENDIF
               WRITE(98,112)NAMER(2,I),M0,M2,CMRCH(13,I),(XS(K),K=1,3)
            ENDIF
         ENDIF
         BODYI = CMRCH(9,I) + CMRCH(12,I)
         SEMI = CMRCH(1,I)
         ECC = 0.0
         PB = DAYS*SEMI*SQRT(ABS(SEMI)/BODYI)
         OSPIN = 365.25*TWOPI/PB
         PB = MIN(PB,999999.9D0)
         PB = LOG10(ABS(PB))
         SEMI = LOG10(ABS(SEMI*SU))
         WRITE(82,123)NAMER(1,I),NAMER(2,I),KW1,KW2,KSTARR(1,I),
     &                ECC,PB,SEMI,M1,M2,ZL1,ZL2,R1,R2,
     &                (XS(K),K=1,3),(VS(K),K=1,3),
     &                CMRCH(10,I),CMRCH(13,I),OSPIN,OSPIN
         IF(EVOUT)THEN
            WRITE(92,1123)KW1,KW2,ECC,SEMI,M01,M0,M1,M2,
     &                    (XS(K),K=1,3),(VS(K),K=1,3),
     &                    CMRCH(10,I),CMRCH(13,I),OSPIN,OSPIN
         ENDIF
         NCNT = NCNT + 2
         IF(KZ(25).GE.3)THEN
            II = II + 1
            RI = XS(1)**2 + XS(2)**2 + XS(3)**2
            RDD(II) = SQRT(RI)
            XDD(II) = M1 + M2
            MLIST(II) = II
            MTOT = MTOT + XDD(II)
         ENDIF
 70   CONTINUE
      NRSAVE = 0
      WRITE(6,*)' HRPLOT CHECK: N NCNT ',N,NCNT
*
      KW1 = -1000
      SEMI = 0.0
      WRITE(82,123)KW1,KW1,KW,KW2,KW,
     &             ECC,PB,SEMI,M1,M2,ZL1,ZL2,R1,R2,
     &             (XS(K),K=1,3),(VS(K),K=1,3),SEMI,SEMI,OSPIN,OSPIN
      WRITE(83,121)KW1,KW,M1,ZL1,R1,(XS(K),K=1,3),(VS(K),K=1,3),
     &             SEMI,OSPIN
 111  FORMAT(F9.1,I8,I4,2F6.2,' HRPLOT')
 112  FORMAT(I6,2F7.3,F9.3,3F10.5)
 121  FORMAT(I6,I3,3F7.3,6F10.5,F10.3,1P,E12.4)
 1121 FORMAT(I3,2F7.3,6F10.5,F10.3,1P,E12.4)
 123  FORMAT(2I6,2I3,I4,F6.3,F10.2,7F7.3,6F10.5,2F10.3,1P,2E12.4)
 1123 FORMAT(2I3,F6.3,5F7.3,6F10.5,2F10.3,1P,2E12.4)
*
      IF(KZ(25).GE.3)THEN
*
*       Sort radii into increasing order.
         CALL SORT1(II,RDD,MLIST)
*
*       Determine the half-mass radius (pc) and relaxation time (Myrs).
         MC = 0.D0
         DO 80 , I = 1,II
            MC = MC + XDD(MLIST(I))
            IF(MC.GE.0.5*MTOT) GOTO 85
 80      CONTINUE
 85      CONTINUE
         RHALF = RDD(I)*RBAR
         TRH = FLOAT(II)/LOG10(0.4D0*FLOAT(II))
         TRH = 0.858D0*TRH*SQRT(RHALF**3/MTOT)
         WRITE(6,125)RHALF,TRH
 125     FORMAT('  HALF-MASS RH TRH ',F6.2,F12.4)
      ENDIF
*       Update next output time.
*     IF(TIME.GT.0.0) DTPLOT = SQRT(2.D0)*DTPLOT
      TPLOT = TPLOT + DTPLOT
*
      CALL FLUSH(82)
      CALL FLUSH(83)
      TPHYS = (TIME+TOFF)*TSTAR
      IF(KZ(25).GE.2.AND.IBS.GT.0) CALL FLUSH(98)
*
      RETURN
      END
***
