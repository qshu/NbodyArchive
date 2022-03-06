      SUBROUTINE ESCAPE
*
*       Escaper detection.
*       ------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX)
      LOGICAL  FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
*
*
*       Adopt twice the tidal radius as escape condition.
      RESC2 = 4.0*RTIDE**2
*       For tidal cutoff check only energy
      IF(KZ(23).GE.3)THEN
          RESC2 = 0.D0
          ETID = ZMASS/RTIDE
          ZMOLD = ZMASS
      END IF
*
      RTIDE2 = RTIDE**2
      NCORR = 0
      NCRIT1 = 0
      NCRIT2 = 0
      DO 1 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
    1 CONTINUE
*
      I = IFIRST
*       Set the distance (squared) with respect to the density centre.
    5 RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                               (X(3,I) - RDENS(3))**2
      IF (RI2.LT.RTIDE2) NCRIT1 = NCRIT1 + 1
*
      IF(KZ(23).LE.2) THEN
*       See whether escape is indicated (retain merger ghost particles).
      IF (RI2.GT.RESC2.AND.RI2.LT.1.0D+10) THEN
          RJMIN2 = 1.0D+10
*       Find distance to the nearest neighbour and calculate potential.
          POTI = 0.0D0
          DO 8 J = IFIRST,NTOT
              IF (J.EQ.I) GO TO 8
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
              POTI = POTI + BODY(J)/SQRT(RIJ2)
              IF (RIJ2.LT.RJMIN2) THEN
                  RJMIN2 = RIJ2
                  JMIN = J
              END IF
    8     CONTINUE
          VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
*       Check escape criterion for tidal case or isolated system.
          EI = 0.5*VI2 - POTI
          IF (KZ(14).GT.0.OR.EI.GT.0.0) GO TO 30
      END IF
*
      ELSE
*       Check tidal cutoff criterion - use last computed phi-value
*       Find distance to the nearest neighbour
          POTI = 0.0D0
          NNB = LIST(1,I)
          RJMIN2 = 1.0D+10
          DO 9 L = 1,NNB
              J = LIST(L+1,I)
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
              IF (RIJ2.LT.RJMIN2) THEN
                  RJMIN2 = RIJ2
                  JMIN = J
              END IF
    9     CONTINUE
          VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
          POTI = -PHIDBL(I)
*       Check escape criterion for tidal case or isolated system.
          EI = 0.5*VI2 - POTI + ETID
*
          IF (EI.GT.0.0.AND.RI2.GT.RTIDE2) THEN
          NCRIT2 = NCRIT2 + 1
          GO TO 30
          END IF
      END IF
*
   10 I = I + 1
   12 IF (I.LE.NTOT) GO TO 5
*
      IF (NCORR.EQ.0) GO TO 25
*
*       Form centre of mass terms.
      DO 15 I = 1,N
          DO 14 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)/ZMASS
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)/ZMASS
   14     CONTINUE
   15 CONTINUE
      CMR(4) = SQRT(CMR(1)**2 + CMR(2)**2 + CMR(3)**2)
*
      if(rank.eq.0)then
      WRITE (6,18)  N, NSESC, NBESC, ZMASS, BE(3), CMR(4), RESC, STEPI,
     &              RSI, ZMASS/FLOAT(N), NCRIT1, (JLIST(J),J=1,NCORR)
   18 FORMAT (/,' ESCAPE   N =',2I5,I4,F8.4,F11.6,2F7.2,F7.3,F7.2,F9.5,
     &                                        I6,2X,8I5,/,5(16X,23I5,/))
*
      IF (KZ(23).EQ.2.OR.KZ(23).EQ.4)THEN
          JLAST = 2*NCORR
          WRITE (6,20)  (ILIST(J),J=1,JLAST)
   20     FORMAT (/,' ESCAPE ANGLES ',11(2I4,2X),9(/,15X,11(2I4,2X)))
      END IF
      end if
*
*       Check updating of global index for chain c.m.
      IF (NCH.GT.0) THEN
          CALL CHFIND
      END IF
*
*       Set phase indicator = -1 for new NLIST in routine INTGRT.
      IPHASE = -1
*
   25 CONTINUE
*
      RETURN
*
   30 A2 = (X(1,JMIN) - RDENS(1))**2 + (X(2,JMIN) - RDENS(2))**2 +
     &                                 (X(3,JMIN) - RDENS(3))**2
*       See whether nearest body satisfies escape condition or RIJ > 10*RMIN.
      IF (A2.GT.RESC2.AND.KZ(14).GT.0) GO TO 40
      IF (RJMIN2.GT.100.0*RMIN2) GO TO 40
*
      IF(KZ(23).GE.3)THEN
          VI2 = XDOT(1,JMIN)**2 + XDOT(2,JMIN)**2 + XDOT(3,JMIN)**2
          EI = 0.5*VI2 + PHIDBL(JMIN) + ETID
          IF(EI.GT.0.D0) GO TO 40
      END IF
*
      A3 = XDOT(1,JMIN)**2 + XDOT(2,JMIN)**2 + XDOT(3,JMIN)**2
      A4 = (XDOT(1,I) - XDOT(1,JMIN))**2 +
     &     (XDOT(2,I) - XDOT(2,JMIN))**2 +
     &     (XDOT(3,I) - XDOT(3,JMIN))**2
      A5 = (BODY(I) + BODY(JMIN))/SQRT(RJMIN2)
*       Check velocity of binary component in case of bound pair.
      A6 = 2.0*A5/A4
      IF (A6.GT.1.0.AND.A3.GT.2.0*ZMASS/SQRT(A2)) GO TO 40
*       Retain escaper if dynamical effect on neighbour is significant.
      IF (A6.GT.0.01) GO TO 10
*
*       Form optional output diagnostics.
   40 X1 = X(1,I) - RDENS(1)
      Y1 = X(2,I) - RDENS(2)
      Z1 = X(3,I) - RDENS(3)
      A3 = ABS(Y1/X1)
      A4 = ABS(Z1)/SQRT(X1**2 + Y1**2)
      ILIST(2*NCORR+1) = DATAN(A3)*180.0/3.14159
      ILIST(2*NCORR+2) = DATAN(A4)*180.0/3.14159
*       Escape angles with respect to the X-axis and XY-plane.
      RESC = SQRT(RI2)
      STEPI = STEP(I)
      RSI = RS(I)
*
*       Accumulate escaper names and save current name in case I > N.
      NCORR = NCORR + 1
      JLIST(NCORR) = NAME(I)
      NAMEI = NAME(I)
      KSTARI = KSTAR(I)
      IF (NAME(I).LT.0) KSTARI = KSTARI + 20
      IF (BODY(I).LE.0.0D0) GO TO 50
*
*       Obtain binding energy of body #I and update optional tidal radius.
      ZK = 0.5D0*BODY(I)*VI2
      IF (KZ(14).GT.0) THEN
          CALL XTRNLV(I,I)
          ZK = ZK + HT
          RTIDE = (ZMASS/TIDAL(1))**0.3333
      END IF
      EI = ZK - BODY(I)*POTI
*
*       Update tidal radius in case of tidal cutoff (R.Sp.)
      IF(KZ(23).GE.3)THEN
          RTOLD = RTIDE
          RTIDE = RTIDE*(ZMASS/ZMOLD)**0.3333
      END IF
*
*       Correct total energy.
      BE(3) = BE(3) - EI
*
*       Update total mass and save energy & number of single escapers.
      IF (I.LE.N) THEN
          ZMASS = ZMASS - BODY(I)
          E(4) = E(4) + EI
          NPOP(4) = NPOP(4) + 1
          NSESC = NSESC + 1
      END IF
*
*       Include optional escape output on unit 11.
      if(rank.eq.0)then
      IF (KZ(23).EQ.2.OR.KZ(23).EQ.4) THEN
          IF (FIRST) THEN
              OPEN (UNIT=11,STATUS='NEW',FORM='FORMATTED',FILE='ESC')
              FIRST = .FALSE.
          END IF
          TESC = TSCALE*TIME
          EESC = EI/BODY(I)
          VB2 = 2.0*ZKIN/ZMASS
*       Distinguish between tidal field and isolated system (ECRIT at RTIDE).
          IF (KZ(14).GT.0) THEN
              ECRIT = -1.5*(TIDAL(1)*ZMASS**2)**0.3333
              EESC = 2.0*(EESC - ECRIT)/VB2
              BESC = ZMBAR*BODY(I)
          ELSE
              EESC = 2.0*EESC/VB2
              BESC = BODY(I)/BODYM
          END IF
          VKM = SQRT(VI2)*VSTAR
          WRITE (11,45)  TESC, BESC, EESC, VKM, KSTARI
   45     FORMAT (F8.1,3F6.1,I4)
      END IF
      end if
*
*       Reduce particle number & total membership.
   50 N = N - 1
      NTOT = NTOT - 1
      NNBOPT = MIN((N - NPAIRS)/2,NNBOPT)
      NNBMAX = NNBOPT
      ZNBMAX = 0.9*NNBMAX
      ZNBMIN = MAX(0.2*NNBOPT,1.0)
*
*       Set indicator for removal of c.m. or KS components.
      KCASE = 1
*
*       Update COMMON arrays to remove the escaper and correct F & FDOT.
      CALL REMOVE(I,1)
*
*       Delete escaper from neighbour lists and reduce higher locations.
   60 DO 150 J = 1,NTOT
          NNB = LIST(1,J)
          IF (NNB.EQ.0) GO TO 150
          L = 2
   70     IF (LIST(L,J).NE.I) GO TO 130
*
*       Move up the remaining list members and reduce membership by one.
          DO 80 K = L,NNB
              LIST(K,J) = LIST(K+1,J)
   80     CONTINUE
          LIST(1,J) = LIST(1,J) - 1
*       Reduce the steps to minimize error effect (do not allow DT < 0).
*         STEP(J) = MAX(0.5D0*STEP(J),TIME - T0(J))
*         STEPR(J) = MAX(0.5D0*STEPR(J),TIME - T0R(J))
*       Add body #J to time-step list unless already a member.
*         IF (T0(J) + STEP(J).LT.TLIST) THEN
*             CALL NLMOD(J,1)
*         END IF
          IF (LIST(1,J).GT.0) GO TO 130
*
*       Add a distant body as neighbour if list only contains escaper.
          K = IFIRST - 1
  100     K = K + 1
          RJK2 = (X(1,J) - X(1,K))**2 + (X(2,J) - X(2,K))**2 +
     &                                  (X(3,J) - X(3,K))**2
          IF (RJK2.LT.RESC2.AND.K.LT.N) GO TO 100
          LIST(1,J) = 1
          LIST(2,J) = K
          GO TO 150
*
*       Reduce higher particle locations by one.
  130     IF (LIST(L,J).GT.I) LIST(L,J) = LIST(L,J) - 1
          L = L + 1
          IF (L.LE.LIST(1,J) + 1) GO TO 70
  150 CONTINUE
*
*       Modify time-step list due to escaper removal (-2 for extra test).
      CALL NLMOD(I,-2)
*
*       Update list of old KS components (remove #I and rename > I).
      NNB = LISTR(1)
      DO 170 L = 2,NNB+1
          IF (LISTR(L).EQ.I) THEN
*       Remove both components of pair and reduce membership by two.
              J = 2*KVEC(L-1)
              DO 165 K = J,NNB-1
                  LISTR(K) = LISTR(K+2)
  165         CONTINUE
              LISTR(1) = LISTR(1) - 2
          END IF
  170 CONTINUE
*
*       Reduce higher particle locations by one (separate loop for pairs).
      DO 180 L = 2,NNB+1
          IF (LISTR(L).GT.I) LISTR(L) = LISTR(L) - 1
  180 CONTINUE
*
*       Update list of high velocity particles (remove #I and rename > I).
      NNB = LISTV(1)
      DO 190 L = 2,NNB+1
          IF (LISTV(L).EQ.I) THEN
*       Remove escaper and reduce the membership.
              DO 185 K = L,NNB
                  LISTV(K) = LISTV(K+1)
  185         CONTINUE
              LISTV(1) = LISTV(1) - 1
          END IF
*       Reduce higher particle locations by one (or three for c.m.).
          IF (LISTV(L).GT.I) THEN
              LISTV(L) = LISTV(L) - 1
              IF (I.GT.N+1) LISTV(L) = LISTV(L) - 2
          END IF
  190 CONTINUE
*
      IF (KCASE.GT.1) GO TO 215
*       See whether the escaper is a single particle or c.m.
      IF (I.LE.N + 1) GO TO 12
*
*       Prepare removal of regularized pair.
      IPAIR = I - N - 1
*
*       Skip correction if ghost is also merged binary (NAMEI = 0 below).
      IF (NAMEI.NE.0) THEN
          ZMB = BODY(2*IPAIR-1) + BODY(2*IPAIR)
          EB = H(IPAIR)*BODY(2*IPAIR-1)*BODY(2*IPAIR)/ZMB
          SEMI = -0.5*ZMB/H(IPAIR)
          ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(SEMI*ZMB)
          ECC = SQRT(ECC2)
          NAME1 = NAME(2*IPAIR-1)
          KSTAR1 = KSTAR(2*IPAIR-1)
          PCRIT = R0(IPAIR)
      ELSE
          EB = 0.0
          ECC = 0.0
      END IF
*
*
*       Update the total energy (or ECOLL) and set scaled escape velocity.
*     IF (KZ(27).GT.0.AND.EB.LT.-1.0) THEN
*       Note problem with energy production in BINOUT if EB added here!
*         ECOLL = ECOLL + EB
*     ELSE
          BE(3) = BE(3) - EB
*     END IF
      VI = SQRT(0.5*VI2*ZMASS/ZKIN)
*
*       Check optional diagnostics for hierarchical systems.
      IF (NAMEI.LE.0.AND.KZ(11).GT.0) THEN
          IPHASE = 7
          CALL HIARCH(IPAIR)
      END IF
*
*       Specify binary type (0: standard; -1: merged binary).
      M = 0
      IF (NAMEI.LE.0) M = -1
*
      if(rank.eq.0)then
      WRITE (6,200)  IPAIR, NAME(2*IPAIR-1), NAME(2*IPAIR),
     &               LIST(2,2*IPAIR), M, BODY(2*IPAIR-1)*ZMBAR,
     &               BODY(2*IPAIR)*ZMBAR, EB, VI, ECC, EI
  200 FORMAT (/,' BINARY ESCAPE    KS =',I4,'  NAME =',2I6,I4,I3,
     &                     '  M =',2F5.1,'  EB =',F8.4,
     &                     '  V/<V> =',F5.2,' E =',F5.2,'  EI =',F8.5)
      end if
*
*       Accumulate escaping binary energies and increase the counter.
      IF (LIST(2,2*IPAIR).EQ.-1) THEN
          E(5) = E(5) + EB
          E(6) = E(6) + EI
          NPOP(5) = NPOP(5) + 1
      ELSE
          E(7) = E(7) + EB
          E(8) = E(8) + EI
          NPOP(6) = NPOP(6) + 1
      END IF
*
      IF (M.EQ.0) THEN
          NBESC = NBESC + 1
      ELSE
          NMESC = NMESC + 1
      END IF
*
*       Reduce particle number, pair index & single particle index.
      N = N - 1
      NPAIRS = NPAIRS - 1
      IFIRST = 2*NPAIRS + 1
*
*       Move up all tables of regularized pairs below IPAIR.
      IF (IPAIR.LE.NPAIRS) THEN
          CALL REMOVE(IPAIR,2)
      END IF
*
  215 KCASE = KCASE + 1
      IF (KCASE.LE.3) THEN
*       Remove COMMON arrays of the second component before the first.
          I = 2*IPAIR + 2 - KCASE
          NTOT = NTOT - 1
*       Reduce NTOT by 3 and N by 2 when KS pair escapes.
          CALL REMOVE(I,3)
          GO TO 60
      END IF
*
*       Include the case of escaping merger.
      IF (NAMEI.GE.0) GO TO 250
*
*       Locate current position in the merger table.
      IM = 0
      DO 220 K = 1,NMERGE
          IF (NAMEM(K).EQ.NAMEI) IM = K
  220 CONTINUE
      IF (IM.EQ.0) GO TO 250
*
*       Copy merger energy to respective energy bins (primordial or new).
      ZMU = CM(1,IM)*CM(2,IM)/(CM(1,IM) + CM(2,IM))
      IF (NBIN0.GT.0) THEN
          E(5) = E(5) + ZMU*HM(IM)
      ELSE
          E(7) = E(7) + ZMU*HM(IM)
      END IF
*
*       Reduce current total & merger energy by dominant contribution.
      BE(3) = BE(3) - ZMU*HM(IM)
      EMERGE = EMERGE - ZMU*HM(IM)
*
*       Reduce membership if last or only case (otherwise done in RESET).
      IF (IM.EQ.NMERGE) THEN
          NMERGE = NMERGE - 1
      END IF
*
*       Identify merged ghost particle (single body or binary c.m.).
      JCOMP = -1
      DO 230 J = IFIRST,NTOT
          IF (BODY(J).EQ.0.0D0.AND.NAME(J).EQ.NAMEG(IM)) JCOMP = J
  230 CONTINUE
*
*       Skip if safety procedure fails to identify the correct ghost.
      IF (JCOMP.GT.0) THEN
*       Remove the ghost particle (NAME = 0 & EB = 0 for second binary).
          I = JCOMP
          NAMEI = 0
          NPOP(7) = NPOP(7) + 1
          SEMI = -0.5*(CM(1,IM) + CM(2,IM))/HM(IM)
          if(rank.eq.0)then
          WRITE (6,240)  NAME1, NAMEG(IM), KSTAR1, KSTAR(JCOMP),
     &                   KSTARM(IM), ZMU*HM(IM), PCRIT/SEMI, SEMI
  240     FORMAT (/,' HIARCH ESCAPE    NM =',2I6,'  K* =',3I3,
     &              '  EB =',F8.4,'  PC/A =',F5.1,'  A =',1P,E8.1)
          end if
*
*       Set large look-up times for ghost KS in case of escaping quadruple.
          IF (I.GT.N) THEN
              JP = I - N
              TEV(2*JP-1) = 1.0D+10
              TEV(2*JP) = 1.0D+10
          END IF
          GO TO 50
      END IF
*
  250 I = IPAIR + N
      GO TO 12
*
      END
