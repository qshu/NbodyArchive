      SUBROUTINE SUBINT(IQ,I10)
*
*
*       Decision-making for subsystems.
*       -------------------------------
*       R. Spurzem / S.J. Aarseth Experimental Version for parallel binaries
*       Sep 2001/July 2002
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      REAL*8  TSLIST(10*KMAX),TXLIST(10*KMAX)
      INTEGER  BSLIST(10*KMAX)
      SAVE  IRUN,LI
      DATA  IRUN /0/
      COMMON/KSSTAT/ISTEPP,ISTEPU,IBINP,IBINU
      DATA ICALL,KBSUM,ISBSUM,KSBSUM /0,0,0,0/
*ifdef PARALLEL
      integer inum(maxpe),ista(maxpe)
      REAL*8 ZMPI(55,KMAX)
      INTEGER IMPI(KMAX)
*endif
*
*     Determine minimum quantized step near binary steps.
      DTBMIN = 1.D30
      DO 999 I = N+1, NTOT
      IPAIR = I - N
      I1 = 2*IPAIR - 1
      DT = STEP(I1)
      CALL STEPK(DT,DTN)
      IF(DTN.LT.DTBMIN) DTBMIN = DTN
 999  CONTINUE
*
      ISTEPP=0
      ISTEPU=0
      IBINP=0
      IBINU=0
      ICALL = ICALL + 1
      FAC = 1.0/LOG(1.9999999)
      JMM = INT(1 - LOG(DTBMIN)*FAC)
      DTBL = TBLOCK - TPREV
      JMN = INT(1 - LOG(DTBL)*FAC)
      PRINT*,' DTBL,DTBMIN,JMM,JMN=',DTBL,DTBMIN,JMM,JMN
*
      IF(JMM.LT.JMN)THEN
      KBLOCK = 1
      ELSE
      KBLOCK = 2**(JMM-JMN)
      END IF
*
      DTBL = DTBL/KBLOCK
      KBSUM = KBSUM + KBLOCK
      IF(MOD(ICALL,1).EQ.0)PRINT*,' SUBINT: Av KBLOCK=',
     * REAL(KBSUM)/REAL(ICALL)
*
*       Determine correct index after restart (NNTB = 0 initially).
      IF (IRUN.EQ.0) THEN
          IRUN = 1
          TI = 1.0D+10
*       Find smallest sum by looping backwards (avoids multiple entries).
          DO 4 K = NNTB,1,-1
              J = KBLIST(K)
              TJ = T0(J) + STEP(J)
              IF (TJ.LT.TI) THEN
                  TI = TJ
                  LI = K
              ELSE
*       Adopt the previous index (increased by 1 below).
                  LI = LI - 1
                  GO TO 1
              END IF
    4     CONTINUE
      END IF
*
    1 CONTINUE
*       See whether to advance any KS solutions at start of block-step.
      IF (NPAIRS.GT.0) THEN
*       Obtain list of all KS pairs due in interval DTB.
          IF (TBLIST.LE.TBLOCK.OR.NPAIRS.NE.NBPREV.OR.
     *       NNTB.GT.KMAX) THEN
          call cputim(ttxx1)
              IF (DTB.EQ.0.0D0) THEN
                  DTB = MAX(DTMIN,TBLOCK - TPREV)
              END IF
    2         TBLIST = TPREV + DTB
              TBLIST = MAX(TBLOCK,TBLIST)
              NNTB = 0
              DO 3 JPAIR = 1,NPAIRS
                  J1 = 2*JPAIR - 1
                  TXLIST(J1) = T0(J1) + STEP(J1)
                  IF (T0(J1) + STEP(J1).LE.TBLIST) THEN
                      NNTB = NNTB + 1
                      KBLIST(NNTB) = J1
                      TSLIST(NNTB) = T0(J1) + STEP(J1)
                  END IF
    3         CONTINUE
*       Increase interval on zero membership.
              IF (NNTB.EQ.0) THEN
                  DTB = 2.0*DTB
                  GO TO 2
              END IF
*       Stabilize interval on membership of 2*SQRT(NPAIRS).
              NBTRY = 2*SQRT(FLOAT(NPAIRS))
              IF (NNTB.GT.NBTRY)  DTB = 0.75*DTB
              IF (NNTB.LT.NBTRY)  DTB = 1.25*DTB
*       Sort the time-step list sequentially in KBLIST and reset pointer.
              IF (NNTB.GT.1) THEN
                  CALL SORT1(NNTB,TSLIST,KBLIST)
              END IF
*
      J1 = KBLIST(1)
      JNNTB = KBLIST(NNTB)
      PRINT*,' KBL ',J1,JNNTB,T0(J1)+STEP(J1),T0(JNNTB)+STEP(JNNTB)
*
              LI = 0
          call cputim(ttxx2)
          ttkbs = ttkbs + (ttxx2-ttxx1)*60.
          END IF
*
          SPREV = TPREV
*
          PRINT*,' KBLOCK,SPREV,DTBL,NNTB,TBLIST ',KBLOCK,SPREV,
     *  DTBL,NNTB,TBLIST
          J1 = KBLIST(1)
          JNNTB = KBLIST(NNTB)
          PRINT*,' 1-NNTB times ',J1,JNNTB,T0(J1)+STEP(J1),
     *             T0(JNNTB)+STEP(JNNTB)
          CALL FLUSH(6)
          PRINT*,'--------------------'
          DO 5 IX = 1,KBLOCK
*
          SBLOCK = SPREV + DTBL
*       Form list of any KS pairs due in the current sub-block-step.
          KSB = 0
          DO 6 LJ = LI+1,NNTB
              J1 = KBLIST(LJ)
              TJ = T0(J1) + STEP(J1)
              IF (T0(J1) + STEP(J1).LE.SBLOCK) THEN
                  IF(LIST(1,J1).GT.0)THEN
                     IBINP = IBINP + 1
                  ELSE
                     IBINU = IBINU + 1
                  END IF
                  KSB = KSB + 1
                  BSLIST(KSB) = J1
              ELSE
*       Stop searching if sub-block level is reached.
                  GO TO 65
              END IF
    6     CONTINUE
   65     CONTINUE
*
*       Continue if no binaries to integrate
          IF (KSB.EQ.0) GO TO 9
*
          ISBSUM = ISBSUM + 1
          KSBSUM = KSBSUM + KSB
        PRINT*,' Start Do Loop KSB, first min 5=',KSB,
     *  ' BSLIST=',(BSLIST(L),L=1,MIN(5,KSB))
        PRINT*,' Average KSB=',REAL(KSBSUM)/REAL(ISBSUM)
        CALL FLUSH(6)
          PRINT*,'--------------------'
*       Advance binaries due in the current sub-block-step.
          DO 7 L = 1,KSB
*
              I1 = BSLIST(L)
              IPAIR = KVEC(I1)
              TIME = T0(I1) + STEP(I1)
              IMULT = 0
*
   10         CONTINUE
          PRINT*,' L, I1, IPAIR, TIME=',L,I1,IPAIR,TIME,IMULT,STEP(I1)
          CALL FLUSH(6)
*
              CALL KSINT(I1)
*
          PRINT*,' Ret KSINT T0+STEP STEP=',T0(I1)+STEP(I1),STEP(I1)
          CALL FLUSH(6)
*       Check for multiple calls of #I1 (saves using CALL INSERT).
              IF (IPHASE.EQ.0) THEN
                  TI = TIME + STEP(I1)
                  IF (TI.LE.SBLOCK) THEN
                      TIME = TI
                      IMULT = 1
                      GO TO 10
                  END IF
              END IF
*
*         DO 71 K = 1,4
*         ZMPI(K,L) = U(K,IPAIR)
*         ZMPI(K+4,L) = U0(K,IPAIR)
*         ZMPI(K+8,L) = UDOT(K,IPAIR)
*         ZMPI(K+12,L) = FU(K,IPAIR)
*         ZMPI(K+16,L) = FUDOT(K,IPAIR)
*         ZMPI(K+20,L) = FUDOT2(K,IPAIR)
*         ZMPI(K+24,L) = FUDOT3(K,IPAIR)
*         ZMPI(K+28,L) = FP0(K,IPAIR)
*         ZMPI(K+32,L) = FD0(K,IPAIR)
*  71     CONTINUE
*         ZMPI(37,L) = H(IPAIR)
*         ZMPI(38,L) = HDOT(IPAIR)
*         ZMPI(39,L) = HDOT2(IPAIR)
*         ZMPI(40,L) = HDOT3(IPAIR)
*         ZMPI(41,L) = HDOT4(IPAIR)
*         ZMPI(42,L) = DTAU(IPAIR)
*         ZMPI(43,L) = TDOT2(IPAIR)
*         ZMPI(44,L) = TDOT3(IPAIR)
*         ZMPI(45,L) = R(IPAIR)
*         ZMPI(46,L) = R0(IPAIR)
*         ZMPI(47,L) = GAMMA(IPAIR)
*         ZMPI(48,L) = H0(IPAIR)
*         DO 72 K = 1,7
*         ZMPI(K+48,L) = SF(K,IPAIR)
*72       CONTINUE
*         IMPI(L) = KSLOW(IPAIR)
*
* Do not forget to broadcast HT 
*
*       Set KS indicator on termination, multiple regularization or merger.
              IF (IPHASE.NE.0) THEN
                  IF (IQ.EQ.0.OR.IPHASE.LT.0) THEN
                      IQ = IPHASE
*       Save KS index until exit (collision treated in situ).
                      IF (IQ.GT.0) THEN
                          I10 = I1
                      END IF
                  END IF
*
*       Reset non-zero decision indicator (continue on positive value).
                  IF (IPHASE.GT.0) THEN
                      IPHASE = 0
                  ELSE
*       Enforce new sorted list on change of KS sequence after collision.
                      IPHASE = 0
                      TBLIST = TIME
                      PRINT*,' WARNING COLLISION '
                      STOP
                      GO TO 1
                  END IF
              END IF
 7         CONTINUE
*
*         DO 75 L = 1,KSB
*
*         I1 = BSLIST(L)
*         IPAIR = KVEC(I1)
*
*         DO 76 K = 1,4
*         U(K,IPAIR) = ZMPI(K,L)
*         U0(K,IPAIR) = ZMPI(K+4,L)
*         UDOT(K,IPAIR) = ZMPI(K+8,L)
*         FU(K,IPAIR) = ZMPI(K+12,L)
*         FUDOT(K,IPAIR) = ZMPI(K+16,L)
*         FUDOT2(K,IPAIR) = ZMPI(K+20,L)
*         FUDOT3(K,IPAIR) = ZMPI(K+24,L)
*         FP0(K,IPAIR) = ZMPI(K+28,L)
*         FD0(K,IPAIR) = ZMPI(K+32,L)
*  76     CONTINUE
*         H(IPAIR) = ZMPI(37,L)
*         HDOT(IPAIR) = ZMPI(38,L)
*         HDOT2(IPAIR) = ZMPI(39,L)
*         HDOT3(IPAIR) = ZMPI(40,L)
*         HDOT4(IPAIR) = ZMPI(41,L)
*         DTAU(IPAIR) = ZMPI(42,L)
*         TDOT2(IPAIR) = ZMPI(43,L)
*         TDOT3(IPAIR) = ZMPI(44,L)
*         R(IPAIR) = ZMPI(45,L)
*         R0(IPAIR) = ZMPI(46,L)
*         GAMMA(IPAIR) = ZMPI(47,L)
*         H0(IPAIR) = ZMPI(48,L)
*         DO 77 K = 1,7
*         SF(K,IPAIR) = ZMPI(K+48,L)
*77       CONTINUE
*         KSLOW(IPAIR) = IMPI(L)
*
*75       CONTINUE

*       Prepare pointer and time for INSERT routine
           LIOLD = LI
           LI = LI + KSB
*
        PRINT*,'--------------------'
        PRINT*,' Before Insert LIOLD,KSB,LI,TBLIST ',LIOLD,KSB,LI,TBLIST
        CALL FLUSH(6)
      J1 = KBLIST(LI)
      JNNTB = KBLIST(NNTB)
      PRINT*,' KBL ',J1,JNNTB,T0(J1)+STEP(J1),T0(JNNTB)+STEP(JNNTB)
*       See whether current pair is due before new KBLIST loop.
*       Insert body #I1 in the correct sequential location.
           DO 8 L = 1,KSB
*
           I1 = BSLIST(L)
*
              IF (T0(I1) + STEP(I1).LT.TBLIST) THEN
                  TXLIST(I1) = T0(I1) + STEP(I1)
                  TIME = TXLIST(I1)
                  PRINT*,' Call Insert I1,TXLIST=',I1,TXLIST(I1)
                  CALL INSERT(I1,LI,TXLIST)
              END IF
*
 8         CONTINUE
*
 9         SPREV = SBLOCK
*
 5         CONTINUE
*       Copy original block time at end of KS treatment.
          TIME = TBLOCK
          NBPREV = NPAIRS
      END IF
*
*       Check time for advancing any triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
   30     TSUB = 1.0D+10
          DO 40 L = 1,NSUB
              IF (TS(L).LT.TSUB) THEN
                  ISUB = L
                  TSUB = TS(L)
              END IF
   40     CONTINUE
*
          IF (TSUB.LE.TBLOCK) THEN
              TIME = TSUB
*       Decide between triple, quad or chain.
              IF (ISYS(ISUB).EQ.1) THEN
*       Update unperturbed size of subsystem and copy c.m. step.
                  CALL EXTEND(ISUB)
                  CALL TRIPLE(ISUB)
              ELSE IF (ISYS(ISUB).EQ.2) THEN
                  CALL EXTEND(ISUB)
                  CALL QUAD(ISUB)
              ELSE
                  IF (STEPS(ISUB).LT.0.0D0) THEN
                      STEPS(ISUB) = 1.0D-10
                      GO TO 50
                  END IF
                  CALL CHAIN(ISUB)
                  IF (ISUB.GT.0.AND.STEPS(ISUB).LT.0.0D0) THEN
                      STEPS(ISUB) = 1.0D-10
                      GO TO 50
                  END IF
              END IF
*
*       Check for termination (set TPREV < TIME and set IQ < 0).
              IF (ISUB.LT.0.OR.IPHASE.LT.0) THEN
                  TPREV = TIME - STEP(NTOT)
                  IQ = -1
              END IF
              GO TO 30
          END IF
   50     TIME = TBLOCK
      END IF
*
      RETURN
*
      END

