      SUBROUTINE DELAY(IDUMMY,IPAIR)
*
*
*       Delay of multiple regularization.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'kscmn6.h'
      SAVE  IPH0,KS0,KS20,JC0,JCL0,JCMAX0,EB0
*
*
*       Check for saving KS pointers from SUBINT or termination at TBLOCK.
      IF (IPAIR.EQ.0) THEN
*       Reset state.
          IPH0 = 0
      ELSE IF (IPAIR.GT.0) THEN
*       Save new state indices.
*
*       Avoid merger (IPHASE = 6) if chain indicator set by KSINT.
          IF (IPH0.EQ.8.AND.IPHASEX(IPAIR).EQ.6) GO TO 10
*
          IPH0 = IPHASEX(IPAIR)
          KS0 = KSPAIRX(IPAIR)
          KS20 = KS2X(IPAIR)
          JC0 = JCOMPX(IPAIR)
          JCL0 = JCLOSEX(IPAIR)
          JCMAX0 = JCMAXX(IPAIR)
          EB0 = EBCH0X(IPAIR)
      ELSE
*       Copy relevant indices for delayed KS termination.
          IPHASE = IPH0
          KSPAIR = KS0
          KS2 = KS20
          JCOMP = JC0      ! JCOMP copied back from KSINT (still = 0).
          JCLOSE = JCL0    ! JCLOSE copied back from labelled common.
          JCMAX = JCMAX0
          EBCH0 = EB0
*      Preserve contents of KSAVE during chain regularization.
          IF (NCH.EQ.0) THEN
              KSAVE(1) = 0
              KSAVE(2) = 0
          END IF
          KSKIP = 0
*
*      Exit in case of new merger or merger termination.
          IF (IPHASE.EQ.6.OR.IPHASE.EQ.7) THEN
              GO TO 10
          END IF
*
*       *** INPUT ***
*       Binary                    KSPAIR>0  KS2=0  JCLOSE=0     JCMAX=0
*       Binary + single           KSPAIR>0  KS2=0  0<JCLOSE<=N  JCMAX=0
*       Binary + binary           KSPAIR>0  KS2>0  JCLOSE>N     JCMAX=0
*       Binary + single + single  KSPAIR>0  KS2=0  0<JCLOSE<=N  0<JCMAX<=N
*       Binary + binary + single  KSPAIR>0  KS2>0  JCLOSE>N     0<JCMAX<=N
*       Binary + binary + binary  KSPAIR>0  KS2>0  JCLOSE>N     JCMAX>N
*       Binary + inert binary     KSPAIR>0  KS2=0  JCLOSE>N     JCMAX=0
*
*       *** OUTPUT ***
*       Binary                    KSPAIR>0  JCLOSE=0     JCMAX=0
*       Binary + single           KSPAIR>0  0<JCLOSE<=N  JCMAX=0
*       Binary + binary           KSPAIR>0  JCLOSE>N     JCMAX=0
*       Binary + single + single  KSPAIR>0  0<JCLOSE<=N  0<JCMAX<=N
*       Binary + binary + single  KSPAIR>0  JCLOSE>N     0<JCMAX<=N
*       Binary + binary + binary  KSPAIR>0  JCLOSE>N     JCMAX>N (unterminated)
*       Binary + inert binary     KSPAIR=0  JCLOSE>N (unterm.)  JCMAX=0
*
*       Include the case of two interacting KS pairs.
          IF (JCLOSE.GT.N) THEN
              IF (IPHASE.EQ.8.AND.KSTAR(N+KSPAIR).NE.0) THEN
                  KSAVE(1) = KSTAR(N+KSPAIR)
                  KSAVE(2) = NAME(2*KSPAIR-1) + NAME(2*KSPAIR)
                  KSKIP = 1
              END IF
*       Terminate smallest pair first and copy second pair index.
              CALL KSTERM
              IF (JCLOSE.GT.N+KSPAIR) JCLOSE = JCLOSE - 1
              IF (JCMAX.GT.N+KSPAIR) JCMAX = JCMAX - 1
              IF (KS2.GT.KSPAIR) KS2 = KS2 - 1
              KSPAIR = KS2
          END IF
*
*       Save KSTAR (> 0) and sum of component names (for chain termination).
          IF (IPHASE.EQ.8.AND.KSTAR(N+KSPAIR).NE.0.AND.
     &        KSKIP.EQ.0) THEN
              KSAVE(1) = KSTAR(N+KSPAIR)
              KSAVE(2) = NAME(2*KSPAIR-1) + NAME(2*KSPAIR)
          END IF
*
*       Terminate binary in triple or widest binary-binary collision pair.
          IF (KSPAIR.GT.0) THEN  ! Include case of inert binary from IMPACT.
              CALL KSTERM
      IF (IPHASE.NE.IPH0) WRITE (6,99)  IPH0, IPHASE
   99 FORMAT (' DELAY!!!    IPH0 IPH ',2I4)
              IPHASE = IPH0
              IF (JCMAX.GT.N+KSPAIR) JCMAX = JCMAX - 1
          END IF
      END IF
*
   10 RETURN
*
      END
