      SUBROUTINE SETSYS
*
*
*       Selection of chain system.
*       --------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX,1)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/INCOND/  X4(3,4),XDOT4(3,4),ANG(3)
*
*
*       Check whether new chain or addition of member(s).
      IF (NCH.GT.0) GO TO 10
*
*       Initialize chain indices, names & masses for largest interaction.
      DO 1 L = 1,4
          JLIST(L) = 2*NPAIRS + L
          NAMEC(L) = NAME(2*NPAIRS+L)
          BODYC(L) = BODY(2*NPAIRS+L)
          M(L) = BODY(2*NPAIRS+L)
    1 CONTINUE

*       Define chain membership for three-body or four-body case.
      IF (JCLOSE.LE.N) THEN
          NCH = 3
          JLIST(3) = JCLOSE
          NAMEC(3) = NAME(JCLOSE)
          BODYC(3) = BODY(JCLOSE)
          M(3) = BODY(JCLOSE)
      ELSE
          NCH = 4
      END IF
      GO TO 50
*
*       Expand membership and save chain variables (single body or KS pair).
   10 IF (JCLOSE.LE.N) THEN
          NCH = NCH + 1
          JLIST(NCH) = JCLOSE
          NAMEC(NCH) = NAME(JCLOSE)
          BODYC(NCH) = BODY(JCLOSE)
          M(NCH) = BODY(JCLOSE)
*       Improve coordinates & velocities of single perturber (c.m. body OK).
          CALL XVPRED(JCLOSE,-1)
      ELSE
*       Set large indicator to skip termination at block time during chain.
          IPHASE = 10
*       Terminate KS pair and copy components.
          KSPAIR = JCLOSE - N
          CALL KSTERM
          DO 20 L = 1,2
              NCH = NCH + 1
              JLIST(NCH) = 2*NPAIRS + L
              NAMEC(NCH) = NAME(2*NPAIRS+L)
              BODYC(NCH) = BODY(2*NPAIRS+L)
              M(NCH) = BODY(2*NPAIRS+L)
   20     CONTINUE
      END IF
*
*       Specify membership for chain COMMON.
   50 NN = NCH
*
      RETURN
*
      END
