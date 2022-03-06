      SUBROUTINE IMF2(BODY10,BODYN)
*
*
*       Mass function for binaries & single stars.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
      REAL  RHO,RAN2
      COMMON/WORK1/  RHO(NMAX)
*
*=========================
* KTG3 MF with alpha1=1.3 :
      DATA  G1,G2,G3,G4  /0.19,1.55,0.050,0.6/
* KTG3 MF with alpha1=1.1 :
c      DATA  G1,G2,G3,G4  /0.28,1.14,0.010,0.1/
*=========================
*
*
*       Generate initial mass function (KZ(20) = 2 or > 2).
      KDUM = IDUM1
      ZMASS = 0.0D0
      DO 10 I = 1,N+NBIN0
    5     XX = RAN2(KDUM)
*
*       Choose between Kroupa et al (M.N. 262, 545) & Eggleton (book).
          IF (KZ(20).EQ.2.OR.KZ(20).EQ.4) THEN
              ZM = 0.08 + (G1*XX**G2 + G3*XX**G4)/(1.0 - XX)**0.58
          ELSE
              ZM = 0.3*XX/(1.0 - XX)**0.55
          END IF
*
*       See whether the mass falls within the specified range.
          IF (ZM.GT.BODYN.AND.ZM.LT.BODY10) THEN
              BODY(I) = ZM
              ZMASS = ZMASS + BODY(I)
          ELSE
              GO TO 5
          END IF
   10 CONTINUE
*
*       See whether to skip mass function for binaries.
      IF (NBIN0.EQ.0) GO TO 50
*
*       Merge binary components in temporary variable for sorting.
      DO 20 I = 1,NBIN0
          RHO(I) = BODY(2*I-1) + BODY(2*I)
          JLIST(I) = I
   20 CONTINUE
*
*       Sort total binary masses in increasing order.
      CALL SORT1(NBIN0,RHO,JLIST)
*
*       Save scaled binary masses in decreasing order for routine BINPOP.
      DO 30 I = 1,NBIN0
          JB = JLIST(NBIN0-I+1)
          BODY0(2*I-1) = MAX(BODY(2*JB-1),BODY(2*JB))/ZMASS
          BODY0(2*I) = MIN(BODY(2*JB-1),BODY(2*JB))/ZMASS
          IF (KZ(20).GT.3) THEN
*       Adopt correlation (m1/m2)' = (m1/m2)**0.4 & constant sum (Eggleton).
              ZMB = BODY0(2*I-1) + BODY0(2*I)
              RATIO = BODY0(2*I-1)/BODY0(2*I)
              BODY0(2*I) = ZMB/(1.0 + RATIO**0.4)
              BODY0(2*I-1) = ZMB - BODY0(2*I)
          END IF
   30 CONTINUE
*
*       Merge binary components into single stars for scaling purposes.
      ZMB = 0.0
      DO 40 I = 1,NBIN0
          BODY(NBIN0-I+1) = RHO(I)
          ZMB = ZMB + RHO(I)
   40 CONTINUE
*
      WRITE (6,45)  NBIN0, BODY(1), BODY(NBIN0), ZMB/FLOAT(NBIN0)
   45 FORMAT (/,12X,'BINARY STAR IMF:   NB =',I5,'  RANGE =',1P,2E10.2,
     &              '  <MB> =',E9.2)
*
*       Move the single stars up to form compact array of N members.
C Addition following nb5 (Aug.1998, P.Kroupa)
   50 IF (N.LE.NBIN0) GO TO 90
      NS = 0
      DO 60 I = NBIN0+1,N
          BODY(I) = BODY(I+NBIN0)
          NS = NS + 1
          RHO(NS) = BODY(I)
          JLIST(NS) = I
   60 CONTINUE
*
*       Sort masses of single stars in increasing order.
      CALL SORT1(NS,RHO,JLIST)
*
*       Copy the masses of single stars to COMMON in decreasing order.
      ZMS = 0.0
      DO 70 I = 1,NS
          BODY(N-I+1) = RHO(I)
          ZMS = ZMS + RHO(I)
   70 CONTINUE
*
      WRITE (6,80)  N-NBIN0, BODY(NBIN0+1), BODY(N), ZMS/FLOAT(N-NBIN0)
   80 FORMAT (/,12X,'SINGLE STAR IMF:   NS =',I5,'  RANGE =',1P,2E10.2,
     &            '  <MS> =',E9.2)
*
*       Replace input value by actual mean mass in solar units.
   90 ZMBAR = ZMASS/FLOAT(N)
*
      IDUM1 = KDUM
*
      RETURN
*
      END
