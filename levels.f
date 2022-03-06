      SUBROUTINE LEVELS
*
*
*       Diagnostics for block time-step levels.
*       ---------------------------------------
*       plus estimation of parallel speed-up (R.Sp.)
*
      INCLUDE 'common6.h'
      INTEGER JHIST,JHISTR
      COMMON/BLKLVL/  JHIST(NMAX),JHISTR(NMAX)
      INTEGER  IHIST(NMAX),IHISTR(NMAX)
      INTEGER IPES,IPROC(9),IY(1024),IYR(1024)
      REAL*8 XSPEED(9),XSPEDR(9)
      DATA IPROC/4,8,16,32,64,128,256,512,1024/
*
*
*       Initialize histogram counters.
      DO 10 J = 1,NMAX
          IHIST(J) = 0
          IHISTR(J) = 0
   10 CONTINUE
*
*       Loop over all single particles & c.m.
      JMAX = 0
      JMAXR = 0
      FAC = 1.0/LOG(1.9999999)
      DO 20 I = IFIRST,NTOT
          J = 1 - LOG(STEP(I))*FAC
          IF(J.GT.NMAX)J = NMAX
          IHIST(J) = IHIST(J) + 1
          JMAX = MAX(J,JMAX)
          J = 1 - LOG(STEPR(I))*FAC
          IF(J.GT.NMAX)J = NMAX
          IHISTR(J) = IHISTR(J) + 1
          JMAXR = MAX(J,JMAXR)
   20 CONTINUE
*
*       Print histograms of block-steps (STEPR with KZ(34) > 1).
      JMAX=MIN(JMAX,NMAX)
      JMAXR=MIN(JMAXR,NMAX)
      if(rank.eq.0)then
	  WRITE (6,30)  (IHIST(J),J=1,JMAX)
   30 FORMAT (' STEP I ',22I5,(/,24I5))
      IF (KZ(34).GT.1)WRITE (6,301)  (IHISTR(J),J=1,JMAXR)
  301 FORMAT (' STEP R ',22I5,(/,24I5))
      end if
*
*       Estimate Parallelism for different Processor Numbers
*
*       IPROC contains a list of possible processor numbers,
*       presently   4 to 1024 in powers of 2.
*
      DO 35 NN = 1,9
      IPES = IPROC(NN)
      DO 36 I = 1,IPES
      IY(I) = 0
 36   IYR(I) = 0
*
*       JHIST contains accumulated frequencies from INTGRT for irr. steps. 
*       JHIST(80) = 5 e.g. means that it occurred 5 times since the last
*                   call of LEVELS that 80 force calculations had
*                   to be done in parallel
*       IY maps this information to a certain number of processors,
*       e.g. in the above case for 64 processors it would mean that
*       five times 64 processors could be fully used, and five times
*       only 16 processors (to get the remaining forces).
*       JHIST(80) = 5 results in incrementing IY(64) and IY(16) by 5.
*
*       (IYE, JHISTR do the same for the regular steps)
*
      DO 37 I = 1,NMAX
      IF(JHIST(I).EQ.0)GOTO 37
      J = MOD(I,IPES)
      K = I/IPES
      IF(J.GT.0)IY(J)=IY(J)+JHIST(I)
      IF(K.GT.0)IY(IPES)=IY(IPES)+K*JHIST(I)
 37   CONTINUE
*
      DO 371 I = 1,NMAX
      IF(JHISTR(I).EQ.0)GOTO 371
      J = MOD(I,IPES)
      K = I/IPES
      IF(J.GT.0)IYR(J)=IYR(J)+JHISTR(I)
      IF(K.GT.0)IYR(IPES)=IYR(IPES)+K*JHISTR(I)
 371  CONTINUE
*
      ISPEED = 0
      ISPEDR = 0
      ITOT = 0
      ITOTR = 0
      DO 38 I = 1,IPES
      ISPEED = ISPEED + I*IY(I)
      ISPEDR = ISPEDR + I*IYR(I)
      ITOT = ITOT + IY(I)
      ITOTR = ITOTR + IYR(I)
 38   CONTINUE
*
*       Estimate theoretical speedup by (sum I*IY) / (sum IY);
*       enumerator is the total number of force calculations necessary;
*       denominator is the number of force calculation `cycles' necessary
*       on a parallel machine. Separately done for regular/irregular
*       steps. Communication not yet included in any way.
*
      IF(ITOT.NE.0)XSPEED(NN) = REAL(ISPEED)/REAL(ITOT)
      IF(ITOTR.NE.0)XSPEDR(NN) = REAL(ISPEDR)/REAL(ITOTR)
 35   CONTINUE
*
	  if(rank.eq.0)then
      WRITE (6,40)  (IPROC(J),XSPEED(J),J=1,9)
 40   FORMAT (' Max Speedup Irr: ',1P,9(I5,D9.2))
      WRITE (6,401)  (IPROC(J),XSPEDR(J),J=1,9)
 401  FORMAT (' Max Speedup Reg: ',1P,9(I5,D9.2))
      end if
*
      DO 50 J = 1,NMAX
          JHIST(J) = 0
          JHISTR(J) = 0
   50 CONTINUE
*
      RETURN
*
      END
