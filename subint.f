      SUBROUTINE SUBINT(IQ,I10)
*
*
*       Decision-making for subsystems.
*       -------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      REAL*8  TSLIST(10*KMAX)
      SAVE  IRUN,LI
      DATA  IRUN /0/
      COMMON/BARR/ ibarcount

*
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
          LI = LI - 1
          DTB = 0.0
      END IF
*
*       See whether to advance any KS solutions at start of block-step.
    1 IF (NPAIRS.GT.0) THEN
*       Obtain list of all KS pairs due in interval DTB.
          IF (TBLIST.LE.TBLOCK.OR.NPAIRS.NE.NBPREV) THEN
              IF (DTB.EQ.0.0D0.OR.DTB.GT.1.0D+06) THEN
                  DTB = MAX(DTMIN,TBLOCK - TPREV)
              END IF
    2         TBLIST = TPREV + DTB
              TBLIST = MAX(TBLOCK,TBLIST)
              NNTB = 0
              DO 3 JPAIR = 1,NPAIRS
                  J1 = 2*JPAIR - 1
                  IF (T0(J1) + STEP(J1).LE.TBLIST) THEN
                      NNTB = NNTB + 1
                      KBLIST(NNTB) = J1
                      TSLIST(NNTB) = T0(J1) + STEP(J1)
*                      JI = JPAIR + N
*                      CALL JPRED(JPAIR+N,TIME,TIME)
*                      IF (LIST(1,JI).GE.1) THEN
*                         JN = LIST(1,JI) + 1
*                         DO JL = 2, JN
*                            JJ = LIST(JL,JI)
*                            call JPRED(JJ,TPREV,TPREV)
*                         END DO
*                      END IF
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
              LI = 0
          END IF
*     --09/04/13 15:10-lwang-debug--------------------------------------**
***** Note:------------------------------------------------------------**
*          if (time.ge.82.2567) then
*          do iii=1,ifirst
*             if(name(iii).eq.499) then
*          write(102+rank,*)'499x',iii,t0(iii),step(iii),ibarcount,time
*             print*,rank,'NNTB ',nntb,time,ibarcount
*             print*,rank,'NLS',(KBLIST(K),K=1,NNTB)
*             call flush(6)
*          end if
*          end do
*           if(time.ge.6.52160644531250000E-002) then
*           end if
*          if(time.ge.6E-002) print*,rank,time
*          if(time.ge.6.6E-02) stop
*     --09/04/13 15:10-lwang-end----------------------------------------*
*       Select members of KBLIST sequentially and set new time.
    5     LI = LI + 1
          IF (LI.GT.NNTB) THEN
*       Enforce KBLIST at TBLOCK in case TIME > TBLOCK (unpert to unpert).
              TBLIST = TBLOCK
              GO TO 1
          END IF
          I1 = KBLIST(LI)
          TIME = T0(I1) + STEP(I1)
*
*       See whether the smallest KS time falls due before next block-step.
          IF (TIME.LE.TBLOCK) THEN
*     --09/04/13 14:37-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*      call mpi_barrier(MPI_COMM_WORLD,ierr)
*      ibarcount = ibarcount + 1
*        if(time.ge.6E-002.and.name(i1).eq.499) then
*        if(tprev.ge.1.41296386718750000E-002.and.name(i1).eq.5287) then
*      if (time.ge.82.25.and.i1.eq.41) then
*          write(102+rank,*),'5287k',i1,t0(i1),step(i1),
*     &         ibarcount,time
*          call flush(6)
*           print*,rank,'5287k',i1,t0(i1),step(i1),time
*       end if
*     --09/04/13 14:37-lwang-end-debug----------------------------------*
   10         CALL KSINT(I1)
*     --09/25/13 14:31-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*              if(name(i1).eq.5287) then
*          write(102+rank,*),'5287kb',i1,t0(i1),step(i1),
*     &         ibarcount,time
*       end if
c$$$      if(tprev.ge.1.41296386718750000E-002.and.name(i1).eq.5287) then
c$$$         print*,rank,'5287kb',i1,t0(i1),step(i1),time
c$$$         call flush(6)
c$$$      end if
*     --09/25/13 14:31-lwang-end----------------------------------------*

*
*       Check for multiple calls of #I1 (saves using CALL INSERT).
              IF (LI.LT.NNTB.AND.IPHASE.EQ.0) THEN
                  TI = TIME + STEP(I1)
                  JX = KBLIST(LI+1)
                  TX = T0(JX) + STEP(JX)
                  IF (TI.LT.TX.AND.TI.LE.TBLOCK) THEN
                      TIME = TI
                      GO TO 10
                  END IF
              END IF
*
*       See whether current pair is due before new KBLIST loop.
              IF (TIME + STEP(I1).LT.TBLIST) THEN
*
*       Form new list if expanded size is too big.
                  IF (NNTB.GE.10*(KMAX-5)) THEN
                      TBLIST = TIME
*       Include test on non-zero indicator (bad luck case).
                      IF (IPHASE.NE.0.AND.IQ.EQ.0) THEN
                          IQ = IPHASE
                          I10 = I1
                          IPHASE = 0
                      END IF
                      GO TO 1
                  END IF
*
*       Insert body #I1 in the correct sequential location (skip T0 < TIME).
                  IF (T0(I1).GE.TIME) THEN
                      CALL INSERT(I1,LI)
                  END IF
              END IF
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
                      GO TO 1
                  END IF
              END IF
*
*       Continue cycle until end of block-step.
              GO TO 5
          END IF
*
*       Copy original block time at end of KS treatment.
          TIME = TBLOCK
          NBPREV = NPAIRS
*       Reduce pointer by 1 for next block-step (otherwise not done).
          LI = LI - 1
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
                  IF (ISUB.GT.0) THEN
                      IF (STEPS(ISUB).LT.0.0D0) THEN
                          STEPS(ISUB) = 1.0D-10
                          GO TO 50
                      END IF
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
