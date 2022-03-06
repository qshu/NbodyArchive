      SUBROUTINE UPDATE(IPAIR)
*
*
*       List modifications after KS termination.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Adjust neighbour lists to new sequence (skip last or only pair).
      ICM = N + IPAIR
      IF (IPAIR.GT.NPAIRS) GO TO 60
*
*       Set largest index of any regularized member.
      ILAST = 2*NPAIRS + 2
*
*       Rename lists containing single regularized components.
      DO 50 J = 1,NTOT
*       Skip renaming if first neighbour exceeds regularized components.
          IF (LIST(2,J).GT.ILAST) GO TO 50
          NNB = LIST(1,J)
          IF (NNB.EQ.0) GO TO 50
          KCOMP = 2
*       Consider the second component before the first.
   10     I = 2*IPAIR + KCOMP - 2
          L = 2
*       Search list for regularized component and rename if > I.
   15     IF (LIST(L,J).NE.I) THEN
              IF (LIST(L,J).GT.I.AND.LIST(L,J).LE.ILAST) THEN
                  LIST(L,J) = LIST(L,J) - 1
              END IF
              L = L + 1
              IF (L.LE.NNB + 1.AND.LIST(L,J).LE.ILAST) GO TO 15
              KCOMP = KCOMP - 1
*       Use same procedure for first component.
              IF (KCOMP.EQ.1) GO TO 10
              GO TO 50
          END IF
*
*       Set new index of old regularized component.
          INEW = 2*NPAIRS + KCOMP
*       Just rename component if last member of the list.
          DO 20 K = L,NNB
              IF (LIST(L+1,J).LE.INEW) THEN
*       Also update the list if INEW happens to be equal to next member.
                  LIST(L,J) = LIST(L+1,J) - KCOMP
*       List member has already been reduced by one if KCOMP = 1.
              END IF
   20     CONTINUE
*
*       Rename old regularized component.
          LIST(L,J) = INEW
   50 CONTINUE
*
*       Replace c.m. by components and reduce subsequent members by one.
   60 DO 80 J = 1,NTOT
          NNB = LIST(1,J)
          L = NNB + 1
*       Check for removal of current c.m. and reduce by one if > N + IPAIR.
   70     IF (L.EQ.1.OR.LIST(L,J).LT.ICM) GO TO 80
          IF (LIST(L,J).NE.ICM) THEN
              LIST(L,J) = LIST(L,J) - 1
              L = L - 1
              GO TO 70
          END IF
*
*       Remove c.m. name by moving up subsequent members.
   72     IF (L.LE.NNB) THEN
              LIST(L,J) = LIST(L+1,J) 
              L = L + 1
              GO TO 72
          END IF
*
          NNB = NNB - 1
*       Expand the list to include both components since c.m. was deleted.
          KCASE = 2
*       Only move neighbours down by one if the list has too many members.
          IF (NNB.GT.LMAX-3) KCASE = 1
          IF (NNB.EQ.0) GO TO 76
*       In this special case L = 2 already.
          L = NNB + 1
*       Take special precaution if last neighbour is a regularized body.
          IF (LIST(L,J).LE.JCOMP) THEN
              L = L + 1
              GO TO 76
          END IF
*
*       Move members down by two (or one) to make room for components.
   74     LIST(L+KCASE,J) = LIST(L,J)
          IF (L.GT.2.AND.LIST(L-1,J).GE.ICOMP) THEN
              L = L - 1
              GO TO 74
          END IF
*
*       Rename deleted c.m. appropriately and increase membership by 2 or 1.
   76     LIST(L,J) = ICOMP
*       Do not over-write the list if NNB > NNBMAX after removal of c.m.
          IF (KCASE.EQ.2) LIST(L+1,J) = JCOMP
          LIST(1,J) = NNB + KCASE
          IF (KCASE.EQ.1) WRITE (6,78)  NNB, J, JCOMP
   78     FORMAT (5X,'WARNING!  TOO LARGE NNB IN UPDATE    NNB J JCOMP',
     &                                                              3I6)
   80 CONTINUE
*
*       Modify the list of previously regularized binaries.
      NNB = LISTR(1) - 1
      L = 0
   91 L = L + 2
   92 IF (L.GT.NNB + 1) GO TO 96
      J = LISTR(L)
      K = LISTR(L+1)
*       First check the current two-body separation of any old pairs.
      RJK2 = (X(1,J) - X(1,K))**2 + (X(2,J) - X(2,K))**2 +
     &                              (X(3,J) - X(3,K))**2
*       Remove pair if RJK > 4*RMIN when special procedure not needed.
      IF (RJK2.LT.16.0*RMIN**2) GO TO 91
      DO 94 K = L,NNB
          LISTR(K) = LISTR(K+2)
   94 CONTINUE
      NNB = NNB - 2
      GO TO 92
*
*       Add ICOMP & JCOMP to LISTR (maximum of MLR/2 - 1 pairs).
   96 IF (NNB.GT.MLR - 4) THEN
*       Note that NNB is one less than the actual membership.
          DO 98 K = 2,NNB
              LISTR(K) = LISTR(K+2)
   98     CONTINUE
          NNB = NNB - 2
*       Removal of the oldest KS pair.
      END IF
      LISTR(NNB+3) = ICOMP
      LISTR(NNB+4) = JCOMP
      LISTR(1) = NNB + 3
*
*       Copy flag index of disrupted pair (set in KSTERM).
      IFLAG = JLIST(3)
*       Add primordial pairs to LISTD (skip new KS pairs).
      IF (IFLAG.EQ.0) GO TO 110
*
*       Check list of disrupted component names.
      NNB = LISTD(1) - 1
      KCOMP = 0
      DO 100 K = 2,NNB+1,2
          IF (LISTD(K).EQ.JLIST(1).AND.LISTD(K+1).EQ.JLIST(2)) KCOMP = 1
  100 CONTINUE
*
*       Include both components unless already members.
      IF (KCOMP.EQ.0) THEN
          IF (NNB.GT.MLD - 4) THEN
              DO 102 K = 2,NNB
                 LISTD(K) = LISTD(K+2)
  102         CONTINUE
              NNB = NNB - 2
          END IF
*       Add most recent names at the end (limit is MLD/2 - 1 pairs).
          LISTD(NNB+3) = JLIST(1)
          LISTD(NNB+4) = JLIST(2)
          LISTD(1) = NNB + 3
      END IF
      IF (IFLAG.NE.-1) WRITE (8,104)  IPAIR, IFLAG, JLIST(1), JLIST(2)
  104 FORMAT (' LISTD INCONSISTENCY!!  IPAIR IFLAG NAMES ',2I5,2I8)
*
*       Remove first component of old KS pair and c.m. body from NLIST.
  110 I = 2*IPAIR - 1
      CALL NLMOD(I,-1)
      CALL NLMOD(ICM,-1)
*
*       Rename any older single KS components and more recent c.m. bodies.
      IF (IPAIR.LE.NPAIRS) THEN
          NNB1 = NLIST(1) + 1
          DO 120 L = 2,NNB1
*       Reduce index of any subsequent c.m. and first components by 1 & 2.
              IF (NLIST(L).GT.ICM) NLIST(L) = NLIST(L) - 1
              IF (NLIST(L).LE.JCOMP.AND.NLIST(L).GT.2*IPAIR - 1)
     &                                           NLIST(L) = NLIST(L) - 2
  120     CONTINUE
      END IF
*
*       Update list of high velocity particles containing c.m. members.
      NNB = LISTV(1)
      DO 130 L = 2,NNB+1
          IF (LISTV(L).EQ.ICM) THEN
*       Remove old c.m. and reduce the membership.
              DO 125 K = L,NNB
                  LISTV(K) = LISTV(K+1)
  125         CONTINUE
              LISTV(1) = LISTV(1) - 1
          END IF
*       Reduce higher particle locations by one.
          IF (LISTV(L).GT.ICM) THEN
              LISTV(L) = LISTV(L) - 1
          END IF
  130 CONTINUE
*
      RETURN
*
      END
