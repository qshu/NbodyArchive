      SUBROUTINE SIEVE(NCBLK,NCPERT,ICBLK,ICPERT)
*
*
*        Merging of block perturber lists.
*        ---------------------------------
*
      INCLUDE 'common4.h'
      INTEGER  ICBLK(KMAX),ICPERT(KMAX)
*
*
*       Increase counters for single & multiple active c.m..
      NCBLK1 = NCBLK1 + 1
      IF (NCBLK.GT.1) NCBLK2 = NCBLK2 + 1
*
*       Determine index & KS membership of c.m. with maximum perturbers.
      DO 1 L = 1,NCBLK
          JL = 2*(ICBLK(L) - N) - 1
          IF (LIST(1,JL).GT.NSPERT) THEN
              LM = L
              NSPERT = LIST(1,JL)
          END IF
    1 CONTINUE
*
*       Ensure that c.m. body #1 has the largest number of perturbers.
      IF (LM.GT.1) THEN
          K = ICBLK(1)
          ICBLK(1) = ICBLK(LM)
          ICBLK(LM) = K
      END IF
*
*       Define list index of KS with maximum membership.
      J1 = 2*(ICBLK(1) - N) - 1
*
*       Copy all perturbers of ICBLK(1) to joint COMMON array.
      DO 5 L = 1,NSPERT
          JPERT(L) = LIST(L+1,J1)
    5 CONTINUE
*
*       Reduce NSPERT to location of last single perturber.
      L2 = NSPERT
   10 IF (JPERT(NSPERT).GT.N) THEN
          NSPERT = NSPERT - 1
          IF (NSPERT.GT.0) GO TO 10
      END IF
*
*       Increase NSPERT to include any unperturbed pairs (sequential).
      DO 12 L = NSPERT+1,L2
          JL = 2*(JPERT(L) - N) - 1
          IF (LIST(1,JL).EQ.0) THEN
              NSPERT = NSPERT + 1
              JPERT(NSPERT) = JPERT(L)
          END IF
   12 CONTINUE
*
*       Save last location for sequential list comparisons.
      LAST = NSPERT
*
*       Add all active c.m. after last single/inactive body (not sequential).
      DO 15 L = 1,NCPERT
          NSPERT = NSPERT + 1
          JPERT(NSPERT) = ICPERT(L)
   15 CONTINUE
*
*       Check all mutual c.m. distances for overlap (inner loop to LL - 1).
      DO 70 LL = 2,NCBLK
          J = ICBLK(LL)
          JPAIR = J - N
*       Set maximum size using (small) apocentre distance or R0, R & 2*RMIN.
          APO = -BODY(J)/H(JPAIR)
          IF (APO.GT.0.0.AND.APO.LT.0.1*RMIN) THEN
              RJ = APO
          ELSE
              RJ = MAX(R0(JPAIR),R(JPAIR),4.0*RMIN)
          END IF
*
*       Specify KS index and membership of any c.m. (> #1) to be checked.
          J1 = 2*JPAIR - 1
          NP1 = LIST(1,J1) + 1
*
*       Examine previous c.m. positions only (total of NC*(NC - 1)/2 terms).
          DO 50 KK = 1,LL-1
              I = ICBLK(KK)
              IPAIR = I - N
              APO = -BODY(I)/H(IPAIR)
              IF (APO.GT.0.0.AND.APO.LT.0.1*RMIN) THEN
                  RI = APO
              ELSE
                  RI = MAX(R0(IPAIR),R(IPAIR),4.0*RMIN)
              END IF
*
*       Include mass factor (BODY1/BODYM)**(2/3) from perturber selection.
              RS2 = BODY23*CMSEP2*(RI + RJ)**2
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2
     &                                    + (X(3,I) - X(3,J))**2
*       Compare c.m. separation with overlapping distance.
              IF (RIJ2.LT.RS2) THEN
                  K1 = 1
                  DO 40 L = 2,NP1
                      JP = LIST(L,J1)
*       Split search into two parts (first list and any subsequent).
   20                 IF (K1.GT.LAST) GO TO 30
                      K = K1
                      IF (JP.EQ.JPERT(K)) THEN
                          K1 = K + 1
                          GO TO 40
                      ELSE IF (JP.GT.JPERT(K)) THEN
                          K1 = K + 1
                          GO TO 20
                      END IF
*       Skip next part if current perturber is contained in JPERT.
   30                 DO 35 K = LAST+1,NSPERT
                          IF (JP.EQ.JPERT(K)) GO TO 40
   35                 CONTINUE
*       Increase membership and save new perturber.
                      NSPERT = NSPERT + 1
                      JPERT(NSPERT) = JP
   40             CONTINUE
                  GO TO 70
              END IF
   50     CONTINUE
*
*       Reduce list index to last single perturber (L2 = 1 is OK).
          L2 = NP1
   55     IF (LIST(L2,J1).GT.N) THEN
              L2 = L2 - 1
              GO TO 55
          END IF
*
*       Add all single particle perturbers of JPAIR to merged list.
          DO 60 L = 2,L2
              J = LIST(L,J1)
              NSPERT = NSPERT + 1
              JPERT(NSPERT) = J
   60     CONTINUE
*
*       Include any inactive c.m. perturbers.
          DO 65 L = L2+1,NP1
              J = LIST(L,J1)
              JL = 2*(J - N) - 1
              IF (LIST(1,JL).EQ.0) THEN
                  NSPERT = NSPERT + 1
                  JPERT(NSPERT) = J
              END IF
   65     CONTINUE
   70 CONTINUE
*
      RETURN
*
      END
