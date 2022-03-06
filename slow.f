      SUBROUTINE SLOW(Y)
*
*
*       Slow-down treatment of chain binary.
*       ------------------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      LOGICAL  KSLOW,KCOLL
      REAL*8  Y(NMX8),KSCH,KSNEW
      COMMON/SLOW1/   TK2(0:NMX),EJUMP,KSCH(NMX),KSLOW,KCOLL
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/KSAVE/   K1,K2
      COMMON/SLOW3/  GCRIT,KZ26
      SAVE
*
*
*       Perform fast perturbation check if slow-down not active.
      IF (.NOT.KSLOW) THEN
*       Determine chain index and largest inverse distance.
          RM = 0.0
          DO I = 1,N-1
              IF (RINV(I).GT.RM) THEN
                  RM = RINV(I)
                  I1 = I
              END IF
          END DO
*
*       Sum the perturbating forces m/r^3 next to #i1.
          sum = 0.0
          do i = 1,n-1
              if (iabs(i-i1).eq.1) then
                  j = i
                  if (i.gt.i1) j = i + 1
                  sum = sum + mc(j)*rinv(i)**3
              end if
          end do
*
*       Skip further search if relative perturbation exceeds limit.
          mb = mc(i1) + mc(i1+1)
          GAMMA = 2.0*sum/(mb*RM**3)
          IF (GAMMA.GT.GCRIT) THEN
              GO TO 50
          END IF
      END IF
*
*       Save the current variables in common (RINV is OK after switching).
      CALL YSAVE(Y)
*
*       Determine chain index for slow-down binary (otherwise from above).
      IF (KSLOW) THEN
          DO I = 1,N-1
              IF (KSCH(I).GT.1.0) i1 = i
          END DO
      END IF
*
*       Check for switching off slow-down at start of iteration.
      if (GCRIT.eq.0.0d0) then
          ksnew = 1.0
          go to 30
      end if
*
*       Save QK & PK and copy current configuration for EREL & TRANSK.
      DO 10 I = 1,N-1
          KS = 4*(I - 1)
          DO 5 J = 1,4
              QK(KS+J) = Q(KS+J)
              PK(KS+J) = P(KS+J)
    5     CONTINUE
   10 CONTINUE
*
*       Evaluate semi-major axis from non-singular variables.
      K1 = INAME(i1)
      K2 = INAME(i1+1)
      CALL EREL(i1,EB,SEMI)
*
*       Exit if no current binary (set KSLOW = .false. just in case).
      IF (SEMI.LE.0.0d0) THEN
          KSLOW = .false.
          TK2(0) = 0.0
          GO TO 50
      END IF
*
*       Sum the perturbations next to #i1 (already known if KSLOW active).
      IF (KSLOW) THEN
          sum = 0.0
          do i = 1,n-1
              if (iabs(i-i1).eq.1) then
                  j = i
                  if (i.gt.i1) j = i + 1
                  sum = sum + mc(j)*rinv(i)**3
              end if
          end do
      END IF
*
*       Form relative perturbation at maximum apocentre.
      rap = 2.0*SEMI
      pert = 2.0*sum*rap**3/(mc(i1) + mc(i1+1))
*
*       Specify the slow-down factor.
      IF (pert.LT.GCRIT) THEN
          ksnew = SQRT(GCRIT/pert)
      ELSE
          ksnew = 1.0
      END IF
*
*     --------------------------------------------------------------
*       Implement discrete scheme (suppressed).
*     i = i1
*     rat = ksch(i)
*     rat = ksnew/rat
*       Check for significant changes (by 2) in slow-down factor.
*     if (rat.gt.0.5.and.rat.le.2.0) then
*         ksnew = Ksch(i)
*     else if (rat.gt.2.0) then
*         ksnew = 2*Ksch(i)
*       Allow an extra factor of 32 at initialization.
*         if (.not.KSLOW) then
*             if (rat.ge.64.0) ksnew = 32*ksnew
*         end if
*     else if (rat.le.0.5) then
*         ksnew = Ksch(i)/2
*         if (rat.le.0.125) ksnew = Ksch(i)/4
*         ksnew = max(ksnew,1.0D0)
*     end if
*     --------------------------------------------------------------
*
*       Check slow-down switch and include any change in binding energy.
   30 if (ksnew.ne.ksch(i1)) then
          if (ksnew.eq.1.0d0) then
              KSLOW = .false.
          else
              KSLOW = .true.
          end if
*       Add change in binary energy and save new slow-down index.
          eb = -0.5d0*MC(i1)*MC(i1+1)/SEMI
          DEB = eb*(1.0/ksnew - 1.0/Ksch(i1))
          Ksch(i1) = ksnew
      else
          DEB = 0.0
      end if
*
*       Accumulate total energy difference due to slow-down.
      EJUMP = EJUMP + DEB
*       Ensure zero EJUMP without slow-down to avoid small residual.
      IF (.NOT.KSLOW) EJUMP = 0.0
*
*       Form modified mass factors.
      do i = 1,n
          tk1(i) = -1.0/MC(I)
      end do
      do i = 1,n-1
          if (Ksch(i).ne.1.0d0) then
              tk1(i) = tk1(i)/Ksch(i)
              tk1(i+1) = tk1(i+1)/Ksch(i)
          end if
      end do
      DO I = 1,N-1
          TKK(I) = 0.5D0*(-tk1(i) - tk1(i+1))
          MKK(I) = MC(I)*MC(I+1)/Ksch(i)
      END DO
      do i = 1,n-1
          m12 = mc(i) + mc(i+1)
          dt12 = 0.5d0*(1.0d0 - 1.0d0/Ksch(i))/m12
          if (i.gt.1) TKK(i-1) = tkk(i-1) + dt12
          if (i.lt.n-1) TKK(i+1) = tkk(i+1) + dt12
          if (i.gt.1.and.i.lt.n-1) TK2(i) = -2.0d0*dt12
      end do
      TK2(0) = 0.0
      TK2(N) = 0.0
*
   50 RETURN
*
      END
