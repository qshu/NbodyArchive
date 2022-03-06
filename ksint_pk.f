      SUBROUTINE KSINT(I1)
*
*
*       Regularized integration.
*       ------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW/  RANGE,ISLOW(10)
      COMMON/GAMDOT/  DGAM
      REAL  W1,W2,W3
      REAL*8  UI(4),UIDOT(4),XI(6),FP(6)
*
*
*       Set second component, pair index & c.m. index.
      I2 = I1 + 1
      IPAIR = KVEC(I2)
      I = N + IPAIR
*
*       Define perturber membership & inverse c.m. mass.
      NNB0 = LIST(1,I1)
      BODYIN = 1.0/BODY(I)
*
*       Update current regularized time (only used as a zero point).
      TAU = T0U(IPAIR) + DTAU(IPAIR)
*
*       Introduce KS slow-down index & modification factor.
      IMOD = KSLOW(IPAIR)
      ZMOD = FLOAT(ISLOW(IMOD))
*
*       Check for further unperturbed motion or reflected restart.
      IF (NNB0.EQ.0) THEN
          KPERT = 0
          IF (TDOT2(IPAIR).LT.0.0D0) THEN
              CALL UNPERT(IPAIR,KPERT)
*       Note KPERT = -1 or +1 for unperturbed (& merger) or perturbed motion.
          END IF
*
*       Form perturber list and restart KS motion (including FREEZE).
          IF (KPERT.GE.0) THEN
              CALL KSLIST(IPAIR)
              IF (LIST(1,I1).GT.0) THEN
                  CALL KSPOLY(IPAIR,1)
              END IF
          END IF
          GO TO 100
      END IF
*
*       Perform KS prediction of U, UDOT & H.
      CALL KSPRED(IPAIR,I1,I,BODYIN,UI,UIDOT,XI)
*
*       Obtain the perturbing force.
      CALL KSPERT(I1,NNB0,XI,FP)
*
*       See if an external force should be added to the perturbation.
      IF (KZ(14).NE.0) THEN
          Q1 = XI(1) - XI(4)
          Q3 = XI(3) - XI(6)
          CALL XTRNLP(Q1,Q3,FP)
      END IF
*
*       Save old radial velocity & relative perturbation and set new GAMMA.
      RDOT = TDOT2(IPAIR)
      G0 = GAMMA(IPAIR)
      GI = SQRT(FP(1)**2 + FP(2)**2 + FP(3)**2)*R(IPAIR)**2*BODYIN
      GAMMA(IPAIR) = GI
*
*       Apply the fourth-order corrector.
      CALL KSCORR(IPAIR,UI,UIDOT,FP,TD2,TDOT4,TDOT5)
*
*       Increase regularization time-step counter and update the time.
      NCOUNT(3) = NCOUNT(3) + 1
      NSTEPU = NSTEPU + 1
      T0(I1) = TIME
*       Check for early return during termination (called from KSTERM).
      IF (IPHASE.NE.0) GO TO 100
*
*       Define useful scalars.
      RI = R(IPAIR)
      HI = H(IPAIR)
*
*       Find the dominant body for large perturbations.
      IF (GI.LT.0.5) GO TO 45
      S = 2.0*STEP(I)
      FMAX = BODY(I)/RI**2
      JCOMP = 0
*       Value of JCOMP is used for optional diagnostics in KSTERM.
      DO 40 L = 2,NNB0+1
          J = LIST(L,I1)
*       Only search bodies within twice the c.m. time-step.
          IF (STEP(J).GT.S) GO TO 40
*       Compare strong perturber and either component with current pair.
          DO 35 K = I1,I2
              RIJ2 = (X(1,J) - X(1,K))**2 + (X(2,J) - X(2,K))**2 +
     &                                      (X(3,J) - X(3,K))**2
              IF (BODY(J) + BODY(K).GT.RIJ2*FMAX) JCOMP = J
   35     CONTINUE
   40 CONTINUE
*
*       Terminate regularization if strong perturber forms dominant pair.
      IF (JCOMP.GT.0.OR.GI.GT.1.0) THEN
*       Check optional binary search.
          IF (KZ(4).GT.0) THEN
              DGAM = GI - G0
              K = KZ(4)
              CALL EVOLVE(IPAIR,K)
          END IF
          GO TO 90
      END IF
*
*       Check termination of hyperbolic encounter (R > R0 or R > RMIN).
   45 IF (HI.GT.0.0D0) THEN
          IF (RI.GT.R0(IPAIR).AND.GI.GT.GMAX) GO TO 90
          IF (RI.GT.RMIN) GO TO 90
      END IF
*
*       Choose basic regularized step using binding energy or lower limit.
      IF (ABS(HI).GT.ECLOSE) THEN
          W1 = 0.5/ABS(HI)
      ELSE
          W1 = R0(IPAIR)*BODYIN
          W2 = 0.5/ABS(HI)
          W1 = MIN(W1,W2)
          IF (RI.GT.R0(IPAIR)) W1 = W1*R0(IPAIR)/RI
*       Maximum square step for soft binaries & weak hyperbolic pairs.
          IF (HI.LT.0.0D0) THEN
*       Include case of hard binary with massive components or merger.
              W2 = -0.5/HI
              W1 = MIN(W2,W1)
              IF (NAME(I).LT.0) W1 = W2
          END IF
      END IF
*
*       Include perturbation factor in predicted step.
      IF (GI.LT.0.0005) THEN
*       Use second-order expansion of cube root for small perturbations.
          W3 = 333.3*GI
          W2 = SQRT(W1)/(1.0 + W3*(1.0 - W3))
      ELSE
          W3 = 1.0 + 1000.0*GI
          W2 = SQRT(W1)/W3**0.3333
      END IF
*
*       Set new regularized step and convert to physical time units.
      DTU = MIN(ETAU*W2,1.2D0*DTAU(IPAIR))
*       Inertial factor of 1.2 used because of smaller initial step.
      STEP(I1) = ((((0.2D0*TDOT5*DTU + 0.5D0*TDOT4)*DTU +
     &                       TDOT3(IPAIR))*ONE6*DTU + TD2)*DTU + RI)*DTU
*
*       Ensure that regularized step is smaller than the c.m. step.
   50 IF (STEP(I1).GT.STEP(I)) THEN
          DTU = 0.5D0*DTU
          STEP(I1) = ((((0.2D0*TDOT5*DTU + 0.5D0*TDOT4)*DTU +
     &                       TDOT3(IPAIR))*ONE6*DTU + TD2)*DTU + RI)*DTU
          GO TO 50
      END IF
*
      DTAU(IPAIR) = DTU
      IF (IMOD.GT.1) THEN
          STEP(I1) = ZMOD*STEP(I1)
      END IF
*
*       Check diagnostic print option.
*     IF (KZ(12).GT.0) THEN
*         WRITE (6,52)  IPAIR, TIME, H(IPAIR), RI, DTAU(IPAIR), GI,
*    &                  STEP(I1), LIST(1,I1), IMOD
*  52     FORMAT (3X,'KS #',I6,2F10.4,2F10.5,1P2E10.2,2I4)
*     END IF
*
*       Employ special termination criterion in merger case.

C New: (Aug.1998, P.Kroupa) Bug-fix with Sverre Aarseth.
      IF (NAME(I).LT.0) THEN
          IF (GI.LT.GMAX) GO TO 60
*       Terminate if apocentre perturbation > 0.25 (R' > 0) or GI > 0.25.
          IF (HI.LT.0.0D0) THEN
              SEMI = -0.5*BODY(I)/HI
              ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
              A0 = SEMI*(1.0 + SQRT(ECC2))/RI
              GA = GI*A0*A0*A0
              IF (GA.GT.0.25.AND.RI.GT.SEMI) GO TO 90
              IF (GI.GT.0.1.AND.RI.GT.RMIN) GO TO 90
              IF (GI.GT.0.25) GO TO 90
          ELSE
              IF (GI.GT.GMAX) GO TO 90
          END IF
          GO TO 60
      END IF

C----------- Old:
c      IF (NAME(I).LT.0) THEN
c          IF (GI.LT.GMAX) GO TO 60
c*       Terminate if apocentre perturbation > 0.25 (R' > 0) or GI > 0.25.
c          IF (TD2.GT.0.0D0) THEN
c              SEMI = -0.5*BODY(I)/HI
c              ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
c              A0 = SEMI*(1.0 + SQRT(ECC2))/RI
c              GA = GI*A0*A0*A0
c              IF (GA.GT.0.25.AND.RI.GT.SEMI) GO TO 90
c          ELSE
c              IF (GI.GT.0.25) GO TO 90
c          END IF
c          GO TO 60
c      END IF
C------------ end old
*
*       Check standard termination criterion.
      IF (RI.GT.R0(IPAIR)) THEN
*
*       See whether termination can be delayed for intermediate energies.
          IF (CMSEP2*RI**2.LT.RS(I)**2) THEN
              A3 = RMIN*ABS(HI)*BODYIN
              GIMAX = A3*A3*A3
              IF (GI.LT.GIMAX) GO TO 60
          END IF
*
          IF (HI.LT.0.0) THEN
              A3 = RMIN*ABS(HI)*BODYIN
              GIMAX = A3*A3*A3
              IF (GI.LT.MIN(GIMAX,GMAX)) GO TO 60
          END IF
*
          IF (KZ(4).GT.0.AND.GI.GT.GPRINT(1)) THEN
              DGAM = GI - G0
              K = KZ(4)
              DO 55 L = 2,K
                  IF (GI.LT.GPRINT(L)) THEN
                      CALL EVOLVE(IPAIR,L-1)
                      GO TO 90
                  END IF
   55         CONTINUE
              CALL EVOLVE(IPAIR,K)
          END IF
          GO TO 90
      END IF
*
*       End integration cycle for hyperbolic motion.
   60 IF (HI.GE.0.0D0) THEN
          IF (RDOT*TD2.LT.0.0D0) THEN
*       Obtain pericentre by Mikkola's algorithm (GAMMA < 0.001).
              IF (GI.LT.0.001) THEN
                  CALL PERI(UI,UIDOT,RI,BODY(I1),BODY(I2),QPERI)
              ELSE
                  QPERI = RI
              END IF
              DMIN2 = MIN(DMIN2,QPERI)
*
C Added Aug.1998 b P.Kroupa:
* If BK(4)=1 then write out relevant data 
* into unit=15 (file=peri_hyperbol.dat) [see nbody6.F]
* for statistics of closest approaches & disk destruction during 
* hyperbolic encounters. (27.08.96, P.Kroupa in nb5)
* write out:  time in Myr, periastron distance in AU
*             names and masses (Msun) of involved stars
*             distance of centre of mass (!) from cluster centre in pc:
*
              if (BK(4).EQ.1) then
                 peri_pk=log10(QPERI*RBAR*206264.84)
                 cldist_pk=SQRT((x(1,I)-rdens(1))**2 + 
     +                       (x(2,I)-rdens(2))**2+
     +                       (x(3,I)-rdens(3))**2)*RBAR
                 write(15,'(2F8.3,2I5,2F8.3,F9.4,2F13.4)') 
     +            time*TSCALE_pk,peri_pk,name(I1),name(I2),
     +            body(I1)*ZMBAR,body(I2)*ZMBAR,cldist_pk,
C Msun/pc^3 in sphere containing 5 nearest stars, potential [(km/s)^2] at 
C position of centre of mass but neglecting the two interacting
C stars:
     +            RHODBL(I)*(ZMBAR/RBAR**3),PHIDBL(I)*(VSTAR_pk**2)
                 call flush(15)
C Also write the stuff into log-file:
                 write(6,*)
                 write(6,*)' Hyperbolic encounter: (see ksint.f)'
                 write(6,*)' time[Myr] log10(Rperi[AU]) star1 star2 ',
     +                  'mass1 mass2[Msun] Rcl[pc]'
                 write(6,'(2F8.3,2I5,2F8.3,F9.4)') 
     +            time*TSCALE_pk,peri_pk,name(I1),name(I2),
     +            body(I1)*ZMBAR,body(I2)*ZMBAR,cldist_pk
                 write(6,*)
                 call flush(6)
              end if
*
*
*       Check optional tidal interaction or stellar collision.
              IF (KZ(27).GT.0) THEN
                  IF (QPERI.LT.3.0*MAX(RADIUS(I1),RADIUS(I2))) THEN
                      IF (QPERI.LT.0.75*(RADIUS(I1) + RADIUS(I2))) THEN
*       Obtain KS variables at pericentre before merging into one body.
                          CALL KSPERI(IPAIR)
                          KSPAIR = IPAIR
                          CALL CMBODY(QPERI,2)
                      ELSE
                          CALL KSTIDE(IPAIR,QPERI)
                      END IF
                  END IF
              END IF
          END IF
          GO TO 100
      END IF
*
*       Check perturbation threshold (H < 0 & GAMMA > GMAX).
      IF (KZ(4).EQ.0.OR.G0.LT.GMAX) GO TO 70
*
      K = KZ(4)
      DO 65 L = 1,K
          IF ((G0 - GPRINT(L))*(GI - GPRINT(L)).LE.0.0) THEN
              IF (L.EQ.1) THEN
                  W1 = -0.5*BODY(I)/HI
                  W2 = W1*BODYIN
                  TK = TWOPI*W1*SQRT(W2)
              END IF
*
*       Estimate smallest permitted output interval at new level.
              DTCRIT = TK*ORBITS(L)
              IF (TIME - TLASTB(L).LT.DTCRIT) GO TO 70
              DGAM = GI - G0
              CALL EVOLVE(IPAIR,L)
              GO TO 70
          END IF
   65 CONTINUE
*
*       Check for partial reflection during approach (NB! only IMOD = 1).
*  70 IF (GI.LT.GMIN.AND.TD2.LT.0.0D0) THEN
*       Skip apocentre position itself.
*         IF (RDOT.LT.0.0D0.AND.IMOD.EQ.1) THEN
*             IF (KZ(25).GT.0) CALL FREEZE(IPAIR)
*             GO TO 100
*         END IF
*     END IF
*
*       Determine new perturbers for binary at apocentre turning point.
   70 IF (RDOT*TD2.GE.0.0D0) GO TO 100
      SEMI = -0.5D0*BODY(I)/HI
*
*       Check minimum two-body separation just after pericentre.
      IF (RDOT.LT.0.0D0) THEN
*       Obtain pericentre by Mikkola's algorithm (GAMMA < 0.001).
          IF (GI.LT.0.001) THEN
              CALL PERI(UI,UIDOT,RI,BODY(I1),BODY(I2),QPERI)
          ELSE
              QPERI = RI
          END IF
          DMIN2 = MIN(DMIN2,QPERI)
*       Check optional tidal interaction or stellar collision.
          IF (KZ(27).GT.0) THEN
              IF (QPERI.LT.2.7*MAX(RADIUS(I1),RADIUS(I2))) THEN
                  IF (QPERI.LT.0.75*(RADIUS(I1) + RADIUS(I2))) THEN
*       Obtain KS variables at pericentre before merging into one body.
                      CALL KSPERI(IPAIR)
                      KSPAIR = IPAIR
                      CALL CMBODY(QPERI,2)
*       Do not evolve synchronous orbit further.
                  ELSE IF (SEMI.GT.RSYNC) THEN
                      CALL KSTIDE(IPAIR,QPERI)
                  END IF
              END IF
          END IF
          GO TO 100
      END IF
*
*       Save maximum separation of persistent binary.
      RMAX = MAX(RMAX,RI)
*
*       Check binary reference radius or merger termination.
      IF (NAME(I).GT.0) THEN
*       Update termination length scale of hard binary if enough perturbers.
          EB = BODY(I1)*BODY(I2)*HI*BODYIN
          IF (EB.LT.-0.5*BODYM*ECLOSE.AND.RI.LT.0.02*RS(I)) THEN
              R0(IPAIR) = MAX(RMIN,2.0*SEMI)
          END IF
      ELSE
*       Check pericentre stability criterion of merged binary (NAME < 0).
          ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          RP = SEMI*(1.0 - SQRT(ECC2))
          IF (RP*(1.0 - 2.0*GI).LT.R0(IPAIR)) GO TO 90
      WRITE (7,777)  LIST(1,I),RI,SEMI,SQRT(ECC2),2.0*RS(I)
  777 FORMAT (' APO:   NB RI A E 2*RS ',I4,2F8.4,2F7.3)
      CALL FLUSH(7)
*       Terminate merger if maximum perturber range > 2*RS.
          IF (CMSEP2*RI**2.GT.4.0*RS(I)**2) GO TO 90
      END IF
*
      IF (KZ(26).GT.0) THEN
*       See whether KS slow-down procedure should be (re)checked.
          IF (KZ(27).GT.0) THEN
              ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
              QPERI = SEMI*(1.0 - SQRT(ECC2))
              RT = 3.0*MAX(RADIUS(I1),RADIUS(I2))
*       Skip modification if predicted pericentre < 5 tidal capture radii.
              IF (QPERI.LT.5.0*RT) GO TO 80
          END IF
          KMOD = RANGE*GMIN/GI
          IF (KMOD.GT.1.OR.IMOD.GT.1) THEN
              CALL KSMOD(IPAIR,KMOD)
              IF (KMOD.LT.0) GO TO 100
          END IF
      END IF
*
*       Set approximate value of next period with perturbation included.
   80 TK = TWOPI*SEMI*SQRT(SEMI*BODYIN)*(1.0 + GI)
      IF (IMOD.GT.1) THEN
          TK = ZMOD*TK
      END IF
*
*       Use old perturber list if next apocentre is before the c.m. step.
      IF (TIME + TK.LT.T0(I) + STEP(I)) THEN
          IF (KZ(26).EQ.0.OR.KZ(27).EQ.0) THEN
              GO TO 100
          ELSE
              IF (SEMI.GT.RSYNC.AND.QPERI.LT.RT) GO TO 100
          END IF
      END IF
*
*       Select new perturbers (adopt unperturbed period if none found).
      CALL KSLIST(IPAIR)
*
*       Check optional search criterion for multiple encounter or merger.
      IF (KZ(15).NE.0.AND.STEP(I).LT.DTMIN) THEN
          CALL IMPACT(I)
      END IF
      GO TO 100
*
*       Terminate regularization of current pair (IPAIR set in KSPAIR).
   90 KSPAIR = IPAIR
*       Set indicator for calling KSTERM in MAIN (permits phase overlay).
      IPHASE = 2
*       Check for rare case of merged binary component.
      IF (NAME(I).LT.0) IPHASE = 7
*
  100 RETURN
*
      END
