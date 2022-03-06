      SUBROUTINE KSINT(I1)
*
*
*       Regularized integration.
*       ------------------------
*
* ---------------------------------------------------
* Modified for PN evolution by J.M.B. Downing 02.2008
* ---------------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX)
      COMMON/SLOW0/  RANGE,ISLOW(10)
      COMMON/GAMDOT/  DGAM
      REAL*8  UI(4),UIDOT(4),XI(6),VI(6),FP(6),FD(6)
*
* PN variables
      REAL*8 FRELP1(6),FRELP2(6),FRELP2_5(6),FRELD1(6),FRELD2(6)
      REAL*8 FRELD2_5(6),VWHOLE,CMVI(3),VCM(6),DTUREL,DTUCLASS
      INTEGER RELORNOT
*
      LOGICAL IQ
*
*
*       Set second component, pair index & c.m. index.
      I2 = I1 + 1
      IPAIR = KVEC(I1)
      I = N + IPAIR
*
*       Define perturber membership & inverse c.m. mass.
      NNB0 = LIST(1,I1)
      BODYIN = 1.0/BODY(I)
*
* ------------------------------------------------------------
* Calculate VOVERC to determine if the binariy is relativistic
* ------------------------------------------------------------
*
      IF (KZ(41).LT.5) THEN
         IF ((VI(1)*VI(1)+VI(2)*VI(2)+VI(3)*VI(3)+VI(4)*VI(4)
     &        +VI(5)*VI(5)+VI(6)*VI(6)).GT.0.0) THEN
*
            CMVI(1) = (BODY(I1)*VI(1)+BODY(I2)*VI(4))
     &           /(BODY(I1)+BODY(I2))
            CMVI(2) = (BODY(I1)*VI(2)+BODY(I2)*VI(5))
     &           /(BODY(I1)+BODY(I2))
            CMVI(3) = (BODY(I1)*VI(3)+BODY(I2)*VI(6))
     &           /(BODY(I1)+BODY(I2))
*
            VCM(1) = VI(1)-CMVI(1)
            VCM(2) = VI(2)-CMVI(2)
            VCM(3) = VI(3)-CMVI(3)
            VCM(4) = VI(4)-CMVI(1)
            VCM(5) = VI(5)-CMVI(2)
            VCM(6) = VI(6)-CMVI(3)
*
            VOVERC = SQRT(MAX(VCM(1)*VCM(1)+VCM(2)*VCM(2)+VCM(3)*VCM(3),
     &           VCM(4)*VCM(4)+VCM(5)*VCM(5)+VCM(6)*VCM(6)))/CLIGHT
         END IF
*
* -------------------------------------------
* Check if the binary requires PN integration
* -------------------------------------------
*
* Relativistic binaries are always perturbed so call to
* unpert.f is skipped
* See define.f for the meanings of KZ(J)
*         IF (KZ(41).LT.5) THEN
*         IF (((KZ(42).EQ.1).OR.((KZ(42).EQ.2).AND.((NAME(I1).EQ.1
*     &        .AND.NAME(I2).EQ.2).OR.(NAME(I1).EQ.2.AND.NAME(I2).EQ.1)))
*     &        .OR.(KZ(42).EQ.3.AND.(NAME(I1).EQ.1.OR.NAME(I2).EQ.1)))
         IF ((KZ(43).EQ.0).OR.(KZ(43).EQ.1.AND.(VOVERC.GT.BETAREL)))
     &        GO TO 3333
      END IF
*
* If not relativistic, check for further unperturbed motion or dissipation 
* at pericentre.
      IF (NNB0.EQ.0.AND.H(IPAIR).LT.0.0) THEN
          CALL UNPERT(IPAIR)
          GO TO 100
      END IF
*
* If relativistic continue without checking unpertubed motion
 3333 CONTINUE
*
*       Perform KS prediction of U, UDOT & H.
      CALL KSPRED(IPAIR,I1,I,BODYIN,UI,UIDOT,XI,VI)
*
*       Obtain the perturbing force & derivative.
      CALL KSPERT(I1,NNB0,XI,VI,FP,FD)
*
* -----------------------------------------------
* Determine and apply PN correction if neccessary
* -----------------------------------------------
*
* Make another determination if binary gets relativistic treatment
* (probably not neccssary if a flag were introduced...)
* Again see define.f
*
      IF (KZ(41).LT.5) THEN
*         IF ((KZ(42).EQ.1).OR.((KZ(42).EQ.2).AND.((NAME(I1).EQ.1
*     &       .AND.NAME(I2).EQ.2).OR.(NAME(I1).EQ.2.AND.NAME(I2).EQ.1)))
*     &       .OR.(KZ(42).EQ.3.AND.(NAME(I1).EQ.1.OR.NAME(I2).EQ.1)))THEN
*
         IF((KZ(43).EQ.0).OR.(KZ(43).EQ.1.AND.(VOVERC.GT.BETAREL)))
     &        THEN
*
* Actual calculation of PN corrections (not yet applied)
            CALL KSREL(IPAIR,I1,XI,VI,FRELP1,FRELP2,FRELP2_5,FRELD1,
     &           FRELD2,FRELD2_5,VWHOLE)
*
* ---------------------------------------
* Check if there is a relativistic merger
* ---------------------------------------
*
            RIJ2 = (XI(1)-XI(4))*(XI(1)-XI(4))
     &         + (XI(2)-XI(5))*(XI(2)-XI(5))+(XI(3)-XI(6))*(XI(3)-XI(6))
            RIJ = SQRT(RIJ2)
*
* Critirion is crude, simple Schwarzchild radius overlap.
            IF (RIJ.LT.(6.0*(BODY(I1)+BODY(I2)))/(CLIGHT*CLIGHT)) THEN
*
               WRITE(*,556) TIME+TOFF,NAME(I1),NAME(I2),BODY(I1)+
     &              BODY(I2),RIJ, 2.0*(BODY(I1)+BODY(I2))/CLIGHT/CLIGHT,
     &              GAMMA(IPAIR),GAMMAREL1(IPAIR),GAMMAREL2(IPAIR),
     &              GAMMAREL2_5(IPAIR),VWHOLE,MAX(BODY(I1)/BODY(I2),
     &              BODY(I2)/BODY(I1)),KSTAR(I1),KSTAR(I2)
 556        FORMAT ('Relativistic merger -->  TIME = ',F12.5,'  I = ',
     &           I6,'  J = ',I6,'  TOT_MASS = ',e15.6,'  RIJ = ',e15.6,
     &           '  R_SCHWARTZ = ',e15.6,' GAMMA = ',e15.6,
     &           '  GAMMA_REL_1PN = ',e15.6,
     &           '  GAMMA_REL_2PN = ',e15.6,
     &           '  GAMMA_REL_2.5PN = ',e15.6,
     &           '  MAX{|V1|,|V2|} = ',e15.6,'  MASS RATIO = ',e15.6,
     &           '  STELLAR_TYPE_I = ',I2,'  STELLAR_TYPE_J = ',I2)
*
               CALL FLUSH(6)
*
* Obtain KS variables at pericenter before merging into one body
               CALL KSPERI(IPAIR)
               KSPAIR = IPAIR
*
* Perform relativistic merger
               CALL CMGREL(QPERI,2)
*
* Since binary is merged there is no application of PN corrections
* and KS regularization is terminated.
               GO TO 100
*
            END IF
         END IF
      END IF
*
*       Save old radial velocity & relative perturbation and set new GAMMA.
      RDOT = TDOT2(IPAIR)
      G0 = GAMMA(IPAIR)
      GI = SQRT(FP(1)**2 + FP(2)**2 + FP(3)**2)*R(IPAIR)**2*BODYIN
      GAMMA(IPAIR) = GI
*
* ------------------------------------------
* Apply PN corrections to regularized motion
* ------------------------------------------
*
      IF (KZ(41).LT.5) THEN
         IF((KZ(43).EQ.0).OR.(KZ(43).EQ.1.AND.(VOVERC.GT.BETAREL)))
     &        THEN
            DO 337 KK = 1,3
               FP(KK) = FP(KK) + FRELP1(KK) + FRELP2(KK) + FRELP2_5(KK)
               FD(KK) = FD(KK) + FRELD1(KK) + FRELD2(KK) + FRELD2_5(KK)
 337        CONTINUE
         END IF
      END IF
*
*       Apply the Hermite corrector.
      CALL KSCORR(IPAIR,UI,UIDOT,FP,FD,TD2,TDOT4,TDOT5,TDOT6)
*
*       Increase regularization time-step counter and update the time.
      NSTEPU = NSTEPU + 1
      T0(I1) = TIME
*
*       Check for early return during termination (called from KSTERM).
      IF (IPHASE.NE.0) GO TO 100
*
*       Define useful scalars.
      RI = R(IPAIR)
      HI = H(IPAIR)
*
*       Initialize termination indicator and check for large perturbation.
      IQ = .FALSE.
      IF (GI.GT.0.5) GO TO 2
      IF (GI.LT.0.2) THEN
          JCOMP = 0
          GO TO 20
      END IF
      CALL FLYBY(I,ITERM)
      IF (ITERM.EQ.1.AND.KZ(15).GT.0) THEN
          CALL IMPACT(I)
          IF (IPHASE.GT.0) GO TO 100
      ELSE IF (ITERM.EQ.2) THEN
          IQ = .TRUE.
          GO TO 20
      ELSE
          GO TO 20
      END IF
*
*       Find the dominant body for large perturbations.
    2 S = 2.0*STEP(I)
      FMAX = BODY(I)/RI**2
*       Initialize JCOMP for optional diagnostics in KSTERM.
      JCOMP = 0
      DO 10 L = 2,NNB0+1
          J = LIST(L,I1)
*       Only search bodies within twice the c.m. time-step.
          IF (STEP(J).GT.S) GO TO 10
*       Compare strong perturber and either component with current pair.
          DO 5 K = I1,I2
              RIJ2 = (X(1,J) - X(1,K))**2 + (X(2,J) - X(2,K))**2 +
     &                                      (X(3,J) - X(3,K))**2
              IF (BODY(J) + BODY(K).GT.RIJ2*FMAX) JCOMP = J
    5     CONTINUE
   10 CONTINUE
*
*       Set termination if strong perturber <= N forms dominant pair.
      IF (JCOMP.GT.0.OR.GI.GT.1.0) THEN
*       Check optional binary search.
          IF (KZ(4).GT.0) THEN
              DGAM = GI - G0
              K = KZ(4)
              CALL EVOLVE(IPAIR,K)
          END IF
          IF (JCOMP.LE.N.OR.GI.GT.2.0) IQ = .TRUE.
      END IF
*
*       Check termination of hyperbolic encounter (R > R0 or R > RMIN).
   20 IF (HI.GT.0.0D0) THEN
          IF ((RI.GT.R0(IPAIR).AND.GI.GT.GMAX).OR.RI.GT.RMIN.OR.
     &        (GI.GT.0.5.AND.TD2.GT.0.0)) THEN
*       Skip termination delay in case of velocity kick (cf. routine KSTERM).
              IF (HI.LT.100.0.OR.GI.GT.0.5) THEN
                  IQ = .TRUE.
              END IF
          END IF
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
*       Modify maximum square step for soft binaries & weak hyperbolic pairs.
          IF (HI.LT.0.0D0) THEN
*       Consider case of hard binary with massive components or merger.
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
*       Form new regularized step.
      DTU = ETAU*W2
      DTU = MIN(1.2*DTAU(IPAIR),DTU)
*
* Relativistic timestep critirion differnet, PN not a tidal force
      IF (KZ(41).LT.5) THEN
         IF((KZ(43).EQ.0).OR.(KZ(43).EQ.1.AND.(VOVERC.GT.BETAREL)))
     &        THEN
            DTUREL = ETAU/(1.0+DTUPN*(SQRT((BODY(I1)+BODY(I2))/RIJ)
     &           /CLIGHT))
*
* Record if reletivistic timestep is used (for diagnositics only)
            IF (DTUREL.LT.DTU) THEN
               RELORNOT = 1
               DTUCLASS = DTU
            ELSE IF (DTUREL.GE.DTU) THEN
               RELORNOT = 0
               DTUCLASS = DTU
            END IF
*
* Determine if classical or PN dynamics dominate the timestep
            DTU = MIN(DTUREL,DTU)
*
* Output timestep data for diagnostic purposes
            IF (KZ(47).GT.0.AND.RFLAG(IPAIR).GT.0) THEN
               WRITE (36,234) TIME+TOFF,DTU,DTUCLASS,DTUREL,RELORNOT
 234           FORMAT (4E15.7,I3)
            END IF
         END IF
      END IF
*
*       Reset reference energy and generate new Stumpff coefficients.
      H0(IPAIR) = H(IPAIR)
   30 Z = -0.5D0*H(IPAIR)*DTU**2
      CALL STUMPF(IPAIR,Z)
      Z5 = SF(6,IPAIR)
      Z6 = SF(7,IPAIR)
      DT12 = ONE12*DTU*Z6
*
*       Convert to physical time units modified by Stumpff coefficients.
      STEP(I1) = (((((TDOT6*DT12 + TDOT5*Z5)*0.2*DTU + 0.5D0*TDOT4)*DTU
     &                     + TDOT3(IPAIR))*ONE6*DTU + TD2)*DTU + RI)*DTU
*
*       Ensure that regularized step is smaller than the c.m. step.
      IF (STEP(I1).GT.STEP(I)) THEN
          DTU = 0.5D0*DTU
          GO TO 30
      END IF
      DTAU(IPAIR) = DTU
*
* Set timestep for relativistic energy calculation (see ksrel.f)
      RELTSTEP(IPAIR) = STEP(I1)
*
*       See whether the KS slow-down procedure is activated. 
      IMOD = KSLOW(IPAIR)
      IF (IMOD.GT.1) THEN
          ZMOD = FLOAT(ISLOW(IMOD))
          STEP(I1) = ZMOD*STEP(I1)
      END IF
*
*       Check diagnostic print option.
      IF (KZ(10).GE.3) THEN
          if(rank.eq.0)
     &    WRITE (6,40)  IPAIR, TIME+TOFF, H(IPAIR), RI, DTAU(IPAIR),
     &                  GI, STEP(I1), LIST(1,I1), IMOD
   40     FORMAT (3X,'KS MOTION',I6,2F10.4,2F10.5,1P,2E10.2,2I4)
      END IF
*
*       Employ special termination criterion in merger case.
      IF (NAME(I).LT.0) THEN
*       Terminate if apocentre perturbation > 0.25 (R' > 0) or GI > 0.25.
          IF (HI.LT.0.0) THEN
              SEMI = -0.5*BODY(I)/HI
              ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
              A0 = SEMI*(1.0 + SQRT(ECC2))/RI
              GA = GI*A0*A0*A0
              IF (GA.GT.0.25.AND.RI.GT.SEMI) IQ = .TRUE.
              IF (GI.GT.0.1.AND.RI.GT.RMIN) IQ = .TRUE.
              IF (GI.GT.0.25) IQ = .TRUE.
          ELSE
              IF (GI.GT.GMAX) IQ = .TRUE.
          END IF
          IF (.NOT.IQ) GO TO 70
      END IF
*
*       Delay termination until end of block for large perturbation.
      IF (IQ) THEN
          DTR = TBLOCK - TIME
*         if(rank.eq.0)
*    &    WRITE (6,45)  IPAIR, TIME+TOFF, GI, RI, DTR, STEP(I1)
*  45     FORMAT (' TERM TEST    KS T G R DTR DT  ',
*    &                           I4,F10.4,F7.3,1P,E10.2,2E9.1)
          IF (DTR.LT.STEP(I1)) GO TO 90
      END IF
*
*       Check standard termination criterion (suppress on IQ = .true.).
      IF (RI.GT.R0(IPAIR).AND..NOT.IQ) THEN
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
*       Check updating of R0 for newly hardened binary orbit.
              IF (HI.LT.-ECLOSE) THEN
                  SEMI = -0.5*BODY(I)/HI
                  R0(IPAIR) = MAX(RMIN,2.0D0*SEMI)
                  GO TO 70
              END IF
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
              IF (rank.eq.0.and.BK(4).EQ.1) THEN
C Added Aug.1998 b P.Kroupa:
* If BK(4)=1 then write out relevant data
* into unit=15 (file=per.15) [see file_init.F]
* for statistics of closest approaches & disk destruction during
* hyperbolic encounters. (27.08.96, P.Kroupa in nb5)
* write out:  time in Myr, periastron distance in AU
*             names and masses (Msun) of involved stars
*             distance of centre of mass (!) from cluster centre in pc:
*
                 peri_pk=log10(QPERI*RBAR*206264.84)
                 cldist_pk=SQRT((x(1,I)-rdens(1))**2 +
     +                       (x(2,I)-rdens(2))**2+
     +                       (x(3,I)-rdens(3))**2)*RBAR
#ifdef PARALLEL
      if(rank.eq.0)then
#endif
                 write(15,'(2F8.3,2I5,2F8.3,F9.4,F13.4)')
     +            time*tscale,peri_pk,name(I1),name(I2),
     +            body(I1)*ZMBAR,body(I2)*ZMBAR,cldist_pk
                 call flush(15)
#ifdef PARALLEL
      end if
#endif
C Also write the stuff into log-file:
*                write(6,*)
*                write(6,*)' Hyperbolic encounter: (see ksint.f)'
*                write(6,*)' time[Myr] log10(Rperi[AU]) star1 star2 ',
*    +                  'mass1 mass2[Msun] Rcl[pc]'
*                write(6,'(2F8.3,2I5,2F8.3,F9.4)')
*    +            time*tscale,peri_pk,name(I1),name(I2),
*    +            body(I1)*ZMBAR,body(I2)*ZMBAR,cldist_pk
*                write(6,*)
*                call flush(6)
              end if
*
*       Check optional tidal interaction or stellar collision.
              IF (KZ(27).GT.0) THEN
                  IF (QPERI.LT.4.0*MAX(RADIUS(I1),RADIUS(I2))) THEN
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
*
*       Check optional tidal interaction or stellar collision (skip merger).
          IF (KZ(27).GT.0) THEN
              ECC = 1.0 - QPERI/SEMI
              ZF = 4.0
              IF (ECC.GT.0.95) ZF = 50.0
              IF (QPERI.LT.ZF*MAX(RADIUS(I1),RADIUS(I2)).AND.
     &                                           NAME(I).GT.0) THEN
                  IF (QPERI.LT.0.75*(RADIUS(I1) + RADIUS(I2))) THEN
*       Obtain KS variables at pericentre before merging into one body.
                      CALL KSPERI(IPAIR)
                      KSPAIR = IPAIR
                      CALL CMBODY(QPERI,2)
*       Do not evolve synchronous orbit further.
                  ELSE IF (KSTAR(I).NE.20) THEN
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
          IF (EB.LT.EBH.AND.RI.LT.0.02*RS(I)) THEN
              R0(IPAIR) = MAX(RMIN,2.0*SEMI)
          ELSE
              R0(IPAIR) = MAX(RMIN,0.03*RS(I))
          END IF
      ELSE 
*       Check pericentre stability criterion of merged binary (NAME < 0).
          ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          ECC1 = SQRT(ECC2)
          RP = SEMI*(1.0 - ECC1)
          IF (RP*(1.0 - GI).LT.R0(IPAIR)) THEN
              Q = BODY(I2)/BODY(I)
              XFAC = (1.0 + Q)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
              DO 72 K = 1,NMERGE
                  IF (NAMEM(K).EQ.NAME(I)) IM = K
   72         CONTINUE
              AIN = -0.5*BODY(I1)/HM(IM)
              PCRIT = 2.8*XFAC**0.4*AIN
              WRITE (6,73)  GI, ECC1, RP, R0(IPAIR), PCRIT, AIN
   73         FORMAT (' BASIC INSTAB    GI E1 RP R0 PCR AIN ',
     &                                  2F7.3,1P,4E10.2)
              GO TO 90
          END IF
*       Include safety check for outer eccentricity change (GI > 0.001).
          IF (GI.GT.0.001) THEN
              Q = BODY(I2)/BODY(I)
              XFAC = (1.0 + Q)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
              DO 75 K = 1,NMERGE
                  IF (NAMEM(K).EQ.NAME(I)) IM = K
   75         CONTINUE
              AIN = -0.5*BODY(I1)/HM(IM)
              PCRIT = 2.8*XFAC**0.4*AIN
*       Adopt 60 % criterion to allow for marginally stable systems.
              IF (RP.LT.0.6*PCRIT) THEN
              WRITE (6,500) GI, ECC1, RP, R0(IPAIR), PCRIT, AIN
  500         FORMAT (' PERT INSTAB!    GI E1 RP R0 PCR AIN ',
     &                                  2F7.3,1P,4E10.2)
              GO TO 90
              END IF
          END IF
*       Terminate merger if maximum perturber range > 2*RS & GI > GMAX.
          IF (CMSEP2*RI**2.GT.4.0*RS(I)**2.AND.GI.GT.GMAX) GO TO 90
      END IF
*
*       See whether KS slow-down procedure should be (re)-checked.
      IF (KZ(26).GT.0) THEN
*       Include case of tidal dissipation.
          IF (KZ(27).GT.0) THEN
              ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
              QPERI = SEMI*(1.0 - SQRT(ECC2))
              RT = 4.0*MAX(RADIUS(I1),RADIUS(I2))
*       Skip modification if predicted pericentre < 5 tidal capture radii.
              IF (QPERI.LT.5.0*RT.AND.GI.LT.RANGE*GMIN) GO TO 80
          END IF
          KMOD = RANGE*GMIN/MAX(GI,1.0D-10)
          IF (KMOD.GT.1.OR.IMOD.GT.1) THEN
              CALL KSMOD(IPAIR,KMOD)
              IF (KMOD.LT.0) GO TO 100
          END IF
      END IF
*
*       Set approximate value of next period with perturbation included.
   80 TK = TWOPI*SEMI*SQRT(SEMI*BODYIN)*(1.0D0 + GI)
      IF (IMOD.GT.1) THEN
          TK = ZMOD*TK
      END IF
*
*       Use old perturber list if next apocentre is before the c.m. step.
      IF (TIME + TK.LT.T0(I) + STEP(I)) THEN
          IF (KZ(26).EQ.0.OR.KZ(27).EQ.0) THEN
              GO TO 100
          ELSE
              IF (QPERI.LT.RT.AND.KSTAR(I).NE.20) GO TO 100
          END IF
      END IF
*
*       Select new perturbers (adopt unperturbed period if none found).
      CALL KSLIST(IPAIR)
*
*       Perform rectification.
      CALL KSRECT(IPAIR)
*
*       Check optional search criterion for multiple encounter or merger.
      IF (KZ(15).GT.0.AND.STEP(I).LT.DTMIN) THEN
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
