      SUBROUTINE KSTIDE(IPAIR,QPERI)
*
*
*       Tidal interaction of KS pair.
*       -----------------------------
*
      INCLUDE 'common6.h'
      REAL  DE(2)
      DATA  ECCMIN  /0.01/
*
*
*       Skip procedure if both stars are highly evolved.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      IF (RADIUS(I1) + RADIUS(I2).LE.0.0D0) GO TO 50
*
*       Obtain kinetic energy loss due to tidal interaction (DE > 0 here).
      ICASE = KZ(27)
      CALL TIDES(QPERI,BODY(I1),BODY(I2),RADIUS(I1),RADIUS(I2),ICASE,DE)
*
*       Set c.m. index & reduced mass.
      I = N + IPAIR
      ZMU = BODY(I1)*BODY(I2)/BODY(I)
*
*       Skip orbit correction if energy loss < 0.01%.
      IF (DE(1) + DE(2).LT.0.0001*ZMU*ABS(H(IPAIR))) GO TO 50
*
*       Determine pericentre variables U & UDOT by backwards reflection.
      CALL KSPERI(IPAIR)
*
*       Form semi-major axis & eccentricity (TDOT2 = 0 at pericentre).
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      ECC = 1.0 - R(IPAIR)/SEMI
      PERI = SEMI*(1.0D0 - ECC)
      AM0 = SEMI*(1.0D0 - ECC**2)
*
*       Calculate the new eccentricity from angular momentum conservation.
      DH = -(DE(1) + DE(2))/ZMU
      ECC2 = ECC**2 + 2.0D0*AM0*DH/BODY(I)
      IF (ECC2.GT.ECCMIN**2) THEN
          ECC1 = SQRT(ECC2)
      ELSE
          ECC1 = ECCMIN
      END IF
*
*       Set new energy and determine semi-major axis & pericentre.
      HI = H(IPAIR)
      H(IPAIR) = H(IPAIR) + DH
      SEMI1 = -0.5D0*BODY(I)/H(IPAIR)
      PERI1 = SEMI1*(1.0D0 - ECC1)
*
*       Form KS coordinate scaling factor from pericentre ratio.
      C1 = SQRT(PERI1/PERI)
*       Specify KS velocity scaling from angular momentum conservation.
      C2 = 1.0/C1
*       See whether circular orbit condition applies.
      AM = SEMI1*(1.0D0 - ECC1**2)
      IF (ECC1.LE.ECCMIN) C2 = SQRT(AM/AM0)/C1
*
*       Transform KS variables to yield the prescribed elements.
      R(IPAIR) = 0.0D0
      DO 10 K = 1,4
          U(K,IPAIR) = C1*U(K,IPAIR)
          UDOT(K,IPAIR) = C2*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
   10 CONTINUE
*
*       Form new perturber list after significant energy loss.
      IF (ABS(SEMI1/SEMI).LT.0.5) THEN
          CALL KSLIST(IPAIR)
      END IF
*
*       Obtain physical variables at pericentre and form new KS polynomials.
      T0(2*IPAIR-1) = TIME
      CALL RESOLV(IPAIR,1)
      IMOD = KSLOW(IPAIR)
      CALL KSPOLY(IPAIR,IMOD)
*
      IF (ECC.GT.0.99) THEN
          WRITE (6,20)  NAME(I1), NAME(I2), SEMI1, ECC, ECC1, HI, QPERI
   20     FORMAT (' NEW KSTIDE   NAM AF E0 EF HI QP ',
     &                                 2I5,1PE10.2,0P2F8.3,F9.1,1PE10.2)
          CALL FLUSH(6)
      END IF
*
*       Check for hierarchical configuration with eccentric inner binary.
      IF (ECC.GT.0.99.AND.HI.LT.0.0) THEN
          NP1 = LIST(1,I1) + 1
          DO 30 L = 2,NP1
              J = LIST(L,I1)
              RIJ2 = 0.0
              VIJ2 = 0.0
              RDOT = 0.0
              DO 25 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
                  VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
                  RDOT = (X(K,I) - X(K,J))*(XDOT(K,I) - XDOT(K,J))
   25         CONTINUE
              RIP = SQRT(RIJ2)
              A1 = 2.0/RIP - VIJ2/(BODY(I) + BODY(J))
              A1 = 1.0/A1
              IF (1.0/A1.GT.0.5/RMIN) THEN
                  ECC2 = (1.0 - RIP/A1)**2 +
     &                                  RDOT**2/(A1*(BODY(I) + BODY(J)))
                  RP = A1*(1.0 - SQRT(ECC2))
                  RA = SEMI*(1.0 + ECC)
                  SR = RP/RA
                  WRITE (6,28)  IPAIR, H(IPAIR), SEMI, A1, RP,
     &                          SQRT(ECC2), SR
   28             FORMAT (' HIERARCHY:   IPAIR H A0 A1 E1 RP SR ',I4,
     &                                         F7.0,1P3E9.1,0PF6.2,F6.1)
                  CALL FLUSH(6)
              END IF
   30     CONTINUE
      END IF
*
*       Do not allow R' < 0 for pericentre procedure in routine KSINT.
      IF (TDOT2(IPAIR).LT.0.0D0) THEN
          TDOT2(IPAIR) = 0.0D0
      END IF
*
*       Increase event counter and update total energy loss.
      NCOUNT(27) = NCOUNT(27) + 1
      NDISS = NDISS + 1
      ECOLL = ECOLL + (DE(1) + DE(2))
*
*       Count any hyperbolic captures.
      IF (SEMI.LT.0.0.AND.SEMI1.GT.0.0) THEN
          NTIDE = NTIDE + 1
      END IF
*
*       Record diagnostics for new synchronous orbit.
      IF (SEMI.GT.RSYNC.AND.SEMI1.LT.RSYNC) THEN
          NSYNC = NSYNC + 1
          ESYNC = ESYNC + ZMU*H(IPAIR)
          WRITE (6,40)  NAME(I1), NAME(I2), SEMI1, ECC, ECC1, HI, QPERI
   40     FORMAT (' END KSTIDE    NAM AF E0 EF HI QP',
     &                                 2I5,1PE10.2,0P2F8.3,F8.0,1PE10.2)
          CALL FLUSH(6)
      END IF
*
   50 RETURN
*
      END
