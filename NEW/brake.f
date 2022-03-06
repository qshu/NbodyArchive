      SUBROUTINE BRAKE(IPAIR,DSEP)
*
*
*       Orbital changes (GR, mass loss and/or tides).
*       ---------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  M1,M2
*
*
*       Set indices of KS components & c.m.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      I = N + IPAIR
      DSEP = DSEP/SU 
*
*       Set mass and radius of each star in solar units.
      M1 = BODY(I1)*ZMBAR
      M2 = BODY(I2)*ZMBAR
      R1 = RADIUS(I1)*SU
      R2 = RADIUS(I2)*SU
*
*       Obtain semi-major axis and eccentricity.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
*
      SEMI1 = SEMI - DSEP 
*
*       Check collision/coalescence (IQCOLL=-2 skips ECOLL updating in COAL).
      RP = SEMI1*SU*(1.D0 - SQRT(ECC2))
      IF(KSTAR(I).LE.0.AND.(RP.LE.(R1+R2)))THEN
         CALL KSPERI(IPAIR)
         IQCOLL = -2
         CALL CMBODY(R(IPAIR),2)
         GOTO 90
      ENDIF
      IF(DSEP.EQ.0.D0) GOTO 90
*
*       Include safety test on new semi-major axis.
      RCHCK = MIN(RADIUS(I1),RADIUS(I2)) 
      IF(SEMI1.LT.RCHCK) SEMI1 = RCHCK 
*
*       Transform to pericentre (R > A & unperturbed).
*     IF(R(IPAIR).GT.SEMI.AND.LIST(1,I1).EQ.0)THEN
*        CALL KSAPO(IPAIR)
*     ENDIF
*
*       Form square of regularized velocity.
      V20 = 0.0
      DO 10 K = 1,4
         V20 = V20 + UDOT(K,IPAIR)**2
 10   CONTINUE
*
*       Update binding energy and collision energy.
      HI = H(IPAIR)
      H(IPAIR) = -0.5*BODY(I)/SEMI1
      ZMU = BODY(I1)*BODY(I2)/BODY(I) 
      ECOLL = ECOLL + ZMU*(HI - H(IPAIR))
*
*       Specify KS coordinate & velocity scaling factors at arbitrary point.
      C2 = SQRT(SEMI1/SEMI)
*     V2 = 0.5*(BODY(I) + H(IPAIR)*SEMI1*(1.0 - SQRT(ECC2)))
      V2 = 0.5*(BODY(I) + H(IPAIR)*R(IPAIR)*(SEMI1/SEMI))
      C1 = SQRT(V2/V20)
*
*       Re-scale KS variables to new energy with constant eccentricity.
      R(IPAIR) = 0.0D0
*     TDOT2(IPAIR) = 0.0D0
      DO 20 K = 1,4
         U(K,IPAIR) = C2*U(K,IPAIR)
         UDOT(K,IPAIR) = C1*UDOT(K,IPAIR)
         U0(K,IPAIR) = U(K,IPAIR)
         R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
*        TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0*U(K,IPAIR)*UDOT(K,IPAIR)
 20   CONTINUE
*
*       Transform back to apocentre for standard unperturbed motion.
*     IF (LIST(1,I1).EQ.0) THEN
*        CALL KSAPO(IPAIR)
*     END IF
*
*       Include new initialization for perturbed orbit.
      IF (LIST(1,I1).GT.0) THEN
         IMOD = KSLOW(IPAIR)
         CALL KSPOLY(IPAIR,IMOD)
      ENDIF
*
*       Include some diagnostic output.
      IF(KSTAR(I).EQ.13)THEN
         NDIAG = NDIAG + 1
         IF(NDIAG.LT.100.OR.MOD(NDIAG,100).EQ.100)THEN
            RCOLL = RADIUS(I1) + RADIUS(I2)
            if(rank.eq.0)
     &      WRITE (6,25)  TTOT, IPAIR, M1, M2, R1, R2, R(IPAIR),
     &                    RCOLL
 25         FORMAT (' BRAKE    T KS M12 R12 R RCOLL ',
     &                         F10.4,I4,2F6.2,2F7.3,1P,2E10.2)
         ENDIF
      ENDIF
*
 90   RETURN
*
      END
