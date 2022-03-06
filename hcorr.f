      SUBROUTINE HCORR(I,DM,RNEW)
*
*
*       Mass loss correction of KS orbit.
*       ---------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Copy pair index and determine secondary component (note I1 = I here).
      IPAIR = KSPAIR
      I1 = I
      I2 = 2*IPAIR - 1
      IF (I2.EQ.I1) THEN
          I2 = I2 + 1
      END IF
*
*       Define c.m. index and save old elements & reduced mass.
      J = N + IPAIR
      SEMI0 = -0.5D0*BODY(J)/H(IPAIR)
*     ECC2 = (1.0 - R(IPAIR)/SEMI0)**2 + TDOT2(IPAIR)**2/(BODY(J)*SEMI0)
*     ECC = SQRT(ECC2)
      ZMU0 = BODY(I1)*BODY(I2)/BODY(J)
*
*       Obtain energy change from MK*A = const and H = -M/(2*A) (DCH 8/96).
      DH = DM/SEMI0*(1.0 - 0.5*DM/BODY(J))
*
*       Reduce mass of c.m. body and subtract energy correction terms.
      BODY(J) = BODY(J) - DM
      ZMU1 = (BODY(I1) - DM)*BODY(I2)/BODY(J)
      EMDOT = EMDOT - ZMU1*DH - (ZMU1 - ZMU0)*H(IPAIR)
*
*       Include diagnostic warning for large relative energy change.
      IF (DH.GT.0.2*ABS(H(IPAIR))) THEN
          WRITE (6,10)  NAME(I), DH, H(IPAIR), R(IPAIR)/SEMI0, DM*SMU
   10     FORMAT (' WARNING!    LARGE CORRECTION    NM DH H R/A DMS ',
     &                                              I6,2F8.2,2F7.2)
      END IF
*
*       Update the binding energy due to mass loss DM and set new radius.
      H(IPAIR) = H(IPAIR) + DH
      RADIUS(I) = RNEW
*
*       Set new mass temporarily to be consistent with routine EXPAND.
      BODY(I) = BODY(I) - DM
*
*       Modify KS variables at constant eccentricity.
      CALL EXPAND(IPAIR,SEMI0)
*
*       Restore the component mass here to allow updating in routine MDOT.
      BODY(I) = BODY(I) + DM
*
      RETURN
*
      END

