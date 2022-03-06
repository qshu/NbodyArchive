      SUBROUTINE REL_BIN_PARAMS (JPAIR)
*
*
*       Output of relativistic binaries.
*       --------------------------------
*
*       By Jonathan M.B. Downing 06.2007
*
      INCLUDE 'common6.h'
      INTEGER NAME1,NAME2
      REAL*8  ORBSEP,VOVERC1,VOVERC2
*
*       Set up the binary counters.
      J2 = JPAIR + JPAIR
      J1 = J2 - 1
*
*       Calculate basic orbital elements.
      ORBSEP = SQRT((X(1,J1)-X(1,J2))**2 + (X(2,J1)-X(2,J2))**2
     &              + (X(3,J1)-X(3,J2))**2)
*
      VR2 = (XDOT(1,J1) - XDOT(1,J2))**2 +
     &      (XDOT(2,J1) - XDOT(2,J2))**2 +
     &      (XDOT(3,J1) - XDOT(3,J2))**2
*
      EREL = 0.5*VR2 - (BODY(J1)+BODY(J2))/ORBSEP
      VOVERC1 = SQRT(XDOT(1,J1)*XDOT(1,J1)+XDOT(2,J1)*XDOT(2,J1)
     &     + XDOT(3,J1)*XDOT(3,J1))/CLIGHT
      VOVERC2 = SQRT(XDOT(1,J2)*XDOT(1,J2)+XDOT(2,J2)*XDOT(2,J2)
     &     + XDOT(3,J2)*XDOT(3,J2))/CLIGHT
*
*       Output for bound orbits.
      IF (EREL.LE.0.0) THEN
*
         SEMI = -0.5*(BODY(J1) + BODY(J2))/EREL
*
         ZN = SQRT((BODY(J1) + BODY(J2))/SEMI**3)
*
         RDOT = (X(1,J1) - X(1,J2))*(XDOT(1,J1) - XDOT(1,J2)) +
     &          (X(2,J1) - X(2,J2))*(XDOT(2,J1) - XDOT(2,J2)) +
     &          (X(3,J1) - X(3,J2))*(XDOT(3,J1) - XDOT(3,J2))
*
         ECC2 = (1.0 - ORBSEP/SEMI)**2 +
     &        RDOT**2/(SEMI*(BODY(J1) + BODY(J2)))
*
         ECC = SQRT(ECC2)
*
         RI = SQRT((X(1,J1) - RDENS(1))**2 + (X(2,J1) - RDENS(2))**2 + 
     &         (X(3,J1) - RDENS(3))**2)
*
*       Output binary parameters to orb.35
         WRITE (35,35)  TIME+TOFF,JPAIR,NAME(J1),NAME(J2),BODY(J1),
     &        BODY(J2),KSTAR(J1),KSTAR(J2),EREL,SEMI,ZN,
     &        ORBSEP,RI,ECC,GAMMA(JPAIR),GAMMAREL1(JPAIR),
     &        GAMMAREL2(JPAIR),GAMMAREL2_5(JPAIR),X(1,J1),
     &        X(2,J1),X(3,J1),X(1,J2),X(2,J2),X(3,J2),XDOT(1,J1),
     &        XDOT(2,J1),XDOT(3,J1),XDOT(1,J2),XDOT(2,J2),XDOT(3,J2),
     &        VOVERC1,VOVERC2,RELENERGY1(JPAIR),RELENERGY2(JPAIR),
     &        RELENERGY2_5(JPAIR)
*
      END IF
*
*       If we also want ouput on hyperbolic encounters.
      IF (KZ(46).EQ.2.AND.EREL.GT.0.0) THEN
*
         SEMI = ABS(-0.5*(BODY(J1) + BODY(J2))/EREL)
*
         ZN = SQRT((BODY(J1) + BODY(J2))/SEMI**3)
*
         RDOT = (X(1,J1) - X(1,J2))*(XDOT(1,J1) - XDOT(1,J2)) +
     &          (X(2,J1) - X(2,J2))*(XDOT(2,J1) - XDOT(2,J2)) +
     &          (X(3,J1) - X(3,J2))*(XDOT(3,J1) - XDOT(3,J2))
*
         ECC2 = (1.0 - ORBSEP/SEMI)**2 +
     &               RDOT**2/(SEMI*(BODY(J1) + BODY(J2)))
*
         ECC = SQRT(ECC2)
*
         RI = SQRT((X(1,J1) - RDENS(1))**2 + (X(2,J1) - 
     &        RDENS(2))**2 + (X(3,J1) - RDENS(3))**2)
*
*       Output binary parameters to orb.35
         WRITE (35,35)  TIME+TOFF,JPAIR,NAME(J1),NAME(J2),BODY(J1),
     &        BODY(J2),KSTAR(J1),KSTAR(J2),EREL,SEMI,ZN,
     &        ORBSEP,RI,ECC,GAMMA(JPAIR),GAMMAREL1(JPAIR),
     &        GAMMAREL2(JPAIR),GAMMAREL2_5(JPAIR),X(1,J1),
     &        X(2,J1),X(3,J1),X(1,J2),X(2,J2),X(3,J2),XDOT(1,J1),
     &        XDOT(2,J1),XDOT(3,J1),XDOT(1,J2),XDOT(2,J2),XDOT(3,J2),
     &        VOVERC1,VOVERC2,RELENERGY1(JPAIR),RELENERGY2(JPAIR),
     &        RELENERGY2_5(JPAIR)
*
      END IF
*
*       Format for binary output on unit 35
   35    FORMAT (F12.6,3I6,2E14.6,2I10,27E14.6)
*
      RETURN
*
      END
