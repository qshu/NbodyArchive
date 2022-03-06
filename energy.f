      SUBROUTINE ENERGY
*
*
*       Total energy.
*       -------------
*
      INCLUDE 'common6.h'

      REAL*8 CMVI(3),RIJ,MU,ETAM,V2,RIJDOT,E1PN,E2PN,E2_5PN,M12
      REAL*8 RIJ2,CL2,CL4,CL5,NX,NY,NZ,NV1,NV2,V12,V22,V1V2,RELCOR
*
*       Sum the total energy of regularized pairs.
      EBIN = 0.0D0
      DO 10 IPAIR = 1,NPAIRS
*       Skip pairs with zero mass of c.m. particle (merged binary ghost).
          IF (BODY(N+IPAIR).GT.0.0D0) THEN
*       Predict coordinates, velocities & binding energy.
              CALL RESOLV(IPAIR,1)
              EBIN = EBIN + BODY(2*IPAIR-1)*BODY(2*IPAIR)*HT/
     &                                                     BODY(N+IPAIR)

          END IF
   10 CONTINUE
*

****************************************RELATIVISTIC ENERGY (Gabor Kupi)
      CL2 = CLIGHT*CLIGHT
      CL4 = CL2*CL2
      CL5 = CL4*CLIGHT

      RELCOR = 0.0
      IF (TIME.GE.0.0.AND.CLIGHT.GT.0.0) THEN


       DO 111 I = 1,N
       DO 222 J = I+1,N

      IF(BODY(I).GT.0.0.AND.BODY(J).GT.0.0) THEN

      M12 = BODY(I)+BODY(J)


      CMVI(1) =(BODY(I)*XDOT(1,I)+BODY(J)*XDOT(1,J))/M12
      CMVI(2) =(BODY(I)*XDOT(2,I)+BODY(J)*XDOT(2,J))/M12
      CMVI(3) =(BODY(I)*XDOT(3,I)+BODY(J)*XDOT(3,J))/M12


         XDOT(1,I) = XDOT(1,I) - CMVI(1)
         XDOT(2,I) = XDOT(2,I) - CMVI(2)
         XDOT(3,I) = XDOT(3,I) - CMVI(3)
         XDOT(1,J) = XDOT(1,J) - CMVI(1)
         XDOT(2,J) = XDOT(2,J) - CMVI(2)
         XDOT(3,J) = XDOT(3,J) - CMVI(3)



         RIJ2 = (X(1,I)-X(1,J))*(X(1,I)-X(1,J))+
     &(X(2,I)-X(2,J))*(X(2,I)-X(2,J))+(X(3,I)-X(3,J))*(X(3,I)-X(3,J))
         RIJ = SQRT(RIJ2)
         MU = BODY(I)*BODY(J)/M12
         ETAM = MU/M12
         V2 = (XDOT(1,I)-XDOT(1,J))*(XDOT(1,I)-XDOT(1,J))+
     &(XDOT(2,I)-XDOT(2,J))*(XDOT(2,I)-XDOT(2,J))+
     &(XDOT(3,I)-XDOT(3,J))*(XDOT(3,I)-XDOT(3,J))
         RIJDOT = ((X(1,I)-X(1,J))*(XDOT(1,I)-XDOT(1,J))+
     &(X(2,I)-X(2,J))*(XDOT(2,I)-XDOT(2,J))+(X(3,I)-X(3,J))
     &*(XDOT(3,I)-XDOT(3,J)))/RIJ

         NX = (X(1,I)-X(1,J))/RIJ
         NY = (X(2,I)-X(2,J))/RIJ
         NZ = (X(3,I)-X(3,J))/RIJ

         NV1 = NX*XDOT(1,I)+NY*XDOT(2,I)+NZ*XDOT(3,I)
         NV2 = NX*XDOT(1,J)+NY*XDOT(2,J)+NZ*XDOT(3,J)

*** L. Blanchet, Phys. Rev. D 54, 1417-1438 (1996), gr-qc/9603048
*** A. Gopakumar, B. I. Iyer and S. Iyer, Phys. Rev. D 55, 6030-6053 (1997), gr-qc/9703075
*** A. Gopakumar, B. I. Iyer and S. Iyer, Phys. Rev. D 57, 6562 (1998)

         E1PN = MU*(0.375*(1.0-3.0*ETAM)*V2*V2+0.5*(3.0+ETAM)
     &*V2*M12/RIJ+0.5*ETAM*M12*RIJDOT*RIJDOT/RIJ+0.5*M12*M12
     &/RIJ2)/CL2

       V1V2 = XDOT(1,I)*XDOT(1,J)+XDOT(2,I)*XDOT(2,J)
&       +XDOT(3,I)*XDOT(3,J)



         E2PN = MU*(0.3125*(1.0-7.0*ETAM+13.0*ETAM*ETAM)*V2*V2*V2
     &+0.125*(21.0-23.0*ETAM-27.0*ETAM*ETAM)*V2*V2*M12/RIJ
     &+0.25*ETAM*(1.0-15.0*ETAM)*M12*V2*RIJDOT*RIJDOT/RIJ
     &-0.375*ETAM*(1.0-3.0*ETAM)*M12*RIJDOT*RIJDOT*RIJDOT*RIJDOT/RIJ
     &-0.25*(2.0+15.0*ETAM)*M12*M12*M12/RIJ2/RIJ
     &+0.125*(14.0-55.0*ETAM+4.0*ETAM*ETAM)*V2*M12*M12/RIJ2
     &+0.125*(4.0+69.0*ETAM+12.0*ETAM*ETAM)*RIJDOT*RIJDOT*M12*M12
     &/RIJ2)/CL4

         E2_5PN = 1.6*BODY(I)*BODY(I)*BODY(J)*BODY(J)*(NV1-NV2)*V2
     &/M12/RIJ2/CL5

         XDOT(1,I) = XDOT(1,I) + CMVI(1)
         XDOT(2,I) = XDOT(2,I) + CMVI(2)
         XDOT(3,I) = XDOT(3,I) + CMVI(3)
         XDOT(1,J) = XDOT(1,J) + CMVI(1)
         XDOT(2,J) = XDOT(2,J) + CMVI(2)
         XDOT(3,J) = XDOT(3,J) + CMVI(3)

         RELCOR = RELCOR + E1PN + E2PN + E2_5PN




       END IF

 222   CONTINUE
 111   CONTINUE
       END IF

      SRELENERGY = 0.0
      DO 500 KK = 1,KMAX
       SRELENERGY = SRELENERGY + RELENERGY(KK)
!!!!!!!!!!!!!!!!!!
!	WRITE(*,999)'ENERG1',TIME,KK,RELENERGY(KK),SRELENERGY
!!!!!!!!!!!!!!!!!!
  500 CONTINUE
  999 FORMAT(A6,F45.25,I6,2F45.25)

	EBIN = EBIN!  - SRELENERGY
      WRITE(*,*) 'TOTAL REL ENERGY',TIME,RELCOR



*       Calculate the potential energy.
      ZKIN = 0.D00
      POT = 0.0
*
      DO 20 I = 1,NTOT
      JMIN = I + 1
      IF (I.LE.2*NPAIRS) THEN
*       Binding energy of regularized pairs is included explicitly above.
          IPAIR = KVEC(I)
          JMIN = 2*IPAIR + 1
      END IF
*
      IPAIR = 0
      IF (I.GT.N)  THEN
*       Binding energy at center of mass position without binary members
          IPAIR = I - N
      END IF
*
      POTJ = 0.D00
      POTI = 0.D00
*       POTI contains potential at particles position to be stored later (R.Sp.)
*
      DO 30 J = 1,N
      IF (J.EQ.I .OR. J.EQ.2*IPAIR-1 .OR. J.EQ.2*IPAIR .OR.
     *    BODY(J).EQ.0.0D0 .OR. BODY(I).EQ.0.0D0)  GO TO 30
          A1 = X(1,I) - X(1,J)
          A2 = X(2,I) - X(2,J)
          A3 = X(3,I) - X(3,J)
      A4 = BODY(J)/DSQRT (A1*A1 + A2*A2 + A3*A3)
      POTI = POTI - A4
*  also J.LT.N?
      IF(J.GE.JMIN)POTJ = POTJ + A4
   30 CONTINUE
*       Store potential in shared vector first (R.Sp.)
      PHIDBL(I) = POTI
      POT = POT + BODY(I)*POTJ
   20 CONTINUE
*
*       Sum the kinetic energy (include c.m. bodies but not components).
      DO 40 I = IFIRST,NTOT
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
   40 CONTINUE
*
      ZKIN = 0.5D0*ZKIN
*
*       Obtain the tidal potential if external field is present.
      ETIDE = 0.0D0
      IF (KZ(14).GT.0) THEN
          CALL XTRNLV(1,N)
      END IF
*
*       Check differential potential energy due to chain subsystem.
      IF (NCH.GT.0) THEN
          CALL CHPOT(DP)
          POT = POT + DP
      END IF
*
*       Total energy = ZKIN - POT + ETIDE + EBIN + ESUB + EMERGE + ECOLL.
*
      RETURN
*
      END
