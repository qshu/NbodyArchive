      SUBROUTINE XUPDATE(I1)
*
*
*       Update common variables from *X ones.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'kscmn6.h'
*
      IPAIR = KVEC(I1)
*     I2 = I1 + 1
      I = N + IPAIR
*
      LKSINT(IPAIR) = .FALSE.
*
      IF (IFLAG(IPAIR).EQ.1) THEN
*       From: BRAKE4
*       Include optional kick velocity of 5*VRMS km/s for coalescence recoil.
          VI20 = 0.0
          DO 40 K = 1,3
              X0(K,I) = X(K,I)
              X0DOT(K,I) = XDOT(K,I)
              VI20 = VI20 + XDOT(K,I)**2
   40     CONTINUE
*
          VF = 5.0*(VRMS/VSTAR)/SQRT(VI20)
          DO 50 K = 1,3
              XDOT(K,I) = VF*XDOT(K,I)
              X0DOT(K,I) = XDOT(K,I)
   50     CONTINUE
          ECD0 = ECDOT
          ECDOT = ECDOT + 0.5*BODY(I)*VI20*(1.0 - VF**2)
          VESC = 5.0*VRMS
          WRITE (6,60)  VF, ECD0-ECDOT, VESC
   60     FORMAT (' COALESCENCE KICK    VF ECDOT VESC ',
     &                                  F7.3,F10.6,F6.1)
*       Form neighbour list and new polynomials.
          RS0 = RS(I)
          CALL NBLIST(I,RS0)
          CALL FPOLY1(I,I,0)
          CALL FPOLY2(I,I,0)
      END IF
*
      RETURN
*
      END
