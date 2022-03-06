      SUBROUTINE XINIT(I1)
*
*
*       Initialize *X arrays.
*       --------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'kscmn6.h'
*
      IPAIR = KVEC(I1)
      I2 = I1 + 1
      I = N + IPAIR
*
      LKSINT(IPAIR) = .TRUE.
      IFLAG(IPAIR)  = 0
*
      IPHASEX(IPAIR) = 0
      JCOMPX (IPAIR) = 0
      JCLOSEX(IPAIR) = 0
      JCMAXX (IPAIR) = 0
      KSPAIRX(IPAIR) = 0
      KS2X   (IPAIR) = 0
      EBCH0X (IPAIR) = 0.d0
*
      T0X(I1) = T0(I1)
*
      DO 10 K = 1,3
          XKS(K,I1)    = X(K,I1)
          XKS(K,I2)    = X(K,I2)
          XDOTKS(K,I1) = XDOT(K,I1)
          XDOTKS(K,I2) = XDOT(K,I2)
   10 CONTINUE
*
      NB = LIST(1,I1)
      DO 20 K = 1,NB+1
          LISTX(K,I1) = LIST(K,I1)
   20 CONTINUE
*
      DO 30 K = 1,4
          U0X(K,IPAIR)     = U0(K,IPAIR)
          UDOTX(K,IPAIR)   = UDOT(K,IPAIR)
          FUX(K,IPAIR)     = FU(K,IPAIR)
          FUDOTX(K,IPAIR)  = FUDOT(K,IPAIR)
          FUDOT2X(K,IPAIR) = FUDOT2(K,IPAIR)
          FUDOT3X(K,IPAIR) = FUDOT3(K,IPAIR)
   30 CONTINUE
*
      HX(IPAIR)     = H(IPAIR)
      HDOTX(IPAIR)  = HDOT(IPAIR)
      HDOT2X(IPAIR) = HDOT2(IPAIR)
      HDOT3X(IPAIR) = HDOT3(IPAIR)
      HDOT4X(IPAIR) = HDOT4(IPAIR)
      DTAUX(IPAIR)  = DTAU(IPAIR)
      TDOT2X(IPAIR) = TDOT2(IPAIR)
      TDOT3X(IPAIR) = TDOT3(IPAIR)
      RX(IPAIR)     = R(IPAIR)
      GAMMAX(IPAIR) = GAMMA(IPAIR)
      KSLOWX(IPAIR) = KSLOW(IPAIR)
*
   40 RETURN
*
      END
