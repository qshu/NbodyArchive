      SUBROUTINE STAR( MASS, TM, TN, TSCLS, LUMS)
*
*
*       Stellar luminosity & evolution time.
*       ------------------------------------
*
      REAL  MASS, TSCLS(6), LUMS(6), TM, TN, AM, CL1, CL2, TGS, LX
*
*       Computes the characteristic luminosities at different stages (LUMS),
*       and the time spent in each stage (TSCLS).
*       Ref: P.P. Eggleton, M.J. Fitchett & C.A. Tout (1989) Ap.J. 347, 998.
*
*       -----------------------------------------------------------
*       Times: 1; Hertz gap          2; Base of GB    3; He burning
*              4; A.G.B.             5; TGS           6; TB
*
*       LUMS:  1; ZAMS               2; MS            3; LG
*              4; He ign, AGB begins 5; He burning    6; Top of AGB
*       -----------------------------------------------------------
*
*
      AM = ALOG10( MASS )
      IF ( MASS .LT.1.0983) THEN
          LX = (1.107* MASS **3 + 240.7* MASS **9)/
     &                                           (281.9* MASS **4 + 1.0)
      ELSE
          LX = 1.399E4* MASS **5 /
     &        ( MASS **4 + 2.151E3* MASS **2 + 3.908E3 * MASS + 9.536E3)
      END IF
*
      IF ( MASS .GT.1.334) THEN
          CL1 = .092088 + .059338*AM
          CL2 = -.05713 + AM*(.3756 - .1744*AM)
      ELSE
          CL1 = .2594 + .1348*AM
          CL2 = .144 - .8333*AM
      END IF
*
      LUMS(1) = LX
      LUMS(2) = LX *10 **(CL1 + Cl2)
      LUMS(3) = (2.15 + 0.22* MASS **3)* MASS * MASS /
     &                    (5.0E-6* MASS **4 + 1.4E-2* MASS * MASS + 1.0)
      LUMS(4) = LUMS(3) + 2.0E3
      LUMS(5) = 0.763* MASS **0.46*LX + 50.0 / MASS **0.1
      LUMS(6) =  MASS *(4.0E3 + 5.0E2* MASS )
*
      TM = (2.5E3 + MASS **2.5*( 6.7E2 + MASS **2 )) /
     &                              (.033* MASS **1.5 + .35* MASS **4.5)
      TGS = 0.15*TM
      TSCLS(1) = 0.54*TM / ( MASS *( MASS - 2.1) + 23.3)
      TSCLS(2) = TGS - TGS * ( LUMS(3) / LUMS(4) )**0.855
      TSCLS(3) = TM*LX / (LUMS(5)*( MASS **0.42 + 0.8))
      TSCLS(5) = TGS
      TSCLS(6) = TGS*(LUMS(3) / LUMS(6))**0.855
      TSCLS(4) = TGS - TSCLS(2) - TSCLS(6)
      TN = TM + TSCLS(1) + TSCLS(2) + TSCLS(3) + TSCLS(4)
*
      RETURN
*
      END
