      IMPLICIT REAL*8 (A-H,M,O-Z)      
*
*     Sambaran's formula.

      M1 = 20.0
      M2 = 20.0
      M = M1 + M2
      ECC = 0.999
      ECC2 = ECC**2
      FE = (1.0 - ECC2)**3.5
      SEMI = 215.0
      SEMI = 4000.0
      TZ = (1.0/M)**3*SEMI**4*FE
      TZ = 150.0*TZ
      WRITE (6,5) ECC, TZ
    5 FORMAT (' ECC TZ  ',F9.5,1P,E10.2)
      
      END 
