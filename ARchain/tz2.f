      
      CLIGHT = 18000.0
      ECC = 0.308
      SEMI = 7.0D-07
      ECC = 0.99
      SEMI = 2.0D-05
      RAU = 2.0*2.0D+05
      SMU = 6.1D+04
      CLIGHT = 25000.0
      ECC = 0.999
      ECC = 0.9999
      ECC = 0.999
      SEMI = 3.0D-08
      SMU = 4.4D+07
      SEMI = 1.0D-04
      SEMI = 1.0D-05
      SEMI = 1.0D-07
      ECC = 0.0
      CLIGHT = 20000.0
      SEMI = 1.98D-05
      ECC = 0.979
      
      SMU = 8.9D+07
      TAUGR = 1.3D+18*RAU**4/SMU**3
      WRITE (6,1) TAUGR
    1 FORMAT (' TAUGR   ',1P,E10.2)
          ECC2 = ECC**2
          FE = 1.0 + (73.0/24.0 + 37.0*ECC2/96.0)*ECC2
          GE = (1.0 - ECC2)**3.5/FE
          ZX = 3.0D-04
      ZX = 8.0D-04
          RATIO = 1.0
      ZX = 3.0D-04
      RATIO = 0.05
      ZX = 1.4D-03
      RATIO = 1.0
*         RATIO = 7.6/0.8
*       Replace physical time-scale by N-body units (cf. Lee 1993).
*         TZ = TAUGR*GE*SEMI**4/(RATIO*(1.0 + RATIO)*ZX**3)
          TZ = GE*SEMI**4/(RATIO*(1.0 + RATIO)*ZX**3)
      WRITE (6,3)  SEMI, ZX, TZ
    3 FORMAT (' SEMI ZX TZ  ',1P,3E10.2)
          TZ = 5.0/64.0*CLIGHT**5*TZ
      WRITE (6,6)  TZ
    6 FORMAT (' TZ  ',1P,E10.2)
      DW = 6.0*3.14*ZX/(SEMI*(1.0 - ECC2)*CLIGHT**2)
      DW2 = 5.0*DW
      WRITE (6,10)  ECC, SEMI, ZX, DW, DW2
   10 FORMAT (' EINSTEIN    E A MX DW DW2 ',F10.5,1P,4E10.2)
      STOP
      END 
