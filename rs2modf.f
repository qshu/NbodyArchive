      SUBROUTINE RS2MODF(fct,massa,massb)
*
*       Subroutine to change neighbour radius according to mass differences
*
      REAL*8 massa, massb, fct
      REAL*8 x
*
      REAL*8 abeta, lbda
      INTEGER itrym
      COMMON/RS2PAR/ abeta,lbda,itrym
*
*
      IF((massa.GT.0.0).AND.(massb.GT.0.0)) THEN
         x = massa / massb
*
      itrym=3
      abeta=0.03
      lbda=1
*
*       Comment either function (Hotshot)
         IF(itrym.GT.0) THEN
            IF(itrym.EQ.1) fct = EXP( abeta*(((x + 1/x)*0.5)**lbda - 1))
            IF(itrym.EQ.2) fct = abeta*((x + 1/x - 1)**lbda - 1) + 1
            IF(itrym.EQ.3) fct = abeta*(((x + 1/x)*0.5)**lbda - 1) + 1
            IF(itrym.EQ.4) THEN
               IF(massa.NE.massb) THEN
                  fct = abeta
               ELSE
                  fct = 1.0
               END IF
            END IF
         ELSE
            fct = 1.0
         END IF
      ELSE
         fct = 1.0
      END IF
*
      END
