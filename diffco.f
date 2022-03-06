      SUBROUTINE DIFFCO(NXTLEN,NXTLST)
*
*     Update Energy Differences
*
      INCLUDE 'common6.h'
      LOGICAL LPOSE(NMAX),LPOSJ(NMAX),LPR
      REAL*4 TRANSE,TRANSJ
      COMMON/DIFF/DCE(NMAX),DCE1(NMAX),DCE2(NMAX),DCE3(NMAX),
     *            DCJ(NMAX),DCJ1(NMAX),DCJ2(NMAX),DCJ3(NMAX),
     *        DCEP(NMAX),DETL(NMAX),DCJP(NMAX),DJTL(NMAX),
     *        T1(NMAX),T2(NMAX),T3(NMAX),
C     ********  extra storage for previous E at maxima in J and vice versa
     $        DCEPP(NMAX),DCJPP(NMAX)
C     ********
      COMMON/POTENT/PHII(NMAX),PHIR(NMAX),PHIR1(NMAX)
      COMMON/STORE/TRANSE(9,1000),TRANSJ(9,1000),ITRE,ITRJ
*
      DO 1000 L = 1,NXTLEN
*
      I = NXTLST(L)
      LPR = I.EQ.96
*
      RDOT2 = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
      VDOT2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
      RDOTV = X(1,I)*XDOT(1,I) + X(2,I)*XDOT(2,I) + X(3,I)*XDOT(3,I)
      DTR = TIME - T0R(I)
      PHINEW = PHII(I) + PHIR(I) + DTR*PHIR1(I)
      ENEW = 0.5D0*VDOT2 + PHINEW
      RNEW = RDOT2*(2.D0*(ENEW-PHINEW)-RDOTV*RDOTV)
      RNEW = DSQRT(RNEW)
*
      IF(T1(I).GT.0.D0)THEN
*
      A1 = 1.D0/(TIME - T1(I))
      B = (ENEW - DCE(I))*A1
      BJ = (RNEW - DCJ(I))*A1
*
      IF(T3(I).GT.0.D0)THEN
*       Construct divided differences
      A2 = 1.D0/(TIME - T2(I))
      A3 = 1.D0/(TIME - T3(I))
      C = (B - DCE1(I))*A2
      D = (C - DCE2(I))*A3
      CJ = (BJ - DCJ1(I))*A2
      DJ = (CJ - DCJ2(I))*A3
*
      DCE(I) = ENEW
      DCE1(I) = B
      DCE2(I) = C
      DCE3(I) = D
      DCJ(I) = RNEW
      DCJ1(I) = BJ
      DCJ2(I) = CJ
      DCJ3(I) = DJ
*
      B = DCE1(I) + DCE2(I)*(TIME-T1(I)) +
     *              DCE3(I)*(TIME-T1(I))*(TIME-T2(I))
      BJ = DCJ1(I) + DCJ2(I)*(TIME-T1(I)) +
     *              DCJ3(I)*(TIME-T1(I))*(TIME-T2(I))
*
      IF(LPR)PRINT*,' TIME,ENEW,B=',TIME,ENEW,B
      IF(LPR)PRINT*,' TIME,RNEW,BJ=',TIME,RNEW,BJ
*
      IF(LPOSE(I).AND.B.LT.0.D0)THEN
*       Check for transition with second derivative
*
      C = 2.D0*(DCE2(I) + DCE3(I)*(2.D0*TIME - T1(I) - T2(I)))
*       Search for maximum
      IF(C.LT.0.D0)THEN
*
      PRINT*,' Maximum Energy found C=',C
      IF(DETL(I).GT.0.D0)THEN
*       Output after second maximum found
*
      ITRE = ITRE + 1
      DELTAT = TIME - DETL(I)
      DELTAE = ENEW - DCEP(I)
c     *************** change in J
      DeltaJE = RNEW - DCJPP(I)
C     *************** 
*
          TRANSE(1,ITRE) = REAL(I)
          TRANSE(2,ITRE) = TIME
          TRANSE(3,ITRE) = RDOT2
          TRANSE(4,ITRE) = DELTAT
          TRANSE(5,ITRE) = DELTAE
          TRANSE(6,ITRE) = ENEW
          TRANSE(7,ITRE) = PHINEW
C     you need to store previous J as well
          TRANSE(8,ITRE)  = DCJPP(I)

C     *************** change in J
          TRANSE(9,ITRE)  = DeltaJE
C     ***************
*
          IF(ITRE.EQ.100)CALL WTRANS(1)
*
      END IF
*       Update values of last transition
      DETL(I) = TIME
      DCEP(I) = ENEW
C     ************** store current J
      DCJPP(I) = Rnew
c     **************
      END IF
*
      END IF
*       Now check for angular momentum transition
      IF(LPOSJ(I).AND.BJ.LT.0.D0)THEN
*       Check for transition with second derivative
*
      CJ = 2.D0*(DCJ2(I) + DCJ3(I)*(2.D0*TIME - T1(I) - T2(I)))
*       Search for maximum
      IF(CJ.LT.0.D0)THEN
*
      PRINT*,' Maximum angular momentum found CJ=',CJ
      IF(DJTL(I).GT.0.D0)THEN
*       Output after second maximum found
*
      ITRJ = ITRJ + 1
      DELTAT = TIME - DJTL(I)
      DELTAJ = RNEW - DCJP(I)
c     ************* change in E
      DeltaEJ = Enew - DCEPP(I)
c     *************
*
          TRANSJ(1,ITRJ) = REAL(I)
          TRANSJ(2,ITRJ) = TIME
          TRANSJ(3,ITRJ) = RDOT2
          TRANSJ(4,ITRJ) = DELTAT
          TRANSJ(5,ITRJ) = DELTAJ
          TRANSJ(6,ITRJ) = RNEW
          TRANSJ(7,ITRJ) = PHINEW
C     you need to store previous E as well
          TRANSJ(8,ITRJ) = DCEPP(I)
c     ************ store change in E
          TRANSJ(9,ITRJ) = DeltaEJ
c     ************
*
          IF(ITRJ.EQ.100)CALL WTRANS(2)
*
      END IF
*       Update values of last transition
      DJTL(I) = TIME
      DCJP(I) = RNEW
C     ************ update E
      DCEPP(I) = ENEW
C     ************
      END IF
*
      END IF
*       End of check for transition
      END IF
*
      LPOSE(I) = B.GT.0.D0
      LPOSJ(I) = BJ.GT.0.D0
*
      END IF    
*
      T3(I) = T2(I)
      T2(I) = T1(I)
      T1(I) = TIME
*
 1000 CONTINUE
*
      RETURN
      END
*
      SUBROUTINE WTRANS(I)
*
      REAL TRANSE,TRANSJ
      COMMON/STORE/TRANSE(9,1000),TRANSJ(9,1000),ITRE,ITRJ
*        Write energy diffusion coefficients
      IF(I.EQ.1)THEN
*
      IF(IRE.EQ.0)THEN
*
 777  CONTINUE
      READ(12,200,END=888,ERR=999)((TRANSE(K,J),K=1,7),J=1,1000)
      GOTO 777
 888  IRE = 1
*
      END IF
      WRITE(12,200)((TRANSE(K,J),K=1,7),J=1,1000)
 200  FORMAT(1X,1P,7(1X,E12.5))
*
      ITRE = 0
      END IF
*        Write angular momentum diffusion coefficients
*        Write energy diffusion coefficients
      IF(I.EQ.2)THEN
*
      IF(IRJ.EQ.0)THEN
*
 778  CONTINUE
      READ(13,200,END=889,ERR=997)((TRANSJ(K,J),K=1,7),J=1,1000)
      GOTO 778
 889  IRJ = 1
*
      END IF
      WRITE(13,200)((TRANSJ(K,J),K=1,7),J=1,1000)
*
      ITRJ = 0
      END IF
*
      RETURN
 999  PRINT*,' Error in Read from Unit 12, routine WTRANS'
      STOP
 997  PRINT*,' Error in Read from Unit 13, routine WTRANS'
      STOP
      END




