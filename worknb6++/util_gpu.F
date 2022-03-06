C--------------------------------------------------------------------
C     Calculate the REG force on local GPU memory
C--------------------------------------------------------------------
      SUBROUTINE CALC_REG_ON_GPU(IREG,istart,iend)

      INCLUDE 'common6.h'

      PARAMETER (maxthr=1024)
      INTEGER*4 IREG(NMAX), istart, iend
      REAL*8 eps2, h2_i(maxthr)
      REAL*8 x_i(3,maxthr), v_i(3,maxthr)
      REAL*8 p_i(maxthr), a_i(3,maxthr), jerk_i(3,maxthr)
      INTEGER LISTGP(LMAX,maxthr)
      REAL*8 FIRR(3),FREG(3),FDR(3)
      LOGICAL LNB
      INTEGER NBFLAG(NMAX)

      NOFL = 0
C---  calculate the pot, a & jerk on the GPU

*      WRITE(6,*) 'LOOPREG rank t istart,iend ', RANK, time,istart,iend
*      WRITE(6,*) ' IREG=',(IREG(K),k=istart,iend)
*      CALL FLUSH(6)

      DO L=istart,iend,maxthr

         ni = maxthr
         IF( (L+ni).GT.iend ) ni = iend - L + 1

C---  fill the GPU "working" arrays

         DO ii=1,ni

         idi       = IREG(L+ii-1)
ccc         ind_i(ii) = NAME(idi)

CCC      Obtain irregular & regular force and determine current neighbours.
CCC      Take volume between inner and outer radius equal to basic sphere.

         RS2      = RS(idi)**2
	 RCRIT2   = 1.59*RS2
         h2_i(ii) = RCRIT2

*      WRITE(6,*) RANK, 'GPU UTIL h2', ii-1, idi, sqrt(h2_i(ii)), 
*    &NAME(idi), LIST(1,idi)

         DO k=1,3
            x_i(k,ii) = X(k,idi)
            v_i(k,ii) = XDOT(k,idi)
         ENDDO				! k

*        p_i(ii) = PHIDBL(idi)
*
*        DO k=1,3
*           a_i(k,ii)    = 2.0D0*F(k,idi)
*           jerk_i(k,ii) = 6.0D0*FDOT(k,idi)
*        ENDDO				! k
*
         ENDDO					! ii


*      WRITE(6,*) ' call gpu L ni lmax nbmax = ',L,ni, lmax,nnbmax
*      WRITE(6,*) ' h2,x,v=',h2_i(1),x_i(1,1),v_i(1,1)
*      CALL FLUSH(6)

C---  GPU first call with first values of p_i, a_i & jerk_i

     `
          CALL gpunb_regf(ni,h2_i,x_i,v_i,a_i,jerk_i,p_i,lmax,
     &                    nnbmax,LISTGP)

        gpu_calls = gpu_calls + 1
C---  Save new neighbour list on intermediate vector for regint
          DO 54 ii = 1,NI
             idi = IREG(L+ii-1)
             NNB = LISTGP(1,ii)
*       Warning if neighbour list overflows.
              IF (NNB.LT.0) THEN
                  if(rank.eq.0)
     *            WRITE (6,56)  NSTEPR, NAME(idi), LIST(1,idi), RS(idi)
   56             FORMAT (' OVERFLOW!   #R NAME NB0 RS ',I11,I6,I4,F8.2)
                  CALL FLUSH(6)
*      Define all particles as neighbours for overflow for full force loop in regint.
                  LIST_G(1,idi) = NNB
                  NOFL = NOFL + 1
                  GO TO 54
              END IF

C---  Copy neighbour list but skip self-interaction.
              DO 553 J = IFIRST,N
 553          NBFLAG(J) = 0
*
              L1 = 1
              DO 53 LL = 2,NNB+1
*       Note GPU address starts from 0 (hence add IFIRST to neighbour list).
                  IF (LISTGP(LL,ii)+IFIRST.NE.idi) THEN
                      L1 = L1 + 1
                      LIST_G(L1,idi) = LISTGP(LL,ii) + IFIRST
                      NBFLAG(LIST_G(L1,idi)) = 1
                  END IF
   53         CONTINUE
              LIST_G(1,idi) = L1 - 1
              PHIDBL(idi) = p_i(ii)

*       Perform fast force loop over single particles.
          DO K = 1,3
          FIRR(K) = 0.D0
          FREG(K) = 0.D0
          FDR(K) = 0.D0
          END DO
          PHIREG = 0.D0
*
          NNBCHK = 0
          RIJMAX = 0.D0
      DO 10 J = IFIRST,N
*RSP
          IF (J.EQ.idi) GO TO 10
*RSP
*
          A1 = X(1,J) - x_i(1,ii)
          A2 = X(2,J) - x_i(2,ii)
          A3 = X(3,J) - x_i(3,ii)
*       Predicted coordinates avoids spurious force differences.
          DV1 = XDOT(1,J) - v_i(1,ii)
          DV2 = XDOT(2,J) - v_i(2,ii)
          DV3 = XDOT(3,J) - v_i(3,ii)
*
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = A1*DV1 + A2*DV2 + A3*DV3
          DRDP = 3.0*DRDV*DR2I
*
*       See whether the distance exceeds the outer shell radius.
          IF (NBFLAG(J).GT.0) THEN
          IF (RIJ2.GT.h2_i(ii))PRINT*,
     * ' t,i,j,util_gpu warning rij2,rcrit2=',TIME,idi,J,RIJ2,h2_i(ii),
     *      ' dr=',RIJ2-h2_i(ii)
          NNBCHK = NNBCHK + 1
          FIRR(1) = FIRR(1) + A1*DR3I
          FIRR(2) = FIRR(2) + A2*DR3I
          FIRR(3) = FIRR(3) + A3*DR3I
          IF(DSQRT(RIJ2).GT.RIJMAX) THEN
             IJMAX = J
             RIJMAX = DSQRT(RIJ2)
             FIRRMAX = DSQRT(FIRR(1)**2+FIRR(2)**2+FIRR(3)**2)
             HMAX = DSQRT(h2_i(ii))
          END IF

          ELSE
*
*       Obtain the regular force.
          FREG(1) = FREG(1) + A1*DR3I
          FREG(2) = FREG(2) + A2*DR3I
          FREG(3) = FREG(3) + A3*DR3I
          FDR(1) = FDR(1) + (DV1 - A1*DRDP)*DR3I
          FDR(2) = FDR(2) + (DV2 - A2*DRDP)*DR3I
          FDR(3) = FDR(3) + (DV3 - A3*DRDP)*DR3I
*       Obtain potential and derivative.
          PHIREG = PHIREG - DR3I*RIJ2
          END IF
   10 CONTINUE
           DFREG = 0.D0
           DFDREG = 0.D0
        DO k=1,3
           FRG(k,idi)    = a_i(k,ii)
           FRGDOT(k,idi) = jerk_i(k,ii)
           DFREG = DFREG + (FRG(K,idi)-FREG(K))**2
           DFDREG = DFDREG + (FRGDOT(K,idi)-FDR(K))**2
        ENDDO                           ! k
           DFREG = DSQRT(DFREG)
           DFDREG = DSQRT(DFDREG)
        IF(DFREG.GT.1.D-4.OR.DFDREG.GT.1.D-4)THEN
        WRITE(6,*)' F Err t i1st n l ii idi nnb1/2 df dfd ',TIME,ifirst,
     * n,l,ii,idi,l1-1,nnbchk,dfreg,dfdreg
        PRINT*,' max i t r h f dr ',IJMAX,TIME,RIJMAX,HMAX,
     * FIRRMAX,HMAX-RIJMAX
       CALL FLUSH(6)

        END IF
   54   CONTINUE

*      WRITE(6,*) 'FORCE REG',RANK,p_i(1),a_i(1,1),jerk_i(1,1)
*      CALL FLUSH(6)

       ENDDO    !        DO L=istart,iend,maxthr



       RETURN
C--------------------------------------------------------------------
C     Calculate the REG force on local GPU memory
C--------------------------------------------------------------------
       END
