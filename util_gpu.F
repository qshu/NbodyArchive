C--------------------------------------------------------------------
C     Calculate the REG force on local GPU memory
C--------------------------------------------------------------------
      SUBROUTINE CALC_REG_ON_GPU(IREG,istart,iend)

      INCLUDE 'common6.h'

      PARAMETER (maxthr=1024)
      INTEGER*4 IREG(NMAX), istart, iend
      REAL*8 a2by18(ID), a1by6(ID), aby2(ID)
      REAL*8 eps2, h2_i(NPIPE)
      REAL*8 x_i(ID,NPIPE), v_i(ID,NPIPE)
      REAL*8 p_i(NPIPE), a_i(ID,NPIPE), jerk_i(ID,NPIPE)
      INTEGER LISTGP(NMAX,maxthr)

      NOFL = 0
C---  calculate the pot, a & jerk on the GPU

      DO L=istart,iend,maxthr

c      WRITE(6,*) 'LOOP REG', RANK, L, LIST(1,1)

         ni = maxthr
         IF( (L+ni).GT.iend ) ni = iend - L + 1

C---  fill the GPU "working" arrays

         DO ii=1,ni

         idi       = IREG(L+ii-1)
ccc         ind_i(ii) = NAME(idi)
         ind_i(ii) = idi

CCC      Obtain irregular & regular force and determine current neighbours.
CCC      Take volume between inner and outer radius equal to basic sphere.

         RS2      = RS(idi)**2
	 RCRIT2   = 1.59*RS2
         h2_i(ii) = RCRIT2

C      WRITE(6,*) RANK, 'GPU UTIL h2', ii-1, idi, sqrt(h2_i(ii)), 
C     &NAME(idi), LIST(1,idi)

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



CCC      WRITE(6,*) ii-1, idi, h2_i(48), NAME(idi)


C---  GPU first call with first values of p_i, a_i & jerk_i

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
     *            WRITE (6,56)  NSTEPR, NAME(I), LIST(1,I), RS(I)
   56             FORMAT (' OVERFLOW!   #R NAME NB0 RS ',I11,I6,I4,F8.2)
                  CALL FLUSH(6)
*      Consider calling regint for non-GPU calculation or reduce RCRIT2.
                  NOFL = NOFL + 1
                  GO TO 54
              END IF

C---  Copy neighbour list but skip self-interaction.
              L1 = 1
              DO 53 L = 2,NNB+1
*       Note GPU address starts from 0 (hence add IFIRST to neighbour list).
                  IF (LISTGP(L,ii)+IFIRST.NE.I) THEN
                      L1 = L1 + 1
                      LIST_G(L1,idi) = LISTGP(L,ii) + IFIRST
                  END IF
   53         CONTINUE
              LIST_G(1,idi) = L1 - 1
              PHIDBL(idi) = p_i(ii)

        DO k=1,3
           FRG(k,idi)    = a_i(k,ii)
           FRGDOT(k,idi) = jerk_i(k,ii)
        ENDDO				! k


       WRITE(6,*) 'FORCE REG', RANK, ii, p_i(ii), a_i(1,ii), jerk_i(1,ii)


       RETURN
C--------------------------------------------------------------------
C     Calculate the REG force on local GPU memory
C--------------------------------------------------------------------
 
