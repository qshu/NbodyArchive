      SUBROUTINE lagrNomass(C)
*
*
*       Lagrangian radii.
*       -----------------
*
      parameter (NLENS=18)
      INCLUDE 'common6.h'
      INCLUDE 'omp_lib.h'
      INTEGER NPARTC(NLENS), NCORE
      REAL*8  R2(NMAX),RSNGL(NMAX),RBIN(NMAX),C(3)
      REAL*8  FLAGR(NLENS),RLAGR(NLENS),AVMASS(NLENS),AVMRC
      REAL*8  RSLAGR(NLENS),RBLAGR(NLENS)
      REAL*8  SIGR2(NLENS),SIGT2(NLENS),SIG2(NLENS)
      REAL*8  SIGR2C,SIGT2C,SIG2C
      REAL*8  VROT(NLENS),VROTC
      REAL*8  VRAVE(NLENS),VTAVE(3,NLENS),VAVE(3,NLENS)
      REAL*8  VTAVEV(NLENS),VAVEV(NLENS)
      REAL*8  VRI(NMAX), VTI(3,NMAX),VTITMP(3),VITMP(3)
      REAL*8  VR_CORE,VT_CORE(3),V_CORE(3),V_COREV,VT_COREV
      INTEGER ISLIST(NMAX),IBLIST(NMAX)
*     --16/05/16 15:56-qishu-debug----------------------------*
***** Note:--------------------------------------------------**
*      RCn : core radius for MLP
      INTEGER NLIST(NMAX), Nnomass
      REAL*8  NLAGR(NLENS), Rnomass(NMAX)
      REAL*8  ZMASS0, RCn
      INTEGER I,J
*     --16/05/16 15:56-qishu-end------------------------------*
*
*     Lagrangian radii fraction of total mass
      DATA FLAGR/0.001D0,0.003D0,0.005D0,0.01D0,0.03D0,0.05D0,0.1D0,
     &     0.2D0,0.3D0,0.4D0,0.5D0,0.6D0,0.7D0,0.8D0,0.9D0,0.95D0,
     &     0.99D0,1.0D0/
*
*     half mass radius index
      DATA IRLAGRH /11/
          
      ZMASS0 = 0.0 
*     Get total mass
      DO I = IFIRST,N
         IF(nomass(I).eq.1) ZMASS0 = ZMASS0 + BODY(I)
      END DO
                               
*     Particle number counts
      NP = 0
      Nnomass = 0
      IF (KZ(8).GT.0) THEN
c         stop 'warnning: KZ(8) >= 0, not finished for lagrNomass.f'
*     Need to exclude massive black hole mass
*     Set square radii of resolved binaries
         IF (KZ(24).EQ.1) THEN
            NPP = 0
            NCMB = -100000
            DO I = 1,IFIRST-1
               IF(NAME(I).NE.NIMBH) THEN
                  NPP = NPP + 1
                  R2(NPP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                 (X(3,I) - C(3))**2
                  JLIST(NPP) = I
               ELSE
                  NCMB = KVEC(I) + N
               END IF
            END DO
*
*     Set square radii of single stalrs
            DO I = IFIRST,N
               IF(NAME(I).NE.NIMBH) THEN
                  NSNGL = NSNGL + 1
                  R2(NPP+NSNGL) = (X(1,I) - C(1))**2 + 
     &                 (X(2,I) - C(2))**2 + (X(3,I) - C(3))**2
                  RSNGL(NSNGL) = R2(NPP+NSNGL)
                  JLIST(NPP+NSNGL) = I
                  ISLIST(NSNGL) = I
               END IF
            END DO
*     Set square radii of c.m.
            DO I = N+1, NTOT
               IF(NAME(I).NE.NCMB) THEN
                  NBIN = NBIN + 1
                  RBIN(NBIN) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                 (X(3,I) - C(3))**2
                  IBLIST(NBIN) = I
               END IF
            END DO
            NP = NPP + NSNGL
         ELSE
*     Set square radii of resolved binaries
!$omp parallel do private(I)
            DO I = 1,IFIRST-1
               R2(I) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &              (X(3,I) - C(3))**2
               JLIST(I) = I
            END DO
!$omp end parallel do

*
*     Set square radii of single stalrs
!$omp parallel do private(I)
            DO I = IFIRST,N
               R2(I) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &              (X(3,I) - C(3))**2
               RSNGL(I-IFIRST+1) = R2(I)
               JLIST(I) = I
               ISLIST(I-IFIRST+1) = I
            END DO
!$omp end parallel do

*     Set square radii of c.m.
!$omp parallel do private(I)
            DO I = N+1, NTOT
               RBIN(I-N) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &              (X(3,I) - C(3))**2
               IBLIST(I-N) = I
            END DO
!$omp end parallel do
            NP = N
            NBIN = NTOT - N
            NSNGL = N - IFIRST + 1
         END IF

*     Sort square distances of all particles with respect to the centre C.
         CALL SORT1(NP,R2,JLIST)
*     Sort square distances of all singles with respect to the centre C.
         CALL SORT1(NSNGL,RSNGL,ISLIST)
*     Sort square distances of c.m. with respect to the centre C.
         CALL SORT1(NBIN,RBIN,IBLIST)
*
      END IF
c      ELSE
*     Only consider singles and resolved K.S.
*     Set square radii
*     exclude massive black hole
         IF (KZ(24).EQ.1) THEN
         stop 'warnning: KZ(24) = 1, not finished for lagrNomass.f'
         ELSE
!$omp parallel do private(I)
            DO I = 1,N
               R2(I) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &              (X(3,I) - C(3))**2
               JLIST(I) = I
            END DO
!$omp end parallel do
            NP = N
*     --16/05/16 15:59-qishu-debug----------------------------*
***** Note:--------------------------------------------------**
*     Set square radii of massless particles
            DO I = IFIRST,N
*               IF(nomass(I).eq.1) THEN
               IF(body(I).le.smallMass) THEN
           if (body(i) .gt. smallMass) stop "lagrN: massive nomass"
                  Nnomass = Nnomass + 1
                  Rnomass(Nnomass) = (X(1,I) - C(1))**2 + 
     &                 (X(2,I) - C(2))**2 + (X(3,I) - C(3))**2
                  NLIST(Nnomass) = I
               END IF
            END DO
*     --16/05/16 15:59-qishu-end------------------------------*
         END IF
*       Sort square distances with respect to the centre C.
         CALL SORT1(NP,R2,JLIST)
*       Sort square distances of massless particles
         CALL SORT1(Nnomass,Rnomass,NLIST)
c      END IF

*     escape calculation if only rscale is needed
      IF (KZ(7).EQ.1) GO TO 100
*     
*  Determine the Lagrangian radii for specified mass fractions.
*     NLAGR = Lagrangian radius for massless particles
*     NPARTC = particle counter within a shell
*     RC = Core radius (calculated in core.f)

*      wait for core.f update
      RCn = 0.0 
      I = 0
      ZM = 0.D0
*
      DO 15 J = 1,NLENS
         IF(KZ(7).EQ.2.OR.KZ(7).EQ.3.OR.J.EQ.1) THEN
            NPARTC(J) = 0
         ELSE
            NPARTC(J) = NPARTC(J-1)
         END IF
*
 20      I = I + 1
*     If reach last particle finish the loop
         IF(I.GT.Nnomass) THEN
            IF(NPARTC(J).NE.0) THEN
               II = I - 1
 31            IM = NLIST(II)
               IF(BODY(IM).EQ.0.0D0) THEN
                  II = II - 1
                  GO TO 31
               ELSE
                  NLAGR(J) = SQRT(Rnomass(II))
                  I = I - 1
                  GO TO 15
               END IF
            ELSE
               NLAGR(J) = NLAGR(J-1)
               I = I - 1
               GO TO 15
            END IF
         END IF
*     Sorted list index
         IM = NLIST(I)
*     escape ghost particle, but add it into particle count
         IF(BODY(IM).EQ.0.0D0) THEN
            NPARTC(J) = NPARTC(J) + 1
            GO TO 20
         END IF
*     Cumulated mass
         ZM = ZM + BODY(IM)
         NPARTC(J) = NPARTC(J) + 1
*
*     Check whether mass within Langrangian radius is complete.
         IF (I.LT.Nnomass.AND.ZM.LT.FLAGR(J)*ZMASS0) GO TO 20

*     Get average within a shell
         RLAGR(J) = SQRT(R2(I))
         NLAGR(J) = SQRT(Rnomass(I))
 15   CONTINUE


      IF (KZ(8).GT.0) then
c      stop "lagrN: KZ(8) = 0 not finished"
*  Determine the Lagrangian radii for singles and binaries separately.
         ZMS = 0.0D0
         I = 0
*
         DO 17 J = 1, NLENS
 21         I = I + 1
            IF(I.GT.NSNGL) THEN
               RSLAGR(J) = RSLAGR(J-1)
               GO TO 17
            END IF
            IM = ISLIST(I)
            IF(BODY(IM).EQ.0.0D0) GO TO 21
            ZMS = ZMS + BODY(IM)
            IF (I.LT.NSNGL.AND.ZMS.LT.FLAGR(J)*ZSMASS0) GO TO 21
            RSLAGR(J) = SQRT(RSNGL(I))
 17      CONTINUE
*
         ZMB = 0.0D0
         I = 0
*     
         DO 18 J = 1, NLENS
 22         I = I + 1
            IF(I.GT.NBIN) THEN
               RBLAGR(J) = RBLAGR(J-1)
               GO TO 18
            END IF
            IM = IBLIST(I)
            IF(BODY(IM).EQ.0.0D0) GO TO 22
            ZMB = ZMB + BODY(IM)
            IF (I.LT.NBIN.AND.ZMB.LT.FLAGR(J)*ZBMASS0) GO TO 22
            RBLAGR(J) = SQRT(RBIN(I))
 18      CONTINUE
*
      END IF
*
*     Write on diagnostics of NLAGR
      if(rank.eq.0)then
         IF (KZ(7).GE.2) THEN
            WRITE (6,40) (FLAGR(K),K=1,NLENS)
 40         FORMAT (/,11X,'TIME  M/MT(massless):',
     &      1P,18(1X,D9.2),2X,'<RC(massless)')
            WRITE (6,41) TTOT, (NLAGR(K),K=1,NLENS),RCn
 41         FORMAT (3X,D12.4,'           RLAGR:',1P,19(1X,D9.2))


            WRITE (6,43) TTOT, (NPARTC(K),K=1,NLENS)
 43         FORMAT (3X,D12.4,'          NPARTC:',18I10)
         END IF
      end if
*
*
 100  RETURN
*
      END
