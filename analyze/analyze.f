
C***********************************************************************
C
C
                        PROGRAM analyze
C
C
C***********************************************************************
C
C
C     Subroutine to analyze the ellipticity of the system
C
C
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)

      PARAMETER  (NMAX=150100)
      PARAMETER  (ID=3)
*
*
*       ------------------------------------------------------
*       NMAX    Maximum number of single bodies & c.m.
*       KMAX    Maximum number of KS solutions.
*       LMAX    Maximum size of neighbour lists.
*       MMAX    Maximum number of merged binaries.
*       MLD     Maximum number of disrupted KS components.
*       MLR     Maximum number of recently regularized KS pairs.
*       MLV     Maximum number of high-velocity particles.
*       MCL     Maximum number of interstellar clouds.
*       NCMAX   Maximum members of chain members.
*       ------------------------------------------------------
*       ID      First dimension of 3d-vectors (usually = 3)
*               Has to be set to 4 on CRAY T3D for hardware reasons
      COMMON/NBODY/BODY(NMAX),X(ID,NMAX),XDOT(ID,NMAX),PHIDBL(NMAX),
     *  RI(NMAX),VI(NMAX)
      COMMON/ROT/ERG(NMAX),XKIN(NMAX),WPT(NMAX),ROT(NMAX),RTRUE(NMAX),
     *     WPTX(NMAX)
*
C   Declaration of local variables.
C   -------------------------------
        DOUBLE PRECISION etoti(1:nmax),mtot,poscm(3),velcm(3),
     &     ti(3,3),tiwork(3,3),dwork(3),ework(3),lam(3)
        INTEGER i,j,k,ief,nbound,nstart,nnext,np,indexev(3)
        INTEGER index(1:nmax)

        INTEGER nef,nef1
        PARAMETER(nef=9,nef1=nef+1)
        DOUBLE PRECISION xf(nef1),ba(nef1),ca(nef1),taue(nef1),
     &            evec(3,3,nef1)
        DOUBLE PRECISION VRAD(3),VP(3),VZ(3)

        PARAMETER (tiny=1.D-30)
        DATA (xf(i),i=1,nef1) /0.01D0,0.05D0,0.1D0,0.2D0,0.5D0,
     &                        0.8D0,0.9D0,0.95D0,1.0D0,99.D0/

C=======================================================================
*       Read particle data
          DO 5 I = 1,NMAX
        READ (10,*,END=555)BODY(I),(X(K,I),K=1,3),(XDOT(K,I),K=1,3)
    5     CONTINUE
*       Read tidal radius if cutoff required
 555    N = I - 2 
*
         ifirst = 1
         ntot = n
          PRINT*,N,' body data read from unit 10 '
*
        REWIND 10
          DO 51 I = 1,N
        READ (10,*)DUMMY
        READ (11,*,END=9999)ERG(I),XKIN(I),WPT(I),ROT(I),RTRUE(I),
     * WPTX(I)
   51     CONTINUE
*
        READ(10,*)RTID,XETID,BETF,YKIN,YP,YP2,YPOT
*
        YPOT = 0.5D0*YPOT
*
        PRINT*,' RTID,1/RTID,XETID=',RTID,1.D0/RTID,XETID
        SXVIR = YKIN/ABS(YPOT)
        PRINT*,' Virial ratio from data=',SXVIR
        PRINT*,' YKIN,YPOT=',YKIN,YPOT
      ZMASS = 0.D0
      DO 45 I = 1,N
   45     ZMASS = ZMASS + BODY(I)
      PRINT*,' ZMASS before scaling=',ZMASS
*       Scale masses to standard units of <M> = 1/N.
      DO 50 I = 1,N
          BODY(I) = BODY(I)/ZMASS
   50 CONTINUE
*
      ZMASS = 1.D0
      IPASS = 0
*
 1111 CONTINUE
* 
*       Calculate the potential energy.
      ZKIN = 0.D00
      POT = 0.D0
*
      DO 20 I = 1,N
      JMIN = I + 1
      POTJ = 0.D00
      POTI = 0.D00
*       POTI contains potential at particles position to be stored later (R.Sp.)
*
      DO 30 J = 1,N
      IF (J.EQ.I)GOTO 30
          A1 = X(1,I) - X(1,J)
          A2 = X(2,I) - X(2,J)
          A3 = X(3,I) - X(3,J)
      A4 = BODY(J)/SQRT(A1*A1 + A2*A2 + A3*A3)
      POTI = POTI - A4
*  also J.LT.N?
      IF(J.GE.JMIN)POTJ = POTJ + A4
   30 CONTINUE
*       Store potential in shared vector first (R.Sp.)
      PHIDBL(I) = POTI
      POT = POT + BODY(I)*POTJ
      ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2+XDOT(2,I)**2+XDOT(3,I)**2)
   20 CONTINUE
*
      ZKIN = 0.5D0*ZKIN
*
      IF(IPASS.EQ.0)THEN
      SXSC = SQRT(0.5D0*ABS(POT)/ZKIN)
*     SXSC = SQRT(YKIN/ABS(YPOT)*ABS(POT)/ZKIN)
      PRINT*,' Scaling factor of velocities=',SXSC
      PRINT*,' YKIN,YPOT,POT,ZKIN=',YKIN,YPOT,POT,ZKIN,IPASS
      END IF
*
      IF(IPASS.EQ.1)THEN
* Scale the kinetic energy to desired virial ratio
      DO 55 I=1,N
      DO 56 K=1,3
 56   XDOT(K,I)=XDOT(K,I)*SXSC
 55   CONTINUE
      END IF
*
      REWIND 20
*       Sum the kinetic energy (include c.m. bodies but not components).
      ZKIN = 0.D0
      EKR = 0.D0
      EKZ = 0.D0
      EKP = 0.D0
      DO 40 I = 1,N
          RI2 = X(1,I)**2+X(2,I)**2
          RI(I) = DSQRT(RI2+X(3,I)**2)
          EROT = 0.5D0*XDOT(2,I)**2
          EXKIN = XDOT(1,I)**2+XDOT(2,I)**2+XDOT(3,I)**2
          VDOTR = X(1,I)*XDOT(1,I)+X(2,I)*XDOT(2,I)
*
      DO 47 K = 1,2
 47       VRAD(K) = VDOTR*X(K,I)/RI2
          VRAD(3) = 0.D0
          VZ(1) = 0.D0
          VZ(2) = 0.D0
          VZ(3) = XDOT(3,I)
      DO 46 K = 1,2
 46       VP(K) = XDOT(K,I) - VRAD(K) - VZ(K)
          VP(3) = 0.D0
*
       EKR = EKR + BODY(I)*(VRAD(1)**2+VRAD(2)**2)
       EKP = EKP + BODY(I)*(VP(1)**2+VP(2)**2)
       EKZ = EKZ + BODY(I)*(VZ(3)**2)
*
          VI(I) = DSQRT(EXKIN)
          EXKIN = 0.5D0*EXKIN
          ETOT = EXKIN + PHIDBL(I)
      WRITE(20,120)I,ERG(I),ETOT,XKIN(I),EXKIN,WPT(I),PHIDBL(I),
     *   ROT(I),EROT,RTRUE(I),RI(I),WPTX(I)
 120  FORMAT(1X,I5,1P,11(E12.3))
*
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
*
   40 CONTINUE
*
          ZKIN = ZKIN *0.5D0
          EKR = EKR * 0.5D0
          EKP = EKP * 0.5D0
          EKZ = EKZ * 0.5D0
*
       PRINT*,' Kin E-Comp R,P,Z=',EKR,EKP,EKZ,' Sum=',EKR+EKP+EKZ
C       sort for particle radius
C       ------------------------

*       CALL indexx(n,ri,index)
*
*     do 557 i=1,20
*     ipo = index(i)
*     print*,' i=',ipo,' pos=',ri(ipo),' phi=',phidbl(ipo),
*    *    ' v=',vi(ipo)
*557  continue
*     do 556 i=n-20,n
*     ipo = index(i)
*     print*,' i=',ipo,' pos=',ri(ipo),' phi=',phidbl(ipo),
*    *    ' v=',vi(ipo)
*556  continue
*
      ETOT = ZKIN - POT
*
      PRINT*,' IPASS=',IPASS,' ETOT, EKIN, POT=',ETOT,ZKIN,-POT
*
      IPASS = IPASS + 1
*
      IF(IPASS.EQ.1)GOTO 1111
*
      E0 = -0.25D0
      SX = E0/ETOT
*
      PRINT*,' Scale with SX=',SX
*
*       Scale coordinates & velocities to the new units.
      DO 70 I = 1,N
          DO 68 K = 1,3
              X(K,I) = X(K,I)/SX
              XDOT(K,I) = XDOT(K,I)*SQRT(SX)
   68     CONTINUE
   70 CONTINUE
      RTID = RTID/SX
*
      IF(IPASS.LE.2)GOTO 1111
*
      rewind 30
      do 777 il=1,n
      write(30,130)body(il),(x(k,il),k=1,3),(xdot(k,il),k=1,3)
  130 FORMAT (7(1pE11.3))
 777  continue
      write(30,131)rtid
 131  FORMAT(1pE11.3)
*
C       calculate specific energy of particles
C       --------------------------------------

        DO 100 i=ifirst,ntot
        etoti(i-ifirst+1) = 0.5D0 * (xdot(1,i)**2 + xdot(2,i)**2 +
     &                         xdot(3,i)**2) + phidbl(i)
100     CONTINUE

C      calculate number of bound particles
C      -----------------------------------

        nbound = 0
        DO 150 i=1,ntot-ifirst+1
           IF(etoti(i).LT.0.D0) nbound = nbound + 1
150     CONTINUE

C       sort for particle energy
C       ------------------------

        CALL indexx(ntot-ifirst+1,etoti,index)

C       initialize tensor of inertia
C       ----------------------------

        DO 210 i=1,3
           DO 200 k=1,3
              ti(i,k) = 0.D0
200        CONTINUE
210     CONTINUE

C       LOOP over fraction of most bound particles and all particles
C       ------------------------------------------------------------

        nstart   = 1
        mtot     = 0.D0
        poscm(1) = 0.D0
        poscm(2) = 0.D0
        poscm(3) = 0.D0
        DO 500 ief=1,nef1

           IF(ief.LE.nef) THEN
C                                  only fraction of bound particles
C                                  --------------------------------
              nnext = MAX(ione,NINT(xf(ief) * nbound))
           ELSE
C                                   all particles
C                                   -------------
              nnext = ntot
           ENDIF

C-----------------------------------------------------------------
C--      at least two particles are required for ellipticity...
C-----------------------------------------------------------------
           IF(nnext.LT.2) THEN
              ba(ief)  = 999.
              ca(ief)  = 999.
              taue(ief) = 999.
              DO 320 k=1,3
                 DO 310 j=1,3
                    evec(k,j,ief) = 0.
310              CONTINUE
320           CONTINUE

           ELSE

C       calculate tensor of inertia
C       ---------------------------
*             print*,' nstart,nnext=',nstart,nnext
              DO 400 i=nstart,nnext
                 ipo = index(i) + ifirst - 1
                 ti(1,1) = ti(1,1) + body(ipo) *
     &                  (x(2,ipo)*x(2,ipo) + x(3,ipo)*x(3,ipo))
                 ti(2,2) = ti(2,2) + body(ipo) *
     &                  (x(1,ipo)*x(1,ipo) + x(3,ipo)*x(3,ipo))
                 ti(3,3) = ti(3,3) + body(ipo) *
     &                  (x(1,ipo)*x(1,ipo) + x(2,ipo)*x(2,ipo))
                 ti(1,2) = ti(1,2) - body(ipo) * x(1,ipo)*x(2,ipo)
                 ti(1,3) = ti(1,3) - body(ipo) * x(1,ipo)*x(3,ipo)
                 ti(2,3) = ti(2,3) - body(ipo) * x(2,ipo)*x(3,ipo)
400           CONTINUE
      
C       correct for center of mass
C       --------------------------

C   A) calculate center of mass data
C   --------------------------------

C--       remove previous correction for center of mass
              ti(1,1) = ti(1,1) + mtot * (poscm(2)**2+poscm(3)**2)
              ti(2,2) = ti(2,2) + mtot * (poscm(1)**2+poscm(3)**2)
              ti(3,3) = ti(3,3) + mtot * (poscm(1)**2+poscm(2)**2)
              ti(1,2) = ti(1,2) - mtot * poscm(1) * poscm(2)
              ti(1,3) = ti(1,3) - mtot * poscm(1) * poscm(3)
              ti(2,3) = ti(2,3) - mtot * poscm(2) * poscm(3)
              poscm(1) = poscm(1) * mtot
              poscm(2) = poscm(2) * mtot
              poscm(3) = poscm(3) * mtot
*             xav = 0.d0
*             yav = 0.d0
*             zav = 0.d0

              DO 405 i=nstart,nnext
                 ipo = index(i) + ifirst - 1
                 poscm(1) = poscm(1) + body(ipo) * x(1,ipo)
                 poscm(2) = poscm(2) + body(ipo) * x(2,ipo)
                 poscm(3) = poscm(3) + body(ipo) * x(3,ipo)
                 mtot     = mtot + body(ipo)
*                xav = xav + abs(x(1,ipo))
*                yav = yav + abs(x(2,ipo))
*                zav = zav + abs(x(3,ipo))

 405           CONTINUE
              poscm(1) = poscm(1) / mtot
              poscm(2) = poscm(2) / mtot
              poscm(3) = poscm(3) / mtot
*             print*,' av=',xav,yav,zav

 
C   B) apply correction
C   -------------------
              ti(1,1) = ti(1,1) - mtot * (poscm(2)**2+poscm(3)**2)
              ti(2,2) = ti(2,2) - mtot * (poscm(1)**2+poscm(3)**2)
              ti(3,3) = ti(3,3) - mtot * (poscm(1)**2+poscm(2)**2)
              ti(1,2) = ti(1,2) + mtot * poscm(1) * poscm(2)
              ti(1,3) = ti(1,3) + mtot * poscm(1) * poscm(3)
              ti(2,3) = ti(2,3) + mtot * poscm(2) * poscm(3)

C       set off-axis values by symmetry
C       -------------------------------

              ti(2,1) = ti(1,2)
              ti(3,1) = ti(1,3)
              ti(3,2) = ti(2,3)
*
*             print*,' mtot,poscm=',mtot,(poscm(k),k=1,3)
C=======================================================
C       determine eigenvalues and axis of inertia
C=======================================================

C------------------------------------------------------
C--          copy tensor of inertia
C------------------------------------------------------
              DO 420 i=1,3
                 DO 410 k=1,3
                    tiwork(i,k) = ti(i,k) 
410              CONTINUE
420           CONTINUE
              np = 3

C------------------------------------------------------
C--          calculate eigenvalues and eigenvectors
C------------------------------------------------------
              CALL tred2(tiwork,np,np,dwork,ework)
              CALL tqli(dwork,ework,np,np,tiwork)

C--               sort for increasing eigenvalues
              CALL indexx(np,dwork,indexev)
C--               find eigenvectors
              DO 450 i=1,np
                 lam(i) = dwork(indexev(i))
                 DO 440 k=1,np
                    evec(k,i,ief) = tiwork(k,indexev(i))
440              CONTINUE
450           CONTINUE

              xhelp    = lam(3) + lam(2) - lam(1)
              xhelp1   = lam(2) - lam(3) + lam(1)
c              IF(xhelp1.LT.0.D0) THEN
c                 PRINT*,' ellan: xhelp1 < 0',xhelp1,tnow
c                 xhelp1 = 0.D0
c             ENDIF
              ba(ief)  = SQRT(MAX(tiny,lam(3)-lam(2)+lam(1)) / xhelp)
              ca(ief)  = SQRT(MAX(tiny,xhelp1) / xhelp) 
              taue(ief) = (ba(ief)-ca(ief)) / MAX(tiny,(1.D0 - ca(ief)))

              nstart = nnext + 1

           ENDIF

500     CONTINUE

C==================================================================
C==         OUTPUT of data
C===================================================================

      	DO 600 ief=1,nef1
           WRITE(60,910) time,xf(ief),ba(ief),ca(ief),taue(ief),
     *     mtot,(poscm(k),k=1,3)
600     CONTINUE

900     FORMAT(1x,1p,e12.5,1x,0p,i9,1x,i9,1x,i2)
910     FORMAT(1x,14(f9.5,1x))

        R E T U R N
9999    PRINT*,' I=',I,' Error in Read 11'
        END

