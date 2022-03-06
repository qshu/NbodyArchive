


      SUBROUTINE BODIES
*
*
*       Output of single bodies or binaries.
*       ------------------------------------
*
      INCLUDE 'common1.h'
      integer k, ibody
      integer min, nwrite, ifreq
      integer i,j,jmin,nprint, nout, iopen
      real*8 fi,ei,rij2,ri,smax,rjmin2,a1,a2,a3,rijmin
      real * 8 vr2,erel,semi,zn,rdot,ecc2,ecc
      character * 132 name
      character * 132 formstr
      data iopen,nprint /0,0/
      nout = time/deltat + 0.5
*
*
*       Check option for printing single bodies.
      IF (KZ(9).EQ.0) GO TO 20
      K = KZ(9)
c      IBODY = MIN(5**K,N)
      IBODY  = N
*
      if(kz(8).eq. 0) kz(8) = 1
      if(kz(9) .ge. 99 .and. mod(nout,kz(8)) .eq. 0) then
          if(iopen .eq. 0) then
             write(6,*)' BODIES in file name'
             read(5,900) name
900          format(A)
             if (kz(9) .eq. 99) then
                formstr= 'formatted'
             else
                formstr= 'unformatted'
             endif
             open(unit=8,file=name,status='new',form=formstr)
             iopen = 1
          endif
          if(kz(9) .eq. 99) then
             write(8,*) N
             write(8,*) 3
             write(8,800) time
 800         format(g21.13)
             do 1000 i = 1, n
                write(8,800) body(i)
 1000           continue
                DO 1010 I = 1,N
                   write (8,810)  (X(K,I),K=1,3)
 810               format(3g21.13)
 1010           CONTINUE
                DO 1020 I = 1,N
                   write (8,810)  (XDOT(K,I),K=1,3)
 1020           CONTINUE
                do 1030 i = 1, n
                   write(8,800) phi(i)
 1030           continue
             else
                ifreq = max(1,kz(9) - 100)
                nwrite = nbh + (n-nbh-1)/ifreq + 1
c                write(8)'BHNBODY1 OUTPUT HEADER'
                write(8)nwrite,n,nbh,ifreq
                write(8) time
                write(8)(body(i),i=1,nbh),(body(i),i=nbh+1,n,ifreq)
                write(8)((x(k,i),k=1,3),i=1,nbh),
     $               ((x(k,i),k=1,3),i=nbh+1,n,ifreq)
                write(8)((xdot(k,i),k=1,3),i=1,nbh),
     $               ((xdot(k,i),k=1,3),i=nbh+1,n,ifreq)
                write(8)(phi(i),i=1,nbh),(body(i),i=nbh+1,n,ifreq)
             endif
             call flush(8)
          goto 20
      endif
      IF (KZ(9).ge.99) GO TO 20
      DO 10 I = 1,IBODY
          FI = 2.0*SQRT(F(1,I)**2 + F(2,I)**2 + F(3,I)**2)
          EI = 0.5*(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
          DO 4 J = 1,N
              IF (J.EQ.I) GO TO 4
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
              EI = EI - BODY(J)/SQRT(RIJ2 + EPS2)
    4     CONTINUE
          RI = SQRT(X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
          WRITE (6,6)  I, BODY(I), STEP(I), EI, RI, (X(K,I),K=1,3),
     &                 (XDOT(K,I),K=1,3)
    6     FORMAT (I7,F7.3,F8.4,F6.1,F6.2,1X,3F7.2,1X,3F6.2,F6.1)
   10 CONTINUE
*
*       Optional search for binaries (frequency NFIX with KZ(6) = 2).
   20 IF (KZ(6).EQ.0) GO TO 50
      IF (KZ(6).EQ.2.AND.NPRINT.NE.1) GO TO 50
      SMAX = 0.02*TCR
*
      DO 40 I = 1,N
          IF (STEP(I).GT.SMAX) GO TO 40
          JMIN = 0
          RJMIN2 = RSCALE**2
          DO 30 J = 1,N
              IF (STEP(J).GT.SMAX.OR.J.EQ.I) GO TO 30
              A1 = X(1,I) - X(1,J)
              A2 = X(2,I) - X(2,J)
              A3 = X(3,I) - X(3,J)
              RIJ2 = A1**2 + A2**2 + A3**2 + EPS2
              IF (RIJ2.LT.RJMIN2) THEN
                  RJMIN2 = RIJ2
                  JMIN = J
              END IF
   30     CONTINUE
          IF (JMIN.LT.I) GO TO 40
          RIJMIN = SQRT(RJMIN2)
          VR2 = (XDOT(1,I) - XDOT(1,JMIN))**2 +
     &          (XDOT(2,I) - XDOT(2,JMIN))**2 +
     &          (XDOT(3,I) - XDOT(3,JMIN))**2
          EREL = 0.5*VR2 - (BODY(I) + BODY(JMIN))/RIJMIN
          IF (EREL.GT.0.0) GO TO 40
          SEMI = -0.5*(BODY(I) + BODY(JMIN))/EREL
*
*       Only print significant binaries.
          IF (SEMI.GT.0.2*RSCALE) GO TO 40
          ZN = SQRT((BODY(I) + BODY(JMIN))/SEMI**3)
          RDOT = (X(1,I) - X(1,JMIN))*(XDOT(1,I) - XDOT(1,JMIN)) +
     &           (X(2,I) - X(2,JMIN))*(XDOT(2,I) - XDOT(2,JMIN)) +
     &           (X(3,I) - X(3,JMIN))*(XDOT(3,I) - XDOT(3,JMIN))
          ECC2 = (1.0 - RIJMIN/SEMI)**2 +
     &                             RDOT**2/(SEMI*(BODY(I) + BODY(JMIN)))
          ECC = SQRT(ECC2)
          RI = SQRT(X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
          WRITE (6,25)  I, JMIN, BODY(I), BODY(JMIN), EREL, SEMI, ZN,
     &                  RIJMIN, RI, ECC
   25     FORMAT (3X,'BINARY ',2I5,2F7.3,F6.1,F8.4,F7.1,F8.4,F6.2,F7.3)
   40 CONTINUE
*
   50 continue
      RETURN
*
      END


