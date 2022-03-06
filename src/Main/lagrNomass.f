      SUBROUTINE lagrNomass(C)
*
*
*       Lagrangian radii.
*       -----------------
*
      parameter (NLENS=18)
      INCLUDE 'common6.h'
      INCLUDE 'omp_lib.h'
      INTEGER NPARTN(NLENS), NCORE
      REAL*8  R2(NMAX),RSNGL(NMAX),RBIN(NMAX),C(3)
      REAL*8  FLAGR(NLENS),RLAGRN(NLENS)
      REAL*8  RSLAGR(NLENS),RBLAGR(NLENS)
      REAL*8  SIGR2(NLENS),SIGT2(NLENS),SIG2(NLENS)
      REAL*8  SIGR2C,SIGT2C,SIG2C
      REAL*8  VROT(NLENS),VROTC
      REAL*8  VRAVE(NLENS),VTAVE(3,NLENS),VAVE(3,NLENS)
      REAL*8  VTAVEV(NLENS),VAVEV(NLENS)
      REAL*8  VRI(NMAX), VTI(3,NMAX),VTITMP(3),VITMP(3)
      REAL*8  VR_CORE,VT_CORE(3),V_CORE(3),V_COREV,VT_COREV
      INTEGER ISLIST(NMAX),IBLIST(NMAX)
*
*     Lagrangian radii fraction of total mass
      DATA FLAGR/0.001D0,0.003D0,0.005D0,0.01D0,0.03D0,0.05D0,0.1D0,
     &     0.2D0,0.3D0,0.4D0,0.5D0,0.6D0,0.7D0,0.8D0,0.9D0,0.95D0,
     &     0.99D0,1.0D0/
*
*     Get total number of massless particles
      NP = 0
      DO I = IFIRST,N
         IF(nomass(I).eq.1) NP = NP + 1
      END DO
*
      IF(KZ(7).EQ.2.OR.KZ(7).EQ.4) THEN
*     Get initial total massless particle number
         NPART0 = NMASS0
      ELSE
         NPART0 = NP
      END IF
*
*     Particle number counts
      NP = 0
*
*     Set square radii of massless particles
!$omp parallel do private(I)
        DO I = 1,N
        IF(nomass(I).EQ.1.AND.BODY(I).GT.smallMass)
     &      STOP " Too massive massless particle "
           IF(BODY(I).LE.smallMass) THEN
               NP = NP + 1
               R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &              (X(3,I) - C(3))**2
               JLIST(NP) = I
           END IF
        END DO
!$omp end parallel do
*     Sort square distances of massless particles
         CALL SORT1(NP,R2,JLIST)
*     
*  Determine the Lagrangian radii for specified particle fractions.
*     RLAGRN = Lagrangian radius for massless particles
*     NPARTN = particle counter within a shell
*     NCORE = particle counter for the core
*     RC = Core radius (calculated in core.f)
*     VAVE = average velocity within a shell
*     VTAVE = average tangential velocity within a shell
*     VRAVE = average radial velocity within a shell
*     SIG2 = velocity dispersion square within a shell
*     SIGR2 = radius velocity dispersion square within a shell
*     SIGT2 = tangential velocity dispersion square within a shell
*     VROT = average rotational velocity projected in x-y plane within a shell
      VR_CORE = 0.0D0
      VT_CORE(1:3) = 0.0D0
      V_CORE(1:3) = 0.0D0
      SIGR2C = 0.0D0
      SIGT2C = 0.0D0
      SIG2C =0.0D0
      VROTC = 0.0D0
      RC2 = RC*RC
      NCORE = 0
      NPART = 0
      I = 0
*
      DO 15 J = 1,NLENS
         IF(KZ(7).EQ.2.OR.KZ(7).EQ.3.OR.J.EQ.1) THEN
            VROT(J) = 0.0D0
            VRAVE(J) = 0.0D0
            VAVE(1:3,J) = 0.0D0
            VTAVE(1:3,J) = 0.0D0
            NPARTN(J) = 0
         ELSE
            VROT(J) = VROT(J-1)
            VRAVE(J) = VRAVE(J-1)
            VAVE(1:3,J) = VAVE(1:3,J-1)
            VTAVE(1:3,J) = VTAVE(1:3,J-1)
            NPARTN(J) = NPARTN(J-1)
         END IF
*
 20      I = I + 1
*     If reach last particle finish the loop
         IF(I.GT.NP) THEN
            IF(NPARTN(J).NE.0) THEN
               II = I - 1
 31            IM = JLIST(II)
               IF(BODY(IM).EQ.0.0D0) THEN
                  II = II - 1
                  GO TO 31
               ELSE
                  RLAGRN(J) = SQRT(R2(II))
                  I = I - 1
                  GO TO 15
               END IF
            ELSE
               RLAGRN(J) = RLAGRN(J-1)
               I = I - 1
               GO TO 15
            END IF
         END IF
*     Sorted list index
         IM = JLIST(I)
*     escape ghost particle, but add it into particle count
         IF(BODY(IM).EQ.0.0D0) THEN
            NPARTN(J) = NPARTN(J) + 1
            GO TO 20
         END IF
*     Cumulated velocity (particle average)
         DO K = 1,3
            VAVE(K,J) = VAVE(K,J) + XDOT(K,IM)
         END DO
*
*     Radial velocity can be calculated by velocity .dot. radial direction vector
         VRI(I) = 0.D0
         DO K = 1,3
            VRI(I) = VRI(I) + XDOT(K,IM)*(X(K,IM)-C(K))/DSQRT(R2(I))
         END DO
*     Cumulate radial velocity (particle average)
         VRAVE(J) = VRAVE(J) + VRI(I)
*     Radial velocity square
         VR2I = VRI(I)*VRI(I)
*     Tangential velocity
         DO K = 1,3
            VTI(K,I) = XDOT(K,IM) - VRI(I)*(X(K,IM)-C(K))/DSQRT(R2(I))
*     Cumulate tangential velocity (particle average)
            VTAVE(K,J) = VTAVE(K,J) + VTI(K,I)
         END DO
*     X-Y plane projected position square
         RR12 = X(1,IM)**2 + X(2,IM)**2
*     X-Y plane radial velocity value * sqrt(RR12) 
         XR12 = 0.D0
         DO K = 1,2
            XR12 = XR12 + XDOT(K,IM)*(X(K,IM)-C(K))
         END DO
*     Rotational velocity
         VROT1 = XDOT(1,IM) - XR12/RR12*X(1,IM)
         VROT2 = XDOT(2,IM) - XR12/RR12*X(2,IM)
*     Rotational direction sign
         XSIGN = VROT1*X(2,IM)/DSQRT(RR12) - VROT2*X(1,IM)/DSQRT(RR12)
         VROTM = DSQRT(VROT1**2+VROT2**2)
*     Cumulate rotational velocity
         IF(XSIGN.GT.0.D0) THEN
            VROT(J) = VROT(J) + VROTM
         ELSE
            VROT(J) = VROT(J) - VROTM
         END IF
*     Cumulate particle number in shell
         NPART = NPART + 1
         NPARTN(J) = NPARTN(J) + 1
*
*     Independent determination of particles in core radius.
         IF (R2(I).LT.RC2) THEN
            VR_CORE = VR_CORE + VRI(I)
            DO K = 1,3
               V_CORE(K) = V_CORE(K) + XDOT(K,IM)
               VT_CORE(K) = VT_CORE(K) + VTI(K,I)
            END DO
            IF(XSIGN.GT.0.D0) THEN
               VROTC = VROTC + VROTM
            ELSE
               VROTC = VROTC - VROTM
            END IF
            NCORE = NCORE + 1
         END IF
*
*     Check whether mass within Langrangian radius is complete.
      IF (I.LT.NP.AND.FLOAT(NPART).LT.FLAGR(J)*FLOAT(NPART0))GO TO 20
*
*     Get average within a shell
         RLAGRN(J) = SQRT(R2(I))
         VRAVE(J) = VRAVE(J)/FLOAT(NPARTN(J))
         VROT(J) = VROT(J)/FLOAT(NPARTN(J))
         DO K = 1,3
            VTAVE(K,J) = VTAVE(K,J)/FLOAT(NPARTN(J))
            VAVE(K,J) = VAVE(K,J)/FLOAT(NPARTN(J))
         END DO
 15   CONTINUE
*
*     Get average within core radius
      VR_CORE = VR_CORE/FLOAT(NCORE)
      V_COREV = 0.0D0
      VT_COREV = 0.0D0
      DO K = 1,3
         VT_CORE(K) = VT_CORE(K)/FLOAT(NCORE)
         VT_COREV = VT_COREV + VT_CORE(K)*VT_CORE(K)
         V_CORE(K) = V_CORE(K)/FLOAT(NCORE)
         V_COREV = V_COREV + V_CORE(K)*V_CORE(K)
      END DO
      VROTC = VROTC/FLOAT(NCORE)
*
*     Get velocity dispersion
      NSTART = 1
      NSHELL = 0
*
      DO 16 J = 1,NLENS
         NSTART = NSHELL + 1
         IF(KZ(7).EQ.2.OR.KZ(7).EQ.3.OR.J.EQ.1) THEN
            SIGR2(J) = 0.0D0
            SIGT2(J) = 0.0D0
            SIG2(J) = 0.0D0
            NSHELL = NSHELL + NPARTN(J)
         ELSE
            SIGR2(J) = SIGR2(J-1)
            SIGT2(J) = SIGT2(J-1)
            SIG2(J) = SIG2(J-1)
            NSHELL = NPARTN(J)
         END IF
*
         DO IK = NSTART, NSHELL
            IM = JLIST(IK)
*     Radial direction
            VR2I = VRI(IK) - VRAVE(J)
            VR2I = VR2I*VR2I
*     Tangential direction and original direction
            VTI2 = 0.0D0
            VI2 = 0.0D0
            DO K = 1,3
               VTITMP(K) = VTI(K,IK) - VTAVE(K,J)
               VTI2 = VTI2 + VTITMP(K)*VTITMP(K)
               VITMP(K) = XDOT(K,IM) - VAVE(K,J)
               VI2 = VI2 + VITMP(K)*VITMP(K)
            END DO
*     Cumulate radial velocity dispersion
            SIGR2(J) = SIGR2(J) + VR2I
*     Cumulate tangential velocity dispersion
            SIGT2(J) = SIGT2(J) + VTI2/2.D0
*     Cumulate velocity dispersion
            SIG2(J) = SIG2(J) + VI2/3.D0
*     Only for core radius
            IF (R2(IK).LT.RC2) THEN
               VR2I = VRI(IK) - VR_CORE
               VR2I = VR2I*VR2I
               SIGR2C = SIGR2C + VR2I

               VTI2 = 0.0D0
               VI2 = 0.0D0
               DO K = 1,3
                  VTITMP(K) = VTI(K,IK) - VT_CORE(K)
                  VTI2 = VTI2 + VTITMP(K)*VTITMP(K)
                  VITMP(K) = XDOT(K,IM) - V_CORE(K)
                  VI2 = VI2 + VITMP(K)*VITMP(K)
               END DO
               SIGT2C = SIGT2C + VTI2/2.D0
               SIG2C = SIG2C + VI2/3.D0
            END IF
         END DO
 16   CONTINUE
*
*     Average velocity dispersion for core region
      SIGR2C = SIGR2C/FLOAT(NCORE)
      SIGT2C = SIGT2C/FLOAT(NCORE)
      SIG2C = SIG2C/FLOAT(NCORE)
*
*     Final average
      DO J = 1, NLENS
         IF(NPARTN(J).GT.0) THEN
            VTAVEV(J) = 0.0D0
            VAVEV(J) = 0.0D0
            DO K = 1,3
               VTAVEV(J) = VTAVEV(J) + VTAVE(K,J)*VTAVE(K,J)
               VAVEV(J) = VAVEV(J) + VAVE(K,J)*VAVE(K,J)
            END DO
            VTAVEV(J) = SQRT(VTAVEV(J))
            VAVEV(J) = SQRT(VAVEV(J))
            SIGR2(J) = SIGR2(J)/FLOAT(NPARTN(J))
            SIGT2(J) = SIGT2(J)/FLOAT(NPARTN(J))
            SIG2(J) = SIG2(J)/FLOAT(NPARTN(J))
         END IF
      END DO
*     Write on diagnostics of RLAGRN
      if(rank.eq.0)then
         IF (KZ(7).GE.2) THEN
            WRITE (6,40) (FLAGR(K),K=1,NLENS)
 40         FORMAT (/,'(massless) TIME   N/NT:',1P,18(1X,D9.2),2X,'<RC')
            WRITE (6,41) TTOT, (RLAGRN(K),K=1,NLENS),RC
 41         FORMAT (3X,D12.4,' RLAGRN:',1P,19(1X,D9.2))
            WRITE (6,43) TTOT, (NPARTN(K),K=1,NLENS),NCORE
 43         FORMAT (3X,D12.4,' NPARTN:',19I10)
            WRITE (6,45) TTOT, (SIGR2(K),K=1,NLENS),SIGR2C
 45         FORMAT (3X,E12.4,' SIGR2N:',1P,19(1X,E9.2))
            WRITE (6,46) TTOT, (SIGT2(K),K=1,NLENS),SIGT2C
 46         FORMAT (3X,E12.4,' SIGT2N:',1P,19(1X,E9.2))
            WRITE (6,47) TTOT, (VROT(K),K=1,NLENS),VROTC
 47         FORMAT (3X,E12.4,' VROTN: ',1P,19(1X,E9.2))
         END IF
      end if
*
*
 100  RETURN
*
      END
