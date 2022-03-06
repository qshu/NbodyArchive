      SUBROUTINE KICK(I,ICASE)
*
*
*       Velocity kick for neutron stars.
*       --------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  VK(4)
      REAL*4  RAN2
      SAVE  IPAIR, KC, VDIS, RI
      DATA  IPAIR,KC /0,0/
*
*
*       Save orbital parameters in case of KS binary (called from KSAPO).
      IF (ICASE.EQ.0) THEN
          IPAIR = I
          KC = KSTAR(N+IPAIR)
*       Identify the correct component (KSTAR reversed in MDOT).
          I1 = 2*IPAIR - 1
          IF (KSTAR(I1).LT.0) THEN
              IN = I1
          ELSE IF (KSTAR(I1+1).LT.0) THEN
              IN = I1 + 1
          END IF
          KSTAR(IN) = -KSTAR(IN)
*
*       Determine mass loss and actual disruption velocity.
          DM = BODY(IN) - 1.4/ZMBAR
          VD2 = 2.0*(BODY(N+IPAIR) - DM)/R(IPAIR)
          VDIS = SQRT(VD2)*VSTAR
*       Set cluster escape velocity (add twice central potential).
          VP2 = 2.0*BODY(N+IPAIR)/R(IPAIR)
          VESC = SQRT(VP2 + 4.0)*VSTAR
          SEMI = -0.5*BODY(N+IPAIR)/H(IPAIR)
          ZM1 = BODY(I1)*SMU
          ZM2 = BODY(I1+1)*SMU
          EB = BODY(I1)*BODY(I1+1)/BODY(N+IPAIR)*H(IPAIR)
          RI = R(IPAIR)
*       Sum whole binding energy (used by BINOUT for net change).
          EKICK = EKICK + EB
          I2 = I1 + 1
          WRITE (6,1)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2), ZM1,
     &                 ZM2, VESC, VDIS, R(IPAIR)/SEMI, EB, R(IPAIR)
    1     FORMAT (' BINARY KICK:    NAM K* M1 M2 VESC VDIS R/A EB R ',
     &                              2I6,2I4,4F7.1,F6.2,F9.4,1P,E10.2)
          NBKICK = NBKICK + 1
          GO TO 30
      END IF
*
*       Generate velocity kick for neutron star (Gordon Drukier Tokyo paper).
*     IT = 0
*     V0 = 330.0
*     VCUT = 1000.0
*   2 VT = VCUT/V0*RAN2(IDUM1)
*     VP = VT*(2.0*RAN2(IDUM1) - 1.0)
*     VN = SQRT(VT**2 - VP**2)
*     FAC = 1.0/0.847*VN**0.3/(1.0 + VN**3.3)
*     IF (FAC.LT.RAN2(IDUM1).AND.IT.LT.10) GO TO 2
*     VKICK = V0*VT
*
*       Adopt the Maxwellian of Hansen & Phinney (MN 291, 569, 1997).
      DISP = 190.0
*       Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
      DO 2 K = 1,2
          X1 = RAN2(IDUM1)
          X2 = RAN2(IDUM1)
*       Generate two velocities from polar coordinates S & THETA.
          S = DISP*SQRT(-2.0*LOG(1.0 - X1))
          THETA = TWOPI*X2
          VK(2*K-1) = S*COS(THETA)
          VK(2*K) = S*SIN(THETA)
    2 CONTINUE
      VKICK = SQRT(VK(1)**2 + VK(2)**2 + VK(3)**2)
      VK(4) = VKICK
*
*       Limit kick velocity to VDIS+10*VSTAR/10*VST for binary/single stars.
      IF (IPAIR.GT.0) THEN
          VBF = SQRT(VDIS**2 + 100.0*VSTAR**2)
          VKICK = MIN(VKICK,VBF)
      ELSE
          VKICK = MIN(VKICK,10.0D0*VSTAR)
      END IF
      VKICK = VKICK/VSTAR
*
*       Randomize the velocity components.
*     A(4) = 0.0
*     DO 5 K = 1,3
*         A(K) = 2.0*RAN2(IDUM1) - 1.0
*         A(4) = A(4) + A(K)**2
*   5 CONTINUE
*
*       Add truncated/full kick velocity and initialize X0DOT.
      VI2 = 0.0
      VF2 = 0.0
      DO 10 K = 1,3
          VI2 = VI2 + XDOT(K,I)**2
*         XDOT(K,I) = XDOT(K,I) + VKICK*A(K)/SQRT(A(4))
          XDOT(K,I) = XDOT(K,I) + VKICK*VK(K)/VK(4)
          X0DOT(K,I) = XDOT(K,I)
          VF2 = VF2 + XDOT(K,I)**2
   10 CONTINUE
*
*       Modify energy loss due to increased velocity of single particle.
      ECDOT = ECDOT - 0.5*BODY(I)*(VF2 - VI2)
      NKICK = NKICK + 1
*
*       Replace final velocity by relative velocity for binary kick.
      IF (IPAIR.GT.0) THEN
          JP = KVEC(I)
          J = I + 1
          IF (I.EQ.2*JP) J = I - 1
          VF2 = 0.0
          DO 15 K = 1,3
              VF2 = VF2 + (XDOT(K,I) - XDOT(K,J))**2
   15     CONTINUE
          HNEW = 0.5*VF2 - (BODY(I) + BODY(J))/RI
          EB1 = BODY(I)*BODY(J)/(BODY(I) + BODY(J))*HNEW
          IF (EB1.LT.0.0) EKICK = EKICK - EB1
          IPAIR = 0
      END IF
*
      IF (NKICK.LT.50.OR.NAME(I).LE.2*NBIN0) THEN
          ZM = BODY(I)*ZMBAR
          WRITE (6,20)  I, NAME(I), KSTAR(I), KC, BODY0(I)*ZMBAR, ZM,
     &                  SQRT(VI2)*VSTAR, VKICK*VSTAR, SQRT(VF2)*VSTAR
   20     FORMAT (' VELOCITY KICK:    I NAM K* KC* M0 M VI VK VF ',
     &                                2I6,2I4,2F7.2,3F7.1)
          KC = 0
      END IF
*
*       Highlight velocities below 4 times rms velocity.
      IF (VKICK.LT.4.0*SQRT(0.5)) THEN
          WRITE (6,25)  I, NAME(I), VKICK*VSTAR, SQRT(VF2)*VSTAR
   25     FORMAT (' LOW KICK:    I NAM VK VF ',2I6,2F6.2)
      END IF
*
*       Include optional list of high-velocity particles.
      IF (KZ(37).GT.0) THEN
          CALL HIVEL(I)
      END IF
*
   30 RETURN
*
      END
