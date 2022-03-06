      SUBROUTINE SWITCH(Y)
*
*       Switching of chain.
*       -------------------
*
      INCLUDE 'COMMON1.CH'
      INCLUDE 'COMMON2.CH'
      REAL*8 Y(1),XCNEW(NMX3)
      INTEGER IOLD(NMX)
*
*       Copy Y-array to COMMON.
      CALL YSAVE(Y)
*
*	First transform to chain coordinates.
      DO I=1,N-1
      L1=3*(I-1)+1
      KS1=4*(I-1)+1
      CALL KSPHYS(Q(KS1),P(KS1),XC(L1),WC(L1))
      END DO
*
      L2=3*(INAME(1)-1)
      DO K=1,3
      X(L2+K)=0.0
      END DO
*
*       Set X for determining new chain indices.
      DO I=1,N-1
      L=3*(I-1)
      L1=L2
      L2=3*(INAME(I+1)-1)
      DO K=1,3
      X(L2+K)=X(L1+K)+XC(L+K)
      END DO
      END DO
*
*	Save the old chain indices.
      DO I=1,N
      IOLD(I)=INAME(I)
      END DO
*
*       Select new indices.
      CALL SELECT
*
*	Transform chain momenta.
      L1=3*(IOLD(1)-1)
      LN=3*(IOLD(N)-1)
      L=3*(N-2)
      DO K=1,3
      PI(L1+K)=-WC(K)
      PI(LN+K)=WC(L+K)
      END DO
      DO I=2,N-1
      L=3*(I-1)
      LI=3*(IOLD(I)-1)
      DO K=1,3
      PI(LI+K)=WC(L+K-3)-WC(L+K)
      END DO
      END DO
      L1=3*(INAME(1)-1)
      LN=3*(INAME(N)-1)
      L=3*(N-2)
      DO K=1,3
      WC(K)=-PI(L1+K)
      WC(L+K)=PI(LN+K)
      END DO
      DO I=2,N-2
      L=3*(I-1)
      LI=3*(INAME(I)-1)
      DO K=1,3
      WC(L+K)=WC(L+K-3)-PI(LI+K)
      END DO
      END DO
*
*       Construct new chain coordinates.
*       Transformation matrix (old to new) has only coefficients -1, 0 or +1.
      DO I=1,3*(N-1)
      XCNEW(I)=0.0
      END DO
      DO ICNEW=1,N-1
*       Find K0 & K1 for iold(k0) = iname(icnew) & iold(k1) = iname(icnew+1).
      LNEW=3*(ICNEW-1)
      DO I=1,N
      IF(IOLD(I).EQ.INAME(ICNEW))K0=I
      IF(IOLD(I).EQ.INAME(ICNEW+1))K1=I
      END DO
      DO ICOLD=1,N-1
      LOLD=3*(ICOLD-1)
      IF( (K1.GT.ICOLD).AND.(K0.LE.ICOLD))THEN
*       Add.
      DO K=1,3
      XCNEW(LNEW+K)=XCNEW(LNEW+K)+XC(LOLD+K)
      END DO
      ELSEIF( (K1.LE.ICOLD).AND.(K0.GT.ICOLD) )THEN
*	Subtract.
      DO K=1,3
      XCNEW(LNEW+K)=XCNEW(LNEW+K)-XC(LOLD+K)
      END DO
      END IF
      END DO
      END DO
*
*	Perform KS-transformations.
      DO I=1,N-1
      L1=3*(I-1)+1
      KS1=4*(I-1)+1
      CALL PHYSKS(XCNEW(L1),WC(L1),Q(KS1),P(KS1))
      END DO
*
*       Define auxiliary quantities.
      MASS=0.0
      DO I=1,N
      L=3*(I-1)
      MC(I)=M(INAME(I))
      MASS=MASS+MC(I)
      END DO
      DO I=1,N-1
      TKK(I)=.5D0*(1./MC(I)+1./MC(I+1))
      TK1(I)=-1./MC(I)
      MKK(I)=MC(I)*MC(I+1)
      DO J=I+1,N
      MIJ(I,J)=MC(I)*MC(J)
      MIJ(J,I)=MIJ(I,J)
      END DO
      END DO
*
*       Copy Y-array from COMMON.
      CALL YCOPY(Y)
*
      RETURN
      END
