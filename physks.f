      SUBROUTINE PHYSKS(XR,PR,Q,P)
*
*       Transformation from physical to KS variables.
*       ---------------------------------------------
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XR(3),PR(3),Q(4),P(4)
*
      R2=XR(1)**2+XR(2)**2+XR(3)**2
      R=SQRT(R2)
      Q(4)=0.
      A=R+ABS(XR(1))
      Q(1)=SQRT(.5D0*A)
      B=1./(2.0D0*Q(1))
      Q(2)=XR(2)*B
      Q(3)=XR(3)*B
      IF(XR(1).lt.0.)THEN
      U1=Q(1)
      Q(1)=Q(2)
      Q(2)=U1
      U3=Q(3)
      Q(3)=Q(4)
      Q(4)=U3
      END IF
      P(1)=2.D0*(+Q(1)*PR(1)+Q(2)*PR(2)+Q(3)*PR(3))
      P(2)=2.D0*(-Q(2)*PR(1)+Q(1)*PR(2)+Q(4)*PR(3))
      P(3)=2.D0*(-Q(3)*PR(1)-Q(4)*PR(2)+Q(1)*PR(3))
      P(4)=2.D0*(+Q(4)*PR(1)-Q(3)*PR(2)+Q(2)*PR(3))
*
      RETURN
      END
