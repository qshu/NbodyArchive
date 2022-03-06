*       kscmn6.
*       -------
*
      LOGICAL  LKSINT
*
      COMMON/KSPAR/  ISTAT(KMAX),IFLAG(KMAX),LKSINT(KMAX),
     &
     &               IPHASEX(KMAX),JCOMPX(KMAX),JCLOSEX(KMAX),
     &               JCMAXX(KMAX),KSPAIRX(KMAX),KS2X(KMAX),
     &               EBCH0X(KMAX),
     &
*       XKS,XDOTKS - for old data
*       X,  XDOT   - for new data
     &               XKS(3,2*KMAX),XDOTKS(3,2*KMAX),
     &
*       T0X,...,LISTX - for old data
*       T0, ...,LIST  - for new data
     &               T0X(2*KMAX),
     &               U0X(4,KMAX),UDOTX(4,KMAX),
     &               FUX(4,KMAX),FUDOTX(4,KMAX),FUDOT2X(4,KMAX),
     &               FUDOT3X(4,KMAX),
     &               HX(KMAX),HDOTX(KMAX),HDOT2X(KMAX),HDOT3X(KMAX),
     &               HDOT4X(KMAX),
     &               DTAUX(KMAX),
     &               TDOT2X(KMAX),TDOT3X(KMAX),
     &               RX(KMAX),GAMMAX(KMAX),
     &               KSLOWX(KMAX),LISTX(LMAX,2*KMAX)
