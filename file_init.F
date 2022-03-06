      SUBROUTINE FILE_INIT(ISTART)
*
*
*       Opening of Files with proper names.
*       ------------------------------------
*       ISTART=0: Open all files after reading input data.
*       ISTART=1: Open only unit 1 before read in mydump.
*       ISTART=2: Open all files except unit 1 after read in mydump.
*
      INCLUDE 'common6.h'
*
*
      CHARACTER*10 FILE(100)
      CHARACTER*4 FRAG(100)
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#ifdef ENSEMBLE
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
#endif
*
      FRAG(1) = 'comm'
      FRAG(2) = 'comm'
      FRAG(3) = 'conf'
      FRAG(4) = 'bdat'
      FRAG(7) = 'lagr'
      FRAG(8) = 'bdat'
      FRAG(9) = 'bdat'
*
      JFMAX = 9
      IF (ISTART.EQ.1) JFMAX = 1
*       Initialize the file names for read and write.
      DO 100 JF = 1,JFMAX
          WRITE (FILE(JF),112) FRAG(JF),JF
 100  CONTINUE
 112  FORMAT(A4,'.',I1,T4)
      FILE(10) = 'dat.10    '
      FILE(11) = 'esc.11    '
      FILE(15) = 'per.15    '
#ifdef ENSEMBLE
      DO 101 JF = 1,JFMAX
      if (isize.lt.10) then
          WRITE (FILE(JF),113) FRAG(JF),rank,JF
 113  FORMAT(A4,'000',I1,'.',I1)
      else if (isize.ge.10.and.isize.lt.100) then
          WRITE (FILE(JF),115) FRAG(JF),rank,JF
 115  FORMAT(A4,'00',I2,'.',I1)
      else if (isize.ge.100.and.isize.lt.1000) then
          WRITE (FILE(JF),117) FRAG(JF),rank,JF
 117  FORMAT(A4,'0',I3,'.',I1)
      else if (isize.ge.1000) then
          WRITE (FILE(JF),119) FRAG(JF),rank,JF
 119  FORMAT(A4,I4,'.',I1)
      end if
 101  CONTINUE
*
      if (isize.lt.10) then
          WRITE (FILE(10),213) rank
          WRITE (FILE(11),313) rank
          WRITE (FILE(15),413) rank
 213  FORMAT('dat000',I1,'.10')
 313  FORMAT('esc000',I1,'.11')
 413  FORMAT('per000',I1,'.15')
      else if (isize.ge.100.and.isize.lt.1000) then
          WRITE (FILE(10),215) rank
          WRITE (FILE(11),315) rank
          WRITE (FILE(15),415) rank
 215  FORMAT('dat00',I2,'.10')
 315  FORMAT('esc00',I2,'.11')
 415  FORMAT('per00',I2,'.15')
      else if (isize.ge.100.and.isize.lt.1000) then
          WRITE (FILE(10),217) rank
          WRITE (FILE(11),317) rank
          WRITE (FILE(15),417) rank
 217  FORMAT('dat0',I3,'.10')
 317  FORMAT('esc0',I3,'.11')
 417  FORMAT('per0',I3,'.15')
      else if (isize.ge.1000) then
          WRITE (FILE(10),219) rank
          WRITE (FILE(11),319) rank
          WRITE (FILE(15),419) rank
 219  FORMAT('dat',I4,'.10')
 319  FORMAT('esc',I4,'.11')
 419  FORMAT('per',I4,'.15')
      end if
#endif
#ifdef PARALLEL
      if(rank.eq.0) then
#endif
      IF ((KZ(1).GT.0.AND.ISTART.EQ.0).OR.(ISTART.EQ.1))
     &OPEN (UNIT=1,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(1))
*
      IF (ISTART.EQ.1) RETURN
*
      IF (KZ(2).GT.0)
     &OPEN (UNIT=2,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(2))
      IF (KZ(3).GT.0)
     &OPEN (UNIT=3,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(3))
      IF (BK(3).GT.0)
     &OPEN (UNIT=4,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(4))
      IF (KZ(7).GE.3)
     &OPEN (UNIT=7,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(7))
      IF (KZ(8).GE.2)
     &OPEN (UNIT=9,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(9))
      IF (KZ(22).GT.0)
     &OPEN (UNIT=10,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(10))
      IF (KZ(23).EQ.2.OR.KZ(23).EQ.4)
     &OPEN (UNIT=11,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(11))
      IF (BK(4).EQ.1)
     &OPEN (UNIT=15,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(15))
#ifdef PARALLEL
      end if
#endif
*
      RETURN
      END
