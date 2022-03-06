# 1 "file_init.F"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "file_init.F"
      SUBROUTINE FILE_INIT(ISTART)
*
*
* Opening of Files with proper names.
* ------------------------------------
* ISTART=0: Open all other files after reading input data.
* ISTART=1: Open only unit 1 for mydump.
* ISTART=2: Open only unit 2 for mydump.
*
      INCLUDE 'common6.h'
*
*
      CHARACTER*10 FILE(100)
      CHARACTER*4 FRAG(100)
*
# 25 "file_init.F"
*
      FRAG(1) = 'comm'
      FRAG(2) = 'comm'
      FRAG(3) = 'conf'
      FRAG(4) = 'bdat'
      FRAG(7) = 'lagr'
      FRAG(8) = 'bdat'
      FRAG(9) = 'bdat'
*
      JFMIN = 1
      JFMAX = 9
      IF (ISTART.EQ.1) THEN
          JFMIN = 1
          JFMAX = 1
      END IF
      IF (ISTART.EQ.2) THEN
          JFMIN = 2
          JFMAX = 2
      END IF
* Initialize the file names for read and write.
      DO 100 JF = JFMIN,JFMAX
          WRITE (FILE(JF),112) FRAG(JF),JF
 100 CONTINUE
 112 FORMAT(A4,'.',I1,T4)
      FILE(10) = 'dat.10    '
      FILE(11) = 'esc.11    '
      FILE(12) = 'hia.12    '
      FILE(13) = 'hid.13    '
      FILE(15) = 'per.15    '
      FILE(82) = 'bev.82    '
      FILE(83) = 'sev.83    '
# 138 "file_init.F"
      IF (KZ(1).GT.0.AND.ISTART.EQ.1)
     &OPEN (UNIT=1,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(1))
*
      IF (KZ(2).GT.0.AND.ISTART.EQ.2)
     &OPEN (UNIT=2,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(2))
*
      IF (ISTART.GT.0) RETURN
*
      IF (KZ(3).GT.0)
     &OPEN (UNIT=3,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(3),
     & ACCESS='APPEND')
      IF (KZ(4).GT.0)
     &OPEN (UNIT=4,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(4),
     & ACCESS='APPEND')
      IF (KZ(7).GE.3)
     &OPEN (UNIT=7,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(7),
     & ACCESS='APPEND')
      IF (KZ(8).GT.0.OR.NBIN0.GT.0)
     &OPEN (UNIT=8,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(8),
     & ACCESS='APPEND')
      IF (KZ(8).GE.2.OR.NBIN0.GT.0)
     &OPEN (UNIT=9,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(9),
     & ACCESS='APPEND')
      IF (KZ(22).GT.0)
     &OPEN (UNIT=10,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(10))
      IF (KZ(23).EQ.2.OR.KZ(23).EQ.4)
     &OPEN (UNIT=11,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(11),
     & ACCESS='APPEND')
      IF (KZ(11).EQ.1.OR.KZ(11).EQ.3)
     &OPEN (UNIT=12,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(12),
     & ACCESS='APPEND')
      IF (KZ(8).GT.3)
     &OPEN (UNIT=13,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(13),
     & ACCESS='APPEND')
      IF (BK(4).EQ.1)
     &OPEN (UNIT=15,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(15),
     & ACCESS='APPEND')
      IF (KZ(12).GT.0)
     &OPEN (UNIT=82,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(82),
     & ACCESS='APPEND')
      IF (KZ(12).GT.0)
     &OPEN (UNIT=83,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(83),
     & ACCESS='APPEND')



*
      RETURN
      END
