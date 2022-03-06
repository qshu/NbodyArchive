      SUBROUTINE FILE_INIT(ISTART)
*
*
*       Opening of Files with proper names.
*       ------------------------------------
*       ISTART=0: Open all files after reading input data.
*       ISTART=1: Open only unit 1 before read in mydump.
*       ISTART=2: Open all files except unit 1 after read in mydump.
*
      INCLUDE 'common4.h'
*
*
      CHARACTER*10 FILE(100)
      CHARACTER*4 FRAG(100)
*
      FRAG(1) = 'comm'
      FRAG(2) = 'comm'
      FRAG(3) = 'conf'
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
      FILE(12) = 'hia.12    '
      FILE(13) = 'hid.13    '
      FILE(14) = 'mdot.14   '
      FILE(15) = 'coll.15   '
      FILE(16) = 'coal.16   '
      FILE(17) = 'degen.17  '
      FILE(33) = 'ns.33     '
      FILE(34) = 'bh.34     '
      FILE(98) = 'bss.98    '
*
      IF ((KZ(1).GT.0.AND.ISTART.EQ.0).OR.(ISTART.EQ.1))
     &OPEN (UNIT=1,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(1))
*
      IF (ISTART.EQ.1) RETURN
*
      IF (KZ(2).GT.0)
     &OPEN (UNIT=2,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(2))
      IF (KZ(3).GT.0)
     &OPEN (UNIT=3,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(3))
      IF (KZ(7).GE.3)
     &OPEN (UNIT=7,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(7))
      IF (KZ(8).GE.2)
     &OPEN (UNIT=9,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(9))
      IF (KZ(22).GT.0)
     &OPEN (UNIT=10,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(10))
      IF (KZ(23).EQ.2.OR.KZ(23).EQ.4)
     &OPEN (UNIT=11,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(11))
      IF (KZ(11).GT.0)
     &OPEN (UNIT=12,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(12))
      IF (KZ(8).GE.2)
     &OPEN (UNIT=13,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(13))
      IF (KZ(21).GT.2)
     &OPEN (UNIT=14,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(14))
      IF (KZ(19).GE.3)
     &OPEN (UNIT=15,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(15))
      IF (KZ(19).GE.3)
     &OPEN (UNIT=16,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(16))
      IF (KZ(8).GE.2)
     &OPEN (UNIT=17,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(17))
      IF (KZ(8).GE.2)
     &OPEN (UNIT=33,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(33))
      IF (KZ(8).GE.2)
     &OPEN (UNIT=34,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(34))
      IF (KZ(25).GT.0.AND.KZ(19).GT.0)
     &OPEN (UNIT=98,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(98))
*
      RETURN
      END
