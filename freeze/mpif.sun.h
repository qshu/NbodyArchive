!
! Copyright (c) 1998, by Sun Microsystems, Inc.
! All rights reserved.
!
! User include file for Fortran MPI programs.
!

!
! MPI version
!
       integer*4 MPI_VERSION, MPI_SUBVERSION

       parameter (MPI_VERSION=1)
       parameter (MPI_SUBVERSION=2)
!
! MPI error codes
!
       integer*4 MPI_SUCCESS, MPI_ERR_BUFFER, MPI_ERR_COUNT
       integer*4 MPI_ERR_TYPE, MPI_ERR_TAG, MPI_ERR_COMM
       integer*4 MPI_ERR_RANK, MPI_ERR_ROOT, MPI_ERR_GROUP
       integer*4 MPI_ERR_OP, MPI_ERR_TOPOLOGY, MPI_ERR_DIMS
       integer*4 MPI_ERR_ARG, MPI_ERR_UNKNOWN, MPI_ERR_TRUNCATE
       integer*4 MPI_ERR_OTHER, MPI_ERR_INTERN, MPI_ERR_IN_STATUS
       integer*4 MPI_ERR_PENDING, MPI_ERR_REQUEST, MPI_ERR_FILE
       integer*4 MPI_ERR_NOT_SAME, MPI_ERR_AMODE
       integer*4 MPI_ERR_UNSUPPORTED_DATAREP
       integer*4 MPI_ERR_UNSUPPORTED_OPERATION
       integer*4 MPI_ERR_NO_SUCH_FILE, MPI_ERR_FILE_EXISTS
       integer*4 MPI_ERR_BAD_FILE, MPI_ERR_ACCESS, MPI_ERR_NO_SPACE
       integer*4 MPI_ERR_QUOTA, MPI_ERR_READ_ONLY, MPI_ERR_FILE_IN_USE
       integer*4 MPI_ERR_DUP_DATAREP, MPI_ERR_CONVERSION, MPI_ERR_IO
       integer*4 MPI_ERR_KEYVAL, MPI_ERR_INFO, MPI_ERR_INFO_KEY
       integer*4 MPI_ERR_INFO_NOKEY, MPI_ERR_INFO_VALUE
       integer*4 MPI_ERR_TIMEDOUT, MPI_ERR_RESOURCES
       integer*4 MPI_ERR_TRANSPORT, MPI_ERR_HANDSHAKE
       integer*4 MPI_ERR_SPAWN, MPI_ERR_LASTCODE

       parameter (MPI_SUCCESS=0)
       parameter (MPI_ERR_BUFFER=1)
       parameter (MPI_ERR_COUNT=2)
       parameter (MPI_ERR_TYPE=3)
       parameter (MPI_ERR_TAG=4)
       parameter (MPI_ERR_COMM=5)
       parameter (MPI_ERR_RANK=6)
       parameter (MPI_ERR_ROOT=7)
       parameter (MPI_ERR_GROUP=8)
       parameter (MPI_ERR_OP=9)
       parameter (MPI_ERR_TOPOLOGY=10)
       parameter (MPI_ERR_DIMS=11)
       parameter (MPI_ERR_ARG=12)
       parameter (MPI_ERR_UNKNOWN=13)
       parameter (MPI_ERR_TRUNCATE=14)
       parameter (MPI_ERR_OTHER=15)
       parameter (MPI_ERR_INTERN=16)
       parameter (MPI_ERR_IN_STATUS=17)
       parameter (MPI_ERR_PENDING=18)
       parameter (MPI_ERR_REQUEST=19)
       parameter (MPI_ERR_FILE=20)
       parameter (MPI_ERR_NOT_SAME=21)
       parameter (MPI_ERR_AMODE=22)
       parameter (MPI_ERR_UNSUPPORTED_DATAREP=23)
       parameter (MPI_ERR_UNSUPPORTED_OPERATION=24)
       parameter (MPI_ERR_NO_SUCH_FILE=25)
       parameter (MPI_ERR_FILE_EXISTS=26)
       parameter (MPI_ERR_BAD_FILE=27)
       parameter (MPI_ERR_ACCESS=28)
       parameter (MPI_ERR_NO_SPACE=29)
       parameter (MPI_ERR_QUOTA=30)
       parameter (MPI_ERR_READ_ONLY=31)
       parameter (MPI_ERR_FILE_IN_USE=32)
       parameter (MPI_ERR_DUP_DATAREP=33)
       parameter (MPI_ERR_CONVERSION=34)
       parameter (MPI_ERR_IO=35)
       parameter (MPI_ERR_KEYVAL=36)
       parameter (MPI_ERR_INFO=37)
       parameter (MPI_ERR_INFO_KEY=38)
       parameter (MPI_ERR_INFO_NOKEY=39)
       parameter (MPI_ERR_INFO_VALUE=40)
       parameter (MPI_ERR_TIMEDOUT=41)
       parameter (MPI_ERR_RESOURCES=42)
       parameter (MPI_ERR_TRANSPORT=43)
       parameter (MPI_ERR_HANDSHAKE=44)
       parameter (MPI_ERR_SPAWN=45)
       parameter (MPI_ERR_LASTCODE=46)
!
! miscellaneous MPI constants
!
       integer*4 MPI_ADDRESS_KIND, MPI_OFFSET_KIND
       integer*4 MPI_PROC_NULL
       integer*4 MPI_ANY_SOURCE, MPI_ANY_TAG
       integer*4 MPI_UNDEFINED, MPI_BSEND_OVERHEAD
       integer*4 MPI_MAX_PROCESSOR_NAME, MPI_MAX_ERROR_STRING
       integer*4 MPI_MAX_OBJECT_NAME, MPI_MAX_INFO_KEY
       integer*4 MPI_MAX_INFO_VAL, MPI_MAX_PORT_NAME
       integer*4 MPI_KEYVAL_INVALID
       integer*4 MPI_GRAPH, MPI_CART
       integer*4 MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED
       integer*4 MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE
       integer*4 MPI_MAX_DATAREP_STRING

       parameter (MPI_ADDRESS_KIND=8)
       parameter (MPI_OFFSET_KIND=8)
       parameter (MPI_PROC_NULL=-1)
       parameter (MPI_ANY_SOURCE=-2)
       parameter (MPI_ANY_TAG=-1)
       parameter (MPI_UNDEFINED=-32766)
       parameter (MPI_BSEND_OVERHEAD=512)
       parameter (MPI_MAX_PROCESSOR_NAME=255)
       parameter (MPI_MAX_ERROR_STRING=255)
       parameter (MPI_MAX_OBJECT_NAME=63)
       parameter (MPI_MAX_INFO_KEY=63)
       parameter (MPI_MAX_INFO_VAL=255)
       parameter (MPI_MAX_PORT_NAME=63)
       parameter (MPI_KEYVAL_INVALID=0)
       parameter (MPI_GRAPH=1)
       parameter (MPI_CART=2)
       parameter (MPI_THREAD_SINGLE=1)
       parameter (MPI_THREAD_FUNNELED=2)
       parameter (MPI_THREAD_SERIALIZED=3)
       parameter (MPI_THREAD_MULTIPLE=4)
       parameter (MPI_MAX_DATAREP_STRING=255)
!
!      IO constants
!
       integer*4 MPI_MODE_RDONLY, MPI_MODE_WRONLY, MPI_MODE_RDWR
       integer*4 MPI_MODE_CREATE, MPI_MODE_EXCL
       integer*4 MPI_MODE_DELETE_ON_CLOSE
       integer*4 MPI_MODE_UNIQUE_OPEN, MPI_MODE_SEQUENTIAL
       integer*4 MPI_MODE_APPEND, MPI_SEEK_SET, MPI_SEEK_CUR
       integer*4 MPI_SEEK_END, MPI_DISTRIBUTE_DFLT_DARG
       integer*4 MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_CYCLIC
       integer*4 MPI_DISTRIBUTE_NONE, MPI_ORDER_C, MPI_ORDER_FORTRAN

       parameter (MPI_MODE_RDONLY=0)
       parameter (MPI_MODE_WRONLY=1)
       parameter (MPI_MODE_RDWR=2)
       parameter (MPI_MODE_CREATE=4)
       parameter (MPI_MODE_EXCL=8)
       parameter (MPI_MODE_DELETE_ON_CLOSE=16)
       parameter (MPI_MODE_UNIQUE_OPEN=32)
       parameter (MPI_MODE_SEQUENTIAL=64)
       parameter (MPI_MODE_APPEND=128)

       parameter (MPI_SEEK_SET=0)
       parameter (MPI_SEEK_CUR=1)
       parameter (MPI_SEEK_END=2)

       parameter (MPI_DISTRIBUTE_DFLT_DARG=0)
       parameter (MPI_DISTRIBUTE_BLOCK=1)
       parameter (MPI_DISTRIBUTE_CYCLIC=2)
       parameter (MPI_DISTRIBUTE_NONE=3)

       parameter (MPI_ORDER_C=128)
       parameter (MPI_ORDER_FORTRAN=129)
!
! MPI status
!
       integer*4 MPI_STATUS_SIZE, MPI_SOURCE, MPI_TAG, MPI_ERROR

       parameter (MPI_STATUS_SIZE=6)
       parameter (MPI_SOURCE=2)
       parameter (MPI_TAG=3)
       parameter (MPI_ERROR=4)
!
! null handles
!
      integer*4 MPI_COMM_NULL, MPI_DATATYPE_NULL, MPI_ERRHANDLER_NULL
      integer*4 MPI_GROUP_NULL, MPI_INFO_NULL, MPI_OP_NULL
      integer*4 MPI_REQUEST_NULL

      parameter (MPI_COMM_NULL=0)
      parameter (MPI_DATATYPE_NULL=0)
      parameter (MPI_ERRHANDLER_NULL=0)
      parameter (MPI_GROUP_NULL=0)
      parameter (MPI_INFO_NULL=0)
      parameter (MPI_OP_NULL=0)
      parameter (MPI_REQUEST_NULL=0)
!
! communicator and group comparisons
!
      integer*4 MPI_IDENT, MPI_CONGRUENT, MPI_SIMILAR, MPI_UNEQUAL

      parameter (MPI_IDENT=0)
      parameter (MPI_CONGRUENT=1)
      parameter (MPI_SIMILAR=2)
      parameter (MPI_UNEQUAL=3)
!
! predefined handles
!
       integer*4 MPI_COMM_WORLD, MPI_COMM_SELF
       integer*4 MPI_GROUP_EMPTY
       integer*4 MPI_ERRORS_ARE_FATAL, MPI_ERRORS_RETURN

       parameter (MPI_COMM_WORLD=1)
       parameter (MPI_COMM_SELF=2)
       parameter (MPI_GROUP_EMPTY=3)
       parameter (MPI_ERRORS_ARE_FATAL=4)
       parameter (MPI_ERRORS_RETURN=5)

       integer*4 MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
       integer*4 MPI_LOGICAL, MPI_CHARACTER, MPI_BYTE
       integer*4 MPI_PACKED, MPI_UB, MPI_LB
       integer*4 MPI_INTEGER1, MPI_INTEGER2, MPI_INTEGER4
       integer*4 MPI_INTEGER8, MPI_REAL4, MPI_REAL8, MPI_REAL16
       integer*4 MPI_COMPLEX, MPI_DOUBLE_COMPLEX, MPI_COMPLEX8
       integer*4 MPI_COMPLEX16, MPI_COMPLEX32
       integer*4 MPI_LOGICAL1, MPI_LOGICAL2, MPI_LOGICAL4
       integer*4 MPI_LOGICAL8, MPI_2DOUBLE_PRECISION, MPI_2INTEGER
       integer*4 MPI_2REAL, MPI_2COMPLEX, MPI_2DOUBLE_COMPLEX

       parameter (MPI_INTEGER=10)
       parameter (MPI_REAL8=11)
       parameter (MPI_DOUBLE_PRECISION=12)
       parameter (MPI_COMPLEX=13)
       parameter (MPI_LOGICAL=14)
       parameter (MPI_CHARACTER=15)
       parameter (MPI_BYTE=16)
       parameter (MPI_PACKED=17)
       parameter (MPI_UB=18)
       parameter (MPI_LB=19)
       parameter (MPI_2REAL=20)
       parameter (MPI_2DOUBLE_PRECISION=21)
       parameter (MPI_2INTEGER=22)
       parameter (MPI_DOUBLE_COMPLEX=23)
       parameter (MPI_2COMPLEX=24)
       parameter (MPI_2DOUBLE_COMPLEX=25)
       parameter (MPI_INTEGER2=26)
       parameter (MPI_INTEGER4=27)
       parameter (MPI_REAL4=28)
       parameter (MPI_REAL=29)
       parameter (MPI_INTEGER8=30)
       parameter (MPI_INTEGER1=31)
       parameter (MPI_REAL16=32)
       parameter (MPI_COMPLEX8=33)
       parameter (MPI_COMPLEX16=34)
       parameter (MPI_COMPLEX32=35)
       parameter (MPI_LOGICAL1=36)
       parameter (MPI_LOGICAL2=37)
       parameter (MPI_LOGICAL4=38)
       parameter (MPI_LOGICAL8=39)

       integer*4 MPI_MAX, MPI_MIN, MPI_SUM, MPI_PROD, MPI_LAND
       integer*4 MPI_BAND, MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR
       integer*4 MPI_MAXLOC, MPI_MINLOC

       parameter (MPI_MAX=40)
       parameter (MPI_MIN=41)
       parameter (MPI_SUM=42)
       parameter (MPI_PROD=43)
       parameter (MPI_LAND=44)
       parameter (MPI_BAND=45)
       parameter (MPI_LOR=46)
       parameter (MPI_BOR=47)
       parameter (MPI_LXOR=48)
       parameter (MPI_BXOR=49)
       parameter (MPI_MAXLOC=50)
       parameter (MPI_MINLOC=51)
!
! predefined keys
!
       integer*4 MPI_TAG_UB, MPI_HOST, MPI_IO
       integer*4 MPI_WTIME_IS_GLOBAL, MPI_APPNUM

       parameter (MPI_TAG_UB=54)
       parameter (MPI_HOST=55)
       parameter (MPI_IO=56)
       parameter (MPI_WTIME_IS_GLOBAL=57)
       parameter (MPI_APPNUM=58)
!
! datatype combiners
!
       integer*4 MPI_COMBINER_NAMED, MPI_COMBINER_DUP
       integer*4 MPI_COMBINER_CONTIGUOUS, MPI_COMBINER_VECTOR
       integer*4 MPI_COMBINER_HVECTOR_INTEGER, MPI_COMBINER_HVECTOR
       integer*4 MPI_COMBINER_INDEXED, MPI_COMBINER_HINDEXED_INTEGER
       integer*4 MPI_COMBINER_HINDEXED, MPI_COMBINER_INDEXED_BLOCK
       integer*4 MPI_COMBINER_STRUCT_INTEGER, MPI_COMBINER_STRUCT
       integer*4 MPI_COMBINER_F90_REAL, MPI_COMBINER_F90_COMPLEX
       integer*4 MPI_COMBINER_F90_INTEGER, MPI_COMBINER_RESIZED
       integer*4 MPI_COMBINER_SUBARRAY, MPI_COMBINER_DARRAY

       parameter (MPI_COMBINER_NAMED=0)
       parameter (MPI_COMBINER_DUP=1)
       parameter (MPI_COMBINER_CONTIGUOUS=2)
       parameter (MPI_COMBINER_VECTOR=3)
       parameter (MPI_COMBINER_HVECTOR_INTEGER=4)
       parameter (MPI_COMBINER_HVECTOR=5)
       parameter (MPI_COMBINER_INDEXED=6)
       parameter (MPI_COMBINER_HINDEXED_INTEGER=7)
       parameter (MPI_COMBINER_HINDEXED=8)
       parameter (MPI_COMBINER_INDEXED_BLOCK=9)
       parameter (MPI_COMBINER_STRUCT_INTEGER=10)
       parameter (MPI_COMBINER_STRUCT=11)
       parameter (MPI_COMBINER_RESIZED=12)
       parameter (MPI_COMBINER_SUBARRAY=13)
       parameter (MPI_COMBINER_DARRAY=14)
       parameter (MPI_COMBINER_F90_INTEGER=15)
       parameter (MPI_COMBINER_F90_REAL=16)
       parameter (MPI_COMBINER_F90_COMPLEX=17)
!
! attribute functions
!
       external MPI_NULL_COPY_FN, MPI_NULL_DELETE_FN
       external MPI_COMM_NULL_COPY_FN, MPI_COMM_NULL_DELETE_FN
       external MPI_TYPE_NULL_COPY_FN, MPI_TYPE_NULL_DELETE_FN
       external MPI_DUP_FN, MPI_COMM_DUP_FN
       external MPI_TYPE_DUP_FN
!
! real functions
!
       real*8 MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK
       external MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK
!
! global variables
!
       integer*4 MPI_BOTTOM, MPI_ARGV_NULL, MPI_ARGVS_NULL
       integer*4 MPI_ERRCODES_IGNORE
       integer*4 MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
       integer*4 MPI_STATUSES_IGNORE(MPI_STATUS_SIZE)
       integer*4 MPI_RESERVE(16)

       common /mpipriv/ MPI_BOTTOM
       common /mpipriv/ MPI_ARGV_NULL
       common /mpipriv/ MPI_ARGVS_NULL
       common /mpipriv/ MPI_ERRCODES_IGNORE
       common /mpipriv/ MPI_STATUS_IGNORE
       common /mpipriv/ MPI_STATUSES_IGNORE
       common /mpipriv/ MPI_RESERVE
!
! make sure the common block is static
!
      save /mpipriv/
