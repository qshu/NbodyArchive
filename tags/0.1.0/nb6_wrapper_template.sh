#!/bin/bash -l

## ATTENTION: This file will be preprocessed
## by sed which means that for example that task_id
## is replaced with the real task id 

# This is a wrapper script for the execution
# of NBODY6(++) Jobs in Globus environments


#if [ -f /etc/profile ]; then
#  echo "Sourcing /etc/profile"
#  source /etc/profile
#fi

#if [ -f /opt/d-grid/grid-env.sh ]; then
#  echo "Souring grid-env"
#  source /opt/d-grid/grid-env.sh 
#fi

function log()
{
  printf "$(date) $1\n"
}


log "Starting NBODY6++ wrapper..."
ORIGIN_WORKDIR=$(pwd)
log "Working directory: $ORIGIN_WORKDIR"


ENABLE_MPI=#ENABLE_MPI
NBODY6_SRC_DIR=nbody6-mpi
EXECUTABLE=nbody6
NBODY_SRC=nbody6_src_#TASK_ID.tar.gz
PARAM_FILE=parameter.in_#TASK_ID
EXEC_DIR=.nbody6_#TASK_ID

OUTPUT_FILES="comm.1 comm.2 conf.3 lagr.7 hia.12 nbody6.out"


if [ x$ENABLE_MPI = 'xyes' ]; then
  echo "Load module MPI"
  module load mpi
fi    


function check_program()
{
 PROGRAM_NAME=$1
 VERSION_ARG=$2

 PROGRAM=$(type -p ${PROGRAM_NAME})

 if [ $? -eq 0 ]; then
   printf "$PROGRAM_NAME  \t\t\t: $PROGRAM\n"
 else
   printf "$PROGRAM_NAME  \t\t\t: NOT FOUND\n"
 fi

}


function check_exit()
{
  if [ $? -eq 0 ]; then
    echo "$* successful"
  else
    echo "$* failed"
    # Make sure that all stage-out files exist
    cd ${ORIGIN_WORKDIR}
    for i in $OUTPUT_FILES; do
      touch ${i}_#TASK_ID 
    done
    exit -1
  fi
}


echo '******************************************************************************'
echo '*********************** System Information ***********************************'

printf "\n\
BASIC INFORMATION for $(uname -n):\n\
------------------------------------------------------------------------------\n\
Machine type           : $(uname -m) \n\
Operating system       : $(uname -o) \n\
Kernel-release         : $(uname -r) \n\
Username               : $(id -un) \n\
Name of the shell      : $SHELL \n"

echo;

printf "\
ENVIRONMENT\n\
------------------------------------------------------------------------------\n\
GLOBUS_LOCATION        : $GLOBUS_LOCATION\n\
GLOBUS_TCP_PORT_RANGE  : $GLOBUS_TCP_PORT_RANGE\n\
GW_LOCATION            : $GW_LOCATION\n\
X509_USER_PROXY        : $X509_USER_PROXY\n\
HOME                   : $HOME\n\
PWD                    : $PWD\n"

echo;

printf "\
COMPILERS AND RELATED STUFF\n\
-------------------------------------------------------------------------------\n\
LIBPATH                : $LIBPATH\n\
LD_LIBRARY_PATH        : $LD_LIBRARY_PATH\n"


check_program gcc --version
check_program g77 --version
check_program gfortran --version
check_program ifort
check_program pgf77
check_program pgf90
check_program mpirun
check_program mpif77 
check_program make
check_program cpp


echo '******************************************************************************'
echo '******************************************************************************'

echo;
echo;

echo "-----------------------"
echo "Configuration Options:"
if [ $ENABLE_MPI = 'yes' ]; then
  echo "MPI is enabled."
else
  echo "MPI is disabled."
fi
echo "-----------------------"


# make sure the executable has executable permissions
# chmod +x $EXECUTABLE

echo Creating $EXEC_DIR
mkdir $EXEC_DIR 

CP_NBODY_SRC="cp $NBODY_SRC $EXEC_DIR"
echo $CP_NBODY_SRC
$CP_NBODY_SRC


#MV_EXECUTABLE="mv $EXECUTABLE $EXEC_DIR"
#echo $MV_EXECUTABLE
#$MV_EXECUTABLE

MV_PARAM_FILE="mv $PARAM_FILE $EXEC_DIR"
echo $MV_PARAM_FILE
$MV_PARAM_FILE


echo Changing working directory to $EXEC_DIR
cd $EXEC_DIR
echo Working directory: $(pwd)

# Uncompressing the source files
UNPACK_NBODY_SRC="tar xzvf $NBODY_SRC"
echo $UNPACK_NBODY_SRC
$UNPACK_NBODY_SRC 1>/dev/null

if [ $? -ne 0 ]; then
  echo "Uncompressing the source files failed"
  exit -1
fi

# make: Compiling and linking
cd $NBODY6_SRC_DIR
echo Changed working directory to $(pwd)
if [ $ENABLE_MPI = 'yes' ]; then
  CONFIGURE_CMD="./configure --enable-mpi --prefix=${ORIGIN_WORKDIR}/${EXEC_DIR} --bindir=${ORIGIN_WORKDIR}/${EXEC_DIR}"
else
  CONFIGURE_CMD="./configure --prefix=${ORIGIN_WORKDIR}/${EXEC_DIR} --bindir=${ORIGIN_WORKDIR}/${EXEC_DIR}"
fi

echo $CONFIGURE_CMD
$CONFIGURE_CMD
if [ $? -ne 0 ]; then
  echo "FALLBACK: Default compiler failed. Trying g77..."
  export F77=g77
  $CONFIGURE_CMD
  check_exit configure
fi

MAKE_CMD="make"
echo $MAKE_CMD
$MAKE_CMD
check_exit make

MAKE_INSTALL_CMD="make install"
echo $MAKE_INSTALL_CMD
$MAKE_INSTALL_CMD
check_exit make install


cd ${ORIGIN_WORKDIR}/${EXEC_DIR}

echo Creating dummy files: $OUTPUT_FILES
touch $OUTPUT_FILES

log "Executing NBODY6++..."
NB6_STARTTIME=$(date +%s)

if [ x$ENABLE_MPI = 'xyes' ]; then
  mpirun -np 4 ./$EXECUTABLE < $PARAM_FILE &> nbody6.out
else
  ./$EXECUTABLE < $PARAM_FILE &> nbody6.out
fi

NB6_STOPTIME=$(date +%s)
exitcode=$?

if [ $exitcode -eq 0 ]; then
 log "NBODY6++ finished successfully"
 echo Execution walltime: $(($NB6_STOPTIME-$NB6_STARTTIME)) seconds
else
 log "NBODY6++ exited with error code $?"
fi

# GridWay/GridGateWay workaround (full paths
# are not supported at the moment so we have
# to move the files from $EXEC_DIR back to $ORIGIN_WORKDIR

echo Moving output files from $EXEC_DIR to $ORIGIN_WORKDIR
for i in $OUTPUT_FILES; do
  CMD="mv $i ${ORIGIN_WORKDIR}/${i}_#TASK_ID"
  echo $CMD
  $CMD
done

# Globus cleanup expects PARAM_FILE in ORIGIN_WORKDIR
mv $PARAM_FILE $ORIGIN_WORKDIR
# Globus cleanup expects EXECUTABLE in ORIGIN_WORKDIR
mv $EXECUTABLE $ORIGIN_WORKDIR

echo Changing working directory to $ORIGIN_WORKDIR
cd $ORIGIN_WORKDIR
echo Removing $EXEC_DIR
rm -r $EXEC_DIR

if [ -d $EXEC_DIR ]; then
echo "$EXEC_DIR could not be removed. Please check this directory."
fi


log "Wrapper finished."

exit $exitcode

