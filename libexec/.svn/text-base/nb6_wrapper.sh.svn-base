#!/bin/bash

function write_nb6wrapper() {
cat > $1 <<EOF 
#!/bin/bash -l

# This is a wrapper script for the execution
# of Nbody6++ jobs in Globus environments.


####################
# VARIABLES
####################

ENABLE_MPI=${ENABLE_MPI}

RESTART=${RESTART}
COMMON_BLOCK_FILE=comm.1_${task_id}

HAVE_INITIAL_DATA_FILE=${HAVE_INITIAL_DATA_FILE}
INITIAL_DATA_FILE=dat.10_${task_id}

NBODY6_SRC_DIR=nbody6-mpi
EXECUTABLE=nbody6
NBODY_SRC=nbody6_src_${task_id}.tar.gz
PARAM_FILE=parameter.in_${task_id}
EXEC_DIR=.nbody6_${task_id}
HOSTS_ENV=hosts.env_${task_id}

ORIGIN_WORKDIR=\$(pwd)
OUTPUT_FILES="comm.1 comm.2 conf.3 lagr.7 hia.12 nbody6.out"

####################
# JOB ENVIRONMENT
####################

NB6_N=$N
NB6_STAGE_OUT_URL=${THIS_URL}/${OUTFILES}
NB6_SUBMISSION_HOST=$(hostname -f)

export NB6_N NB6_STAGE_OUT_URL NB6_SUBMISSION_HOST

###################
# FUNCTIONS
###################

function log()
{
  printf "\$(date) \$1\n"
}


function load_hostenv() {
  if [ -f \$HOSTS_ENV ]; then
    echo ""
    echo "Sourcing \$HOSTS_ENV"
    echo ""
    source \$HOSTS_ENV
  else
    echo "WARNING: Can not source \$HOSTS_ENV. (Not found)"
  fi
}


function init() {
  log "Starting Nbody6++ wrapper..."
  log "Current working directory: \$ORIGIN_WORKDIR"

  load_hostenv
}


function check_program()
{
 PROGRAM_NAME=\$1
 VERSION_ARG=\$2

 PROGRAM=\$(type -p \${PROGRAM_NAME})

 if [ \$? -eq 0 ]; then
   printf "\$PROGRAM_NAME  \t\t\t: \$PROGRAM\n"
   return 0
 else
   printf "\$PROGRAM_NAME  \t\t\t: NOT FOUND\n"
   return 1
 fi
}


function check_exit()
{
  if [ \$? -eq 0 ]; then
    echo "\$* successful"
  else
    echo "\$* failed"
    # Make sure that all stage-out files exist
    cd \${ORIGIN_WORKDIR}
    for i in \$OUTPUT_FILES; do
      touch \${i}_${task_id} 
    done
    exit -1
  fi
}


function print_environment() {

echo '******************************************************************************'
echo '*********************** System Information ***********************************'

printf "\n\
BASIC INFORMATION for \$(uname -n):\n\
------------------------------------------------------------------------------\n\
Machine type           : \$(uname -m) \n\
Operating system       : \$(uname -o) \n\
Kernel-release         : \$(uname -r) \n\
Username               : \$(id -un) \n\
Name of the shell      : \$SHELL \n"

echo;

printf "\
ENVIRONMENT\n\
------------------------------------------------------------------------------\n\
GLOBUS_LOCATION        : \$GLOBUS_LOCATION\n\
GLOBUS_TCP_PORT_RANGE  : \$GLOBUS_TCP_PORT_RANGE\n\
GW_LOCATION            : \$GW_LOCATION\n\
X509_USER_PROXY        : \$X509_USER_PROXY\n\
HOME                   : \$HOME\n\
PWD                    : \$PWD\n"

echo;

printf "\
COMPILERS AND RELATED STUFF\n\
-------------------------------------------------------------------------------\n\
LIBPATH                : \$LIBPATH\n\
LD_LIBRARY_PATH        : \$LD_LIBRARY_PATH\n"


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
check_program gnuplot
check_program globus-url-copy

echo '******************************************************************************'
echo;
}


function print_configuration() {
  echo ""
  echo "Configuration Options:"
  echo "-----------------------"
  
  echo -n -e "MPI\t\t"
  if [ \$ENABLE_MPI = 'yes' ]; then
    echo "enabled"
  else
    echo "disabled"
  fi
  echo -n -e "RESTART\t\t"

  if [ \$RESTART = 'yes' ]; then
    echo "enabled."
  else
    echo "no"
  fi  
  echo "-----------------------"
  echo ""
}


function create_sandbox() {
  echo "*****************************************************************************"
  echo "************************ Create Sandbox *************************************"
  
  # make sure the executable has executable permissions
  # chmod +x \$EXECUTABLE

  echo Creating \$EXEC_DIR
  mkdir \$EXEC_DIR 

  CP_NBODY_SRC="cp \$NBODY_SRC \$EXEC_DIR"
  echo \$CP_NBODY_SRC
  \$CP_NBODY_SRC

  MV_PARAM_FILE="mv \$PARAM_FILE \$EXEC_DIR"
  echo \$MV_PARAM_FILE
  \$MV_PARAM_FILE

  if [ \$RESTART = 'yes' ]; then
    MV_COMMON_BLOCK_FILE="mv \$COMMON_BLOCK_FILE \$EXEC_DIR/comm.1"
    echo \$MV_COMMON_BLOCK_FILE
    \$MV_COMMON_BLOCK_FILE
  fi

  if [ \$HAVE_INITIAL_DATA_FILE = 'yes' ]; then
    MV_INITIAL_DATA_FILE="mv \$INITIAL_DATA_FILE \$EXEC_DIR/dat.10"
    echo \$MV_INITIAL_DATA_FILE
    \$MV_INITIAL_DATA_FILE
  fi
 
  # Copy the plugins to the EXEC_DIR
  echo "Copying the plugins to \${EXEC_DIR}..."
  echo "Working dir: \$(pwd)"
  echo "PLUGINS: nb6_plugin_*_${task_id}"
  for i in nb6_plugin_*_${task_id}; do
    CP_PLUGIN="cp \$i \$EXEC_DIR"
    echo "\$CP_PLUGIN"
    \$CP_PLUGIN
  done

  echo Changing working directory to \$EXEC_DIR
  cd \$EXEC_DIR
  echo Working directory: \$(pwd)

  # Uncompressing the source files
  UNPACK_NBODY_SRC="tar xzvf \$NBODY_SRC"
  echo \$UNPACK_NBODY_SRC
  \$UNPACK_NBODY_SRC 1>/dev/null

  if [ \$? -ne 0 ]; then
    echo "Uncompressing the source files failed"
    exit -1
  fi
  echo ""
}



function compile() {
  echo "*****************************************************************************"
  echo "********************** Compile and Install **********************************"

  # make: Compiling and linking
  cd \$NBODY6_SRC_DIR
  echo Changed working directory to \$(pwd)
  if [ \$ENABLE_MPI = 'yes' ]; then
    CONFIGURE_CMD="./configure --enable-mpi --prefix=\${ORIGIN_WORKDIR}/\${EXEC_DIR} --bindir=\${ORIGIN_WORKDIR}/\${EXEC_DIR}"
  else
    CONFIGURE_CMD="./configure --prefix=\${ORIGIN_WORKDIR}/\${EXEC_DIR} --bindir=\${ORIGIN_WORKDIR}/\${EXEC_DIR}"
  fi

  echo \$CONFIGURE_CMD
  \$CONFIGURE_CMD
  if [ \$? -ne 0 ]; then
    echo "FALLBACK: Default compiler failed. Trying g77..."
    export F77=g77
    \$CONFIGURE_CMD
    check_exit configure
  fi

  MAKE_CMD="make"
  echo \$MAKE_CMD
  \$MAKE_CMD
  check_exit make

  MAKE_INSTALL_CMD="make install"
  echo \$MAKE_INSTALL_CMD
  \$MAKE_INSTALL_CMD
  check_exit make install
  echo ""
}


# Executes the plugins
# Parameters:
 #1  First parameter to the plugin
function plugins_execution()
{
  if [ -f nb6_plugin_*_${task_id} ]; then
    for i in nb6_plugin_*_${task_id}; do
      [ ! -x \$i ] && chmod +x \$i
      EXEC_PLUGIN="./\$i \$1"
      echo "****************************"
      echo "EXECUTING PLUGIN \$i"
      echo "****************************"
      \$EXEC_PLUGIN
    done
  fi
}

function execute() {
  echo "*****************************************************************************"
  echo "******************************* Execute *************************************"
  # These dummy files needs to be created
  # to prevent stage-out failures.

  for i in \$OUTPUT_FILES; do
    touch \${i}_${task_id}
  done

  cd \${ORIGIN_WORKDIR}/\${EXEC_DIR}

  plugins_execution pre_execution

  echo Creating dummy files: \$OUTPUT_FILES
  touch \$OUTPUT_FILES

  log "Executing NBODY6++..."
  NB6_STARTTIME=\$(date +%s)

  if [ x\$ENABLE_MPI = 'xyes' ]; then
    mpirun -np 4 ./\$EXECUTABLE < \$PARAM_FILE &> nbody6.out
  else
    ./\$EXECUTABLE < \$PARAM_FILE 2>&1 | tee nbody6.out
  fi

  NB6_STOPTIME=\$(date +%s)
  exitcode=\$?

  if [ \$exitcode -eq 0 ]; then
   log "NBODY6++ finished successfully"
   echo Execution walltime: \$((\$NB6_STOPTIME-\$NB6_STARTTIME)) seconds
  else
    log "NBODY6++ exited with error code \$?"
  fi
  echo ""

  plugins_execution post_execution
}



function post_execute() {
  echo "*****************************************************************************"
  echo "************************** Post-execute *************************************"
  # GridWay/GridGateWay workaround (full paths
  # are not supported at the moment so we have
  # to move the files from \$EXEC_DIR back to \$ORIGIN_WORKDIR

  echo Moving output files from \$EXEC_DIR to \$ORIGIN_WORKDIR
  for i in \$OUTPUT_FILES; do
    CMD="mv \$i \${ORIGIN_WORKDIR}/\${i}_${task_id}"
    echo \$CMD
    \$CMD
  done

  # Globus cleanup expects PARAM_FILE in ORIGIN_WORKDIR
  mv \$PARAM_FILE \$ORIGIN_WORKDIR
  # Globus cleanup expects EXECUTABLE in ORIGIN_WORKDIR
  mv \$EXECUTABLE \$ORIGIN_WORKDIR

  echo Changing working directory to \$ORIGIN_WORKDIR
  cd \$ORIGIN_WORKDIR
  echo Removing \$EXEC_DIR
  rm -r \$EXEC_DIR

  if [ -d \$EXEC_DIR ]; then
    echo "\$EXEC_DIR could not be removed. Please check this directory."
  fi
  
  log "Wrapper finished."
  exit \$exitcode
}


function deploy() {

  init
  print_environment  
  print_configuration

  create_sandbox
  compile
  execute
  post_execute
}


####################
# MAIN
####################

deploy

EOF
}

