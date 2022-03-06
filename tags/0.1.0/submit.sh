#!/bin/sh

# NBODY6++ submission script for AstroGrid-D
# Thomas Bruesemeister <tbruese@ari.uni-heidelberg.de>
#
# $Id: submit.sh 2955 2007-10-04 11:32:27Z tbruese $

####################
# DEFAULT VALUES 
####################

# Execution host, format: [protocol://]{hostname|hostaddr}[:port][/service]
HOST=hydra  

# Run in batch mode? ('yes', 'no')
BATCH='yes'

# The job-manager which should be used by Globus
FACTORY_TYPE=GW   # Fork, GW, PBS, SGE, Condor 

# Enable MPI? ('yes', 'no') 
ENABLE_MPI='no'

# Enable job statistics? ('yes', 'no')
ENABLE_JOB_STATS='no'

# Delegate a full credential?
DELEGATE_CREDENTIAL='no'

########################################################################

# This URL must point to *this* location
THIS_URL=gsiftp://$(hostname -f)$(pwd)

PARAMETER_INPUT_FILE='not_specified'

#echo "THIS_URL is $THIS_URL"


####################
# INCLUDES
####################

. common.sh
. stats/stats.sh


####################
# FUNCTIONS
####################

function print_usage()
{
 
 echo "Submits Nbody6++ jobs to Globus nodes."
 
 echo "Usage: $0 [options] <parameter-file>"
 echo ""
 echo "Options: "
 echo "  -d                  Delegate full credential."
 echo "  -g host             Submits the job to <host>." [$HOST]
 echo "  -h                  Print this help."
 echo "  -n                  Disable batch mode."
 echo "  -m                  Enable MPI. (Experimental)."
 echo "  -s                  Enable job statistics (Experimental)."
 echo "  -t job-manager      Use <job-manager> as Globus Job Manager. [$FACTORY_TYPE]"
 echo ""
 echo "$PACKAGE_NAME, $VERSION"
 echo ""
 exit $1
}


function parse_parameters()
{
  if [ $# -lt 1 ]; then
   print_usage -1
  fi

  while [ $# -gt 0 ]
  do
    case "$1" in 
      -d) DELEGATE_CREDENTIAL='yes';;
      -g) HOST="$2"; shift;;
      -h) print_usage 0;;
      -t) FACTORY_TYPE=$2; shift;;
      -n) BATCH='no';;
      -m) ENABLE_MPI='yes';;
      -s) ENABLE_JOB_STATS='yes';;
      -*) print_usage -1;;
       *) PARAMETER_INPUT_FILE=$1;;   # teminate while loop
    esac
    shift
  done

  if [ $PARAMETER_INPUT_FILE = 'not_specified' ]; then
   echo "Argument <parameter-file> is missing."
   echo ""
   print_usage -1
  fi
}


function print_parameters()
{
  echo "DELEGATE_CREDENTIAL: $DELEGATE_CREDENTIAL"
  echo "HOST: $HOST"
  echo "FACTORY_TYPE: $FACTORY_TYPE"
  echo "PARAMTERE_INPUT_FIL: $PARAMETER_INPUT_FILE"
  echo "ENABLE_JOB_STATS: $ENABLE_JOB_STATS"
  echo "BATCH: $BATCH"
}


function init()
{
  get_lock

  task_id=0
  get_task_id task_id

  trap graceful_exit SIGINT SIGTERM

  OUTFILES="outfiles/$task_id"
  echo "------------------------"
  echo "NBODY6++ Job #${task_id}"
  echo "------------------------"

  USER_IDENTITY=$(grid-proxy-info -identity)
  
  # CREATE DIRECTORY FOR OUTPUT FILES
  mkdir $OUTFILES

  if [ $? -ne 0 ]; then
    echo "Could not directory ${OUTPUTFILES}. Exit"
    graceful_exit 1
  fi

  if [ -f $PARAMETER_INPUT_FILE ]; then
    echo "Using Parameter input file: $PARAMETER_INPUT_FILE"
  else
    echo "ERROR: Parameter input file $PARAMETER_INPUT_FILE not found."
    graceful_exit -1
  fi
}


function package_source_code()
{
  printf "Packaging the source code..."
  # Package the src files
  if [ -d "nbody6src" ]; then
    if [ -f "nbody6src/src/nbody6.F" ]; then
      cd nbody6src
     if [ ! -f Makefile ]; then
        ./configure 1>/dev/null
        if [ $? -ne 0 ]; then
          echo "Packaging ('configure') failed. exit."
          graceful_exit -1
        fi
      fi
      make dist 1>/dev/null
      if [ $? -ne 0 ]; then
        echo "Packaging ('make dist') failed. exit."
        graceful_exit -1
      fi
      mv nbody6-mpi.tar.gz ../tmp/nbody6_src_${task_id}.tar.gz
      cd ..
      # tar czvf tmp/nbody6_src_${task_id}.tar.gz nbody6src 1>/dev/null
    else
      echo; 
      echo "nbody6.F is missing in src/. exit."
      graceful_exit -1
    fi
  else
    echo;
    echo "nbody6src directory is missing. exit."
    graceful_exit -1
  fi
  echo "done."
}


function generate_job_description()
{
  # PREPROCESS WRAPPER_TEMPLATE and RSL_TEMPLATE
  sed -e "s/#TASK_ID/$task_id/g" $WRAPPER_TEMPLATE > ${WRAPPER_TEMPLATE}.tmp
  sed -e "s/#ENABLE_MPI/$ENABLE_MPI/g" ${WRAPPER_TEMPLATE}.tmp > ${WRAPPER}_${task_id}
  sed -e "s/#TASK_ID/$task_id/g" ${RSL_TEMPLATE} > ${RSL}.tmp_${task_id}
  sed -e "s|#PARAMETER_INPUT_FILE|$PARAMETER_INPUT_FILE|g" ${RSL}.tmp_${task_id} > ${RSL}.tmp2_${task_id} 
  sed -e "s|#THIS_URL|$THIS_URL|g" ${RSL}.tmp2_${task_id} > ${RSL}_${task_id}
}


function check_wrapper()
{
  # Make sure that the wrapper has executable permissions 
  if [ -f ${WRAPPER}_${task_id} ]; then
    chmod +x ${WRAPPER}_${task_id}  
  else
    echo Wrapper does not exist. exit.
    graceful_exit -1
  fi
}


function delegate_credential()
{
  CREDENTIAL_DELEG_HOST=$HOST
  if [ -f $DELEG_EPR_FILE ]; then
    echo "Warning: Removing stale $DELEG_EPR_FILE"
  fi
  globus-credential-delegate -h $CREDENTIAL_DELEG_HOST $DELEG_EPR_FILE
  if [ $? -ne 0 ]; then
    echo "Delegating the credentials failed. exit"
    graceful_exit -1
  fi
}


function submit()
{
  if [ $BATCH = 'yes' ]; then
    echo "setbatch"
    _SETBATCH="-b"
  fi

  if [ -f $DELEG_EPR_FILE ]; then
    CMD="globusrun-ws -submit $_SETBATCH -F $HOST -Ft $FACTORY_TYPE -S -Jf $DELEG_EPR_FILE \
         -o $OUTFILES/job.epr -Io $OUTFILES/$UUID_FILE -f ${RSL}_${task_id}" 
    echo $CMD
  else
    CMD="globusrun-ws -submit $_SETBATCH -F $HOST -Ft $FACTORY_TYPE -S -o $OUTFILES/job.epr \
         -Io $OUTFILES/$UUID_FILE -f ${RSL}_${task_id}" 
    echo $CMD
  fi

  $CMD
}


function send_job_statistics()
{
  if [ $? -eq 0 ]; then
  # Submit job statistics
  if [ x$ENABLE_JOB_STATS = 'xyes' ]; then
    # Generate job statistics RDF
    echo "Job statistics enabled"
    USER_IDENTITY=$(grid-proxy-info -identity)
    if [ -f $OUTFILES/$UUID_FILE ]; then
      UUID=$(cat $OUTFILES/$UUID_FILE | awk -F: '{ print $2 }')
      read_parameter_input_file $PARAMETER_INPUT_FILE
      generate_rdf
      submit_job_rdf
    fi
  fi
  job_log "SUCCESS: Submitted job #${task_id}, $PARAMETER_INPUT_FILE to $HOST, job-manager $FACTORY_TYPE" 
  else
    job_log "FAILED: Submitting job #${task_id}, $PARAMETER_INPUT_FILE to $HOST, job-manager $FACTORY_TYPE"
    echo "STOP."
    graceful_exit -1
  fi
}


#####################################################################

####################
# MAIN
####################

# Checking Globus Environment, see common.sh
check_globus
# Checking the proxy, see common.sh
check_proxy

# Checking for required programs
check_program "globusrun-ws"
check_program "grid-proxy-info"
check_program "globus-credential-delegate"
check_program "sed"
check_program "awk"

# Parse command-line parameters
parse_parameters $@

# Get the lock and create job directory
init

package_source_code

generate_job_description

# Check if wrapper has been generated
# and has executable permissions
check_wrapper

if [ x$DELEGATE_CREDENTIAL = 'xyes' ]; then
  delegate_credential
fi

submit
send_job_statistics

graceful_exit 0

