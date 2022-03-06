#!/bin/sh

# Nbody6++ submission script for AstroGrid-D
# Thomas Bruesemeister <tbruese@ari.uni-heidelberg.de>
#
# $Id: submit.sh 35 2007-10-19 10:47:57Z tbruese $

####################
# DEFAULT VALUES 
####################

## Possible verbosity levels are
# 0 (QUIET)
# 1 (ERROR)
# 2 (WARNING)
# 3 (INFO)
# 4 (DEBUG), not implemented yet
VERBOSITY=3

# Execution host, format: [protocol://]{hostname|hostaddr}[:port][/service]
HOST="hydra.ari.uni-heidelberg.de"  

# Run in batch mode? ('yes', 'no')
BATCH='yes'

# The job-manager which should be used by Globus
FACTORY_TYPE='GW'   # Fork, GW, PBS, SGE, Condor 

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

RESTART='no'
COMMON_BLOCK_FILE='not_specified'

INITIAL_DATA_FILE='not_specified'
HAVE_INITIAL_DATA_FILE='no'

JOB_QUEUE='not_specified'

#echo "THIS_URL is $THIS_URL"


####################
# INCLUDES
####################

. scripts/common.sh
. libexec/nb6_wrapper.sh
. scripts/jobdescription.sh
. scripts/stats/stats.sh


####################
# FUNCTIONS
####################

function print_usage()
{
 
 echo "Submits Nbody6++ jobs to Globus nodes."
 
 echo "Usage: $0 [options] <parameter-file>"
 echo ""
 echo "Options: "
 echo "  -d                  Delegate full credential. [$DELEGATE_CREDENTIAL]"
 echo "  -g host             Submits the job to <host>." [$HOST]
 echo "  -h                  Print this help."
 echo "  -m                  Enable MPI. (Experimental)."
 echo "  -n                  Disable batch mode."
 echo "  -q queue            Use queue <queue>."
 echo "  -s                  Enable job statistics (Experimental)."
 echo "  -t job-manager      Use <job-manager> as Globus Job Manager. [$FACTORY_TYPE]"
 echo ""
 echo "Stage-in Options:"
 echo "  -fr file            Stage-in a common-block file for restart."
 echo "  -fd file            Stage-in a file for initial data of m,r,v."
 echo ""
 echo "Example: ./submit.sh -d var/in1000.comment"
 echo "          (Option -d must be used to provide a proxy for GridWay.)"
 echo ""
 echo "$PACKAGE_NAME, $VERSION"
 echo ""
 RC=$1
 if [ $# -gt 1 ]; then
   shift; echo -e "$@"
 fi
 exit $RC
}


function parse_parameters()
{
  if [ $# -lt 1 ]; then
   print_usage -1 "You have to specify at least a parameter file."
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
      -fd) HAVE_INITIAL_DATA_FILE='yes'; INITIAL_DATA_FILE=$2; shift;;
      -fr) RESTART='yes'; COMMON_BLOCK_FILE=$2; shift;;
       -s) ENABLE_JOB_STATS='yes';;
       -q) JOB_QUEUE=$2; shift;;
       -*) print_usage -1 "Option $1 is not known.";;
        *) PARAMETER_INPUT_FILE=$1;;   # teminate while loop
    esac
    shift
  done

  if [ $PARAMETER_INPUT_FILE = 'not_specified' ]; then
   print_usage -1 "Argument <parameter-file> is missing."
  fi
}


function get_plugin_configuration_hints()
{
  if [ -f ${PLUGINS_DIR}/* ]; then
    for i in "${PLUGINS_DIR}/*"; do
      grep "@CREDENTIAL_DELEGATE" ${i} 1>/dev/null
      if [ $? -eq 0 ];then
        pecho 3 "  Plugin $(basename $i) forces credential delegation."
        DELEGATE_CREDENTIAL='yes'  
      fi
    done
  fi
}


function print_configuration()
{

  echo;
  echo -e "$BOLD_ON Configuration: $BOLD_OFF
------------------------------------------------------"
  echo ""
  get_plugin_configuration_hints
  echo""
  
  echo "  Parameter input file: $PARAMETER_INPUT_FILE"
  if [ x$RESTART = 'xyes' ]; then
    echo "  Common block file   : $COMMON_BLOCK_FILE"
  else
    echo "  Common block file   : n/a"
  fi
  if [ x$HAVE_INITIAL_DATA_FILE = 'xyes' ]; then
    echo "  Initial data file   : $INITIAL_DATA_FILE"
  else
    echo "  Initial data file   : n/a"
  fi
  echo  "  Host                : $HOST
  Job-manager         : $FACTORY_TYPE
  Delegate credential : $DELEGATE_CREDENTIAL
  Job statistics      : $ENABLE_JOB_STATS
  Batch mode          : $BATCH"
  if [ -f $PLUGINS_DIR/* ]; then
    echo -n "  Plugins             : "
    for i in $PLUGINS_DIR/*; do echo -n "$(basename $i) "; done; echo
  fi
}


function init()
{
  get_lock

  task_id=0
  get_task_id task_id

  trap "graceful_exit 0" SIGINT SIGTERM

  OUTFILES="outfiles/$task_id"
  echo "------------------------"
  echo -e "$BOLD_ON Nbody6++ Job #${task_id} $BOLD_OFF"
  echo "------------------------"

  USER_IDENTITY=$(grid-proxy-info -identity)
  
  # CREATE DIRECTORY FOR OUTPUT FILES
  mkdir $OUTFILES

  if [ $? -ne 0 ]; then
    perror 1 "Could not create directory ${OUTPUTFILES}. Exit"
    graceful_exit 1
  fi

  if [ ! -f $PARAMETER_INPUT_FILE ]; then
    perror 1 "ERROR: Parameter input file $PARAMETER_INPUT_FILE not found."
    graceful_exit -1
  fi
}


function package_source_code()
{
  echo;
  echo -e "$BOLD_ON Job Preparation: $BOLD_OFF"
  echo "------------------------------------------------------"
  printf "Packaging the source code..."
  # Package the src files
  if [ -d "nbody6src" ]; then
    if [ -f "nbody6src/src/nbody6.F" ]; then
      cd nbody6src
     if [ ! -f Makefile ]; then
        ./configure 1>/dev/null
        if [ $? -ne 0 ]; then
          perror 1 "Packaging ('configure') failed. exit."
          graceful_exit -1
        fi
      fi
      make dist 1>/dev/null
      if [ $? -ne 0 ]; then
        perror 1 "Packaging ('make dist') failed. exit."
        graceful_exit -1
      fi
      mv nbody6-mpi.tar.gz ../tmp/nbody6_src_${task_id}.tar.gz
      cd ..
      # tar czvf tmp/nbody6_src_${task_id}.tar.gz nbody6src 1>/dev/null
    else
      echo; 
      perror 1 "nbody6.F is missing in src/. exit."
      graceful_exit -1
    fi
  else
    echo;
    echo "The source directory is missing."
    read -p "Should the sources be downloaded from SVN [y/n]? " get_sources
    if [ x$get_sources = "xy" ]; then
      read -p "Which branch should be downloaded [${DEFAULT_BRANCH}]? " nbody_branch
      if [ x$nbody_branch = "x" ]; then
	get_nbody_source  # see common.sh  
      else
        get_nbody_source $nbody_branch   # see common.sh
      fi
      package_source_code
      return
    else
      graceful_exit 0
    fi
    graceful_exit 1
  fi
  echo "done."
}


function generate_job_description()
{
  pecho_n 3 "Generating wrapper..."
  # see libexec/nb6_wrapper.sh
  write_nb6wrapper ${WRAPPER}_${task_id}
  pecho 3 "done."

  pecho_n 3 "Generating job description..."
  # see jobdescription.sh
  write_rsl ${RSL}_${task_id}
  pecho 3 "done."
}


function check_wrapper()
{
  # Make sure that the wrapper has executable permissions 
  if [ -f ${WRAPPER}_${task_id} ]; then
    chmod +x ${WRAPPER}_${task_id}  
  else
    perror 1 Wrapper does not exist. exit.
    graceful_exit -1
  fi
}


function ask_continue()
{
  if [ $VERBOSITY -ge 2 ]; then 
    read -p "Do you really want to continue [y/n]?" cont
    if [ ! x$cont = 'xy' ]; then graceful_exit 1; fi
  fi
}


function check_parameter_input_file()
{
  # Checking if the file exists is done in init()

  ECHO_WARNING="${BOLD_ON}Warning:${BOLD_OFF}"

  read_parameter_input_file $PARAMETER_INPUT_FILE
  echo ""

  # When using a restart (common block) file then
  # check if KSTART >= 2
  if [ x$RESTART = 'xyes' ] && [ $KSTART -lt 2 ];then
    pecho 2 "$ECHO_WARNING You have selected to continue a run by submitting"
    pecho 2 "a common block file but KSTART in the parameter input file is $KSTART"
    ask_continue
  fi

  # When given the instruction to stage a initial data file
  # check if KZ(22) >= 2
  if [ x$HAVE_INITIAL_DATA_FILE = 'xyes' ] && [ $KZ_22 -lt 2 ]; then
    pecho 2 "$ECHO_WARNING You have given the instruction to stage-in a"
    pecho 2 "data input file but KZ(22) in the parameter input file is $KZ_22"
    ask_continue
  fi

}


function delegate_credential()
{
  if [ x$HOST = "xlocalhost" ] || [ x$HOST = "x127.0.0.1" ]; then
    pecho 2 "Using $HOST as host causes problems when trying to"
    pecho 2 "delegate a credential"
    ask_continue
  fi
  
  CREDENTIAL_DELEG_HOST=$HOST
  if [ -f $DELEG_EPR_FILE ]; then
    pecho 2 "Warning: Removing stale $DELEG_EPR_FILE"
  fi
  globus-credential-delegate -h $CREDENTIAL_DELEG_HOST $DELEG_EPR_FILE
  if [ $? -ne 0 ]; then
    perror 1 "Delegating the credentials failed. exit"
    graceful_exit -1
  fi
}


function submit()
{
  echo;
  echo -e "$BOLD_ON Job Submission: $BOLD_OFF"
  echo "------------------------------------------------------"
  if [ $BATCH = 'yes' ]; then
    _SETBATCH="-b"
  fi

  if [ -f $DELEG_EPR_FILE ]; then
    CMD="globusrun-ws -submit $_SETBATCH -F $HOST -Ft $FACTORY_TYPE -S -Jf $DELEG_EPR_FILE \
         -o $OUTFILES/job.epr -Io $OUTFILES/$UUID_FILE -f ${RSL}_${task_id}" 
    pecho 3 $CMD
  else
    CMD="globusrun-ws -submit $_SETBATCH -F $HOST -Ft $FACTORY_TYPE -S -o $OUTFILES/job.epr \
         -Io $OUTFILES/$UUID_FILE -f ${RSL}_${task_id}" 
    pecho 3 $CMD
  fi

  $CMD

  if [ $? -eq 0 ]; then
    job_log "SUCCESS: Submitted job #${task_id}, $PARAMETER_INPUT_FILE to $HOST, job-manager $FACTORY_TYPE"
  else
    job_log "FAILED: Submitting job #${task_id}, $PARAMETER_INPUT_FILE to $HOST, job-manager $FACTORY_TYPE"
    perror 1 "Failed submitting job #${task_id}"
    graceful_exit -1
  fi
}


function send_job_statistics()
{
  # Submit job statistics
  if [ x$ENABLE_JOB_STATS = 'xyes' ]; then
    pecho_n 3 "Sending job statistics..."
    # Generate job statistics RDF
    USER_IDENTITY=$(grid-proxy-info -identity)
    if [ -f $OUTFILES/$UUID_FILE ]; then
      UUID=$(cat $OUTFILES/$UUID_FILE | awk -F: '{ print $2 }')
      generate_rdf
      submit_job_rdf
      if [ $? -eq 0 ]; then
        pecho 3 "done."
      else
        perror 1 "failed."
      fi
    else
      perror 1 "failed (UUID file missing)."
    fi
  fi
}


#####################################################################

####################
# MAIN
####################

# Parse command-line parameters
parse_parameters $@

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

# Get the lock and create job directory
init

print_configuration

check_parameter_input_file

# Job Preparation
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

