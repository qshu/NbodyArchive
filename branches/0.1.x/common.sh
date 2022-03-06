###################
# VARIABLES
###################

VERSION='v0.1.1 ($Revision: 2964 $)'
PACKAGE_NAME='Nbody6++ deployment package for AstroGrid-D'

NBODY_REPOSITORY="http://svn.ari.uni-heidelberg.de/repos/nbody"
DEFAULT_BRANCH="trunk"
CHECKOUT_TARGET="nbody6src"

WRAPPER_TEMPLATE=nb6_wrapper_template.sh
WRAPPER=tmp/nb6_wrapper.sh
RSL_TEMPLATE=nbody6_template.rsl
RSL=tmp/nbody6.rsl
DELEG_EPR_FILE=deleg.epr
UUID_FILE=job.uuid
JOB_LOG="var/job.log"


###################
# COMMON FUNCTIONS
###################

function job_log()
{
  printf "$(date) $@\n" >> $JOB_LOG
}


function get_task_id()
{
  if [ ! -f task.id ]; then
    touch task.id
  fi
  
  TASK_ID_FILE=task.id
  read task_id < $TASK_ID_FILE
  task_id=$((task_id + 1))
  echo $task_id > $TASK_ID_FILE
  
  return $task_id
}


function get_lock()
{
 while [ -f .lock ];do
   echo .lock exists. Waiting...[Force exit using CTRL-c]
   sleep 3
 done
 touch .lock
}


function graceful_exit()
{
  if [ -f $DELEG_EPR_FILE ]; then
    rm $DELEG_EPR_FILE
  fi
  
  if [ -f .lock ]; then
    echo "Removing .lock"
    rm .lock
  fi
  
  echo Finalizing done.
  if [ $1 -ne 0 ];then
    echo ""
    echo >&2 "ERROR.STOP."
  fi
  
  exit $1
}


function check_globus()
{
  if [ -z $GLOBUS_LOCATION ]; then
    echo GLOBUS_LOCATION is not set.
    graceful_exit -1
  fi

}


function check_proxy()
{
  grid-proxy-info -exists 2> /dev/null

  if [ "x$?" = "x0" ]; then
    export X509_USER_PROXY=`grid-proxy-info -path`
  else
    echo No valid proxy certificate found.
    graceful_exit -1
  fi
}

function check_program()
{
  type $1 &>/dev/null

  if [ $? -ne 0 ]; then
    echo "Error: Program $1 is not installed or is not in PATH"
    graceful_exit -1
  fi  
}

function get_nbody_source()
{
  check_program svn
  
  echo -n "Getting source code from SVN server..."
  if [ x${1} = "x" ]; then
    CMD="svn export ${NBODY_REPOSITORY}/${DEFAULT_BRANCH} ${CHECKOUT_TARGET}"
  else
    CMD="svn export ${NBODY_REPOSITORY}/${1} ${CHECKOUT_TARGET}"
  fi

  job_log "Downloading source code from ${NBODY_REPOSITORY}"
  job_log "$CMD"

  $CMD

  if [ $? -eq 0 ]; then
    echo "done."
  else
    echo "failed."
    exit 1
  fi
}

