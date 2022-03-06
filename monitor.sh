#!/bin/sh

# Shows the status of Nbody6++ jobs.
# Thomas Bruesemeister <tbruese@ari.uni-heidelberg.de>
#
# $Id: monitor.sh 35 2007-10-19 10:47:57Z tbruese $

. scripts/common.sh

function print_usage()
{
  echo "Attach to a running job."
  echo ""
  echo "Usage: $0 jobid"
  echo ""
  echo "$PACKAGE_NAME, $VERSION"
  echo ""
  exit 1
}


if [ $# -lt 1 ]; then
  print_usage 1
fi

check_globus
check_proxy


echo "ATTENTION: Please *DO NOT* use CTRL-c to 
      stop your monitoring session otherwise
      you kill your job. Use CTRL-z to 
      detach from that session. This is 
      bad behaviour from globusrun-ws.
      Monitoring jobs through the GW job
      manager does not work at the moment."

read -p "Continue [y/n]?" cont

if [ x$cont = 'xy' ]; then

  if [ -f outfiles/$1/job.epr ]; then
    globusrun-ws -monitor -s -j outfiles/$1/job.epr # 2>&1 | tee monitor.out
    exit $?
  else
    echo "outfiles/$1/job.epr does not exist."
  fi

fi


