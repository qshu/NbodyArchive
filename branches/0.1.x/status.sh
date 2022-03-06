#!/bin/sh

# Shows the status of Nbody6++ jobs.
# Thomas Bruesemeister <tbruese@ari.uni-heidelberg.de>
#
# $Id: status.sh 2955 2007-10-04 11:32:27Z tbruese $

. common.sh

function print_usage()
{
  echo "Shows the status of Nbody6++ jobs."
  echo ""
  echo "Usage: $0 { jobid | -a }"
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

if [ $1 = "-a" ]; then
  cd outfiles
  printf "Job \t Status\n"
  echo "---------------"
  for i in $(ls -1 | sort -n); do
    if [ -d $i ]; then
      if [ -f $i/job.epr ]; then 
        printf "%3.d\t" $i
	STATUS=$(globusrun-ws -status -j $i/job.epr 2>/dev/null)
	if [ $? -eq 0 ]; then
           PRINT_STATUS=$(echo $STATUS | awk -F: '{ print $2 }')
	   printf "$PRINT_STATUS\n"
        else
           printf "ERROR \n"
	fi
      fi
    fi
  done
  echo ""
else
  if [ -f outfiles/$1/job.epr ]; then
    globusrun-ws -status -j outfiles/$1/job.epr
    exit $?
  else
    echo "outfiles/$1/job.epr does not exist."
  fi
fi

