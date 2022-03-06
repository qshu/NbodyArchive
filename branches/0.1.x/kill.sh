#!/bin/sh


# Kills Nbody6++ jobs.
# Thomas Bruesemeister <tbruese@ari.uni-heidelberg.de>
#
# $Id: kill.sh 2955 2007-10-04 11:32:27Z tbruese $


####################
# INCLUDES
####################

. common.sh


####################
# FUNCTIONS
####################

function print_usage()
{
  echo "Kills Nbody6++ jobs."
  echo ""
  echo "Usage: $0 jobid1 jobid2 ..jobidn"
  echo ""
  echo "$PACKAGE_NAME, $VERSION"
  echo ""
  exit $1
}


####################
# MAIN
####################

if [ $# -lt 1 ]; then
  print_usage 1
fi

check_globus
check_proxy


for i in $@; do

if [ -f outfiles/$i/job.epr ]; then
  globusrun-ws -kill -j outfiles/$i/job.epr
  if [ $? -eq 0 ]; then
    job_log "SUCCESS: Job #$i killed."
  else
    job_log "FAILURE: Could not kill #$i."
  fi
else
  echo "outfiles/$i/job.epr does not exist."
fi

if [ -d outfiles/$i ] ; then
  printf "Should the directory outfiles/$i be deleted? (n/y)"
   read kb_input
   if [ x$kb_input = "xy" ]; then
     rm -r outfiles/$i
     if [ $? -eq 0 ]; then
       echo "Directory outfiles/$i has been deleted."
     else
       echo "Error deleting directory outfiles/$i."
     fi
   else
    echo "Leaving directory untouched."
   fi
fi

done

