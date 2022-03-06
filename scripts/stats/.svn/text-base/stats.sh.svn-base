#!/bin/sh 

# Script to generate and submit Nbody6++ metadata to
# Stellaris
# Thomas Bruesemeister <tbruese@ari.uni-heidelberg.de>
#
# $Id: stats.sh 3048 2007-10-11 17:52:30Z engro $

####################
# VARIABLES
####################
RDF_SERVER="mintaka.aip.de"
RDF_PORT="24060" 
URL="http://${RDF_SERVER}:${RDF_PORT}/context/nbody6/jobs?action=add"

RDF_DOCUMENT=tmp/jobstat.rdf

PUT_CMD="curl --connect-timeout 10 -d @${RDF_DOCUMENT} ${URL}"
PUT_CMD2="wget -q --post-file=${RDF_DOCUMENT} ${URL}"


####################
# INCLUDES
####################

. scripts/stats/paramread.sh

####################
# FUNCTIONS
####################

function url_encode()
{
  result=$(echo "$@" | od -t x1 -A n | tr " " %)
}

function put_rdf()
{
  echo "$@" >> $RDF_DOCUMENT
}


# KSTART, TCOMP, TCRITp, isernb,iserreg                         (nbody6.F)
# N, NFIX, NCRIT, NRAND, NNBOPT, NRUN                           (input.F)
# ETAI, ETAR, RS0, DTADJ, DELTAT, TCRIT, QE, RBAR, ZMBAR        (input.F)
# (KZ(J),J=1,40)                                                (input.F)
# (BK(J),J=1,10)                                                (input.F)
# DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX                         (input.F)
# ALPHA, BODY1, BODYN, NBIN0, ZMET, EPOCH0, DTPLOT              (data.F)
# Q, VXROT, VZROT, RSPH2                                        (scale.F)
# NBIN, SEMI0, ECC0, RATIO,RANGE, NSKIP, IDORM                  (binpop.f)


function generate_rdf()
{
  if [ -f $RDF_DOCUMENT ]; then
    rm $RDF_DOCUMENT
  fi

  put_rdf '<?xml version="1.0"?>'
  put_rdf '<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"'
  put_rdf '         xmlns:nb6="http://www.gac-grid.org/nb6/jobs#">'
  put_rdf "  <rdf:Description rdf:about=\"http://www.gac-grid.org/nb6/jobs/$UUID\">"
  put_rdf "    <nb6:userdn>$USER_IDENTITY</nb6:userdn>"
  put_rdf "    <nb6:jobname>n/a</nb6:jobname>"
  put_rdf "    <nb6:description>n/a</nb6:description>"
  put_rdf "    <nb6:localJobId>$task_id</nb6:localJobId>"
  put_rdf "    <nb6:starttime>$(date)</nb6:starttime>"
  put_rdf "    <nb6:stoptime></nb6:stoptime>"
  put_rdf "    <nb6:parameterfile>$PARAMETER_INPUT_FILE</nb6:parameterfile>"
  put_rdf "    <nb6:host>$HOST</nb6:host>"
  put_rdf "    <nb6:job-manager>$FACTORY_TYPE</nb6:job-manager>"
  
  put_rdf "    <nb6:KSTART>$KSTART</nb6:KSTART>"
  put_rdf "    <nb6:TCOMP>$TCOMP</nb6:TCOMP>"
  put_rdf "    <nb6:TCRITp>$TCRITp</nb6:TCRITp>"
  put_rdf "    <nb6:isernb>$isernb</nb6:isernb>"
  put_rdf "    <nb6:iserreg>$iserreg</nb6:iserreg>"
  put_rdf "    <nb6:N>$N</nb6:N>"
  put_rdf "    <nb6:NFIX>$NFIX</nb6:NFIX>"
  put_rdf "    <nb6:NCRIT>$NCRIT</nb6:NCRIT>"
  put_rdf "    <nb6:NRAND>$NRAND</nb6:NRAND>"
  put_rdf "    <nb6:NNBOPT>$NNBOPT</nb6:NNBOPT>"
  put_rdf "    <nb6:NRUN>$NRUN</nb6:NRUN>"
  put_rdf "    <nb6:ETAI>$ETAI</nb6:ETAI>"
  put_rdf "    <nb6:ETAR>$ETAR</nb6:ETAR>"
  put_rdf "    <nb6:RS0>$RS0</nb6:RS0>"
  put_rdf "    <nb6:DTADJ>$DTADJ</nb6:DTADJ>"
  put_rdf "    <nb6:DELTAT>$DELTAT</nb6:DELTAT>"
  put_rdf "    <nb6:TCRIT>$TCRIT</nb6:TCRIT>"
  put_rdf "    <nb6:QE>$QE</nb6:QE>"
  put_rdf "    <nb6:RBAR>$RBAR</nb6:RBAR>"
  put_rdf "    <nb6:ZMBAR>$ZMBAR</nb6:ZMBAR>"
 
  put_rdf "    <nb6:KZ1>$KZ_1</nb6:KZ1>"
  put_rdf "    <nb6:KZ2>$KZ_2</nb6:KZ2>"
  put_rdf "    <nb6:KZ3>$KZ_3</nb6:KZ3>"
  put_rdf "    <nb6:KZ4>$KZ_4</nb6:KZ4>"
  put_rdf "    <nb6:KZ5>$KZ_5</nb6:KZ5>"
  put_rdf "    <nb6:KZ6>$KZ_6</nb6:KZ6>"
  put_rdf "    <nb6:KZ7>$KZ_7</nb6:KZ7>"
  put_rdf "    <nb6:KZ8>$KZ_8</nb6:KZ8>"
  put_rdf "    <nb6:KZ9>$KZ_9</nb6:KZ9>"
  put_rdf "    <nb6:KZ10>$KZ_10</nb6:KZ10>"
  
  put_rdf "    <nb6:KZ11>$KZ_11</nb6:KZ11>"
  put_rdf "    <nb6:KZ12>$KZ_12</nb6:KZ12>"
  put_rdf "    <nb6:KZ13>$KZ_13</nb6:KZ13>"
  put_rdf "    <nb6:KZ14>$KZ_14</nb6:KZ14>"
  put_rdf "    <nb6:KZ15>$KZ_15</nb6:KZ15>"
  put_rdf "    <nb6:KZ16>$KZ_16</nb6:KZ16>"
  put_rdf "    <nb6:KZ17>$KZ_17</nb6:KZ17>"
  put_rdf "    <nb6:KZ18>$KZ_18</nb6:KZ18>"
  put_rdf "    <nb6:KZ19>$KZ_19</nb6:KZ19>"
  put_rdf "    <nb6:KZ20>$KZ_20</nb6:KZ20>"
  
  put_rdf "    <nb6:KZ21>$KZ_21</nb6:KZ21>"
  put_rdf "    <nb6:KZ22>$KZ_22</nb6:KZ22>"
  put_rdf "    <nb6:KZ23>$KZ_23</nb6:KZ23>"
  put_rdf "    <nb6:KZ24>$KZ_24</nb6:KZ24>"
  put_rdf "    <nb6:KZ25>$KZ_25</nb6:KZ25>"
  put_rdf "    <nb6:KZ26>$KZ_26</nb6:KZ26>"
  put_rdf "    <nb6:KZ27>$KZ_27</nb6:KZ27>"
  put_rdf "    <nb6:KZ28>$KZ_28</nb6:KZ28>"
  put_rdf "    <nb6:KZ29>$KZ_29</nb6:KZ29>"
  put_rdf "    <nb6:KZ30>$KZ_30</nb6:KZ30>"
  
  put_rdf "    <nb6:KZ31>$KZ_31</nb6:KZ31>"
  put_rdf "    <nb6:KZ32>$KZ_32</nb6:KZ32>"
  put_rdf "    <nb6:KZ33>$KZ_33</nb6:KZ33>"
  put_rdf "    <nb6:KZ34>$KZ_34</nb6:KZ34>"
  put_rdf "    <nb6:KZ35>$KZ_35</nb6:KZ35>"
  put_rdf "    <nb6:KZ36>$KZ_36</nb6:KZ36>"
  put_rdf "    <nb6:KZ37>$KZ_37</nb6:KZ37>"
  put_rdf "    <nb6:KZ38>$KZ_38</nb6:KZ38>"
  put_rdf "    <nb6:KZ39>$KZ_39</nb6:KZ39>"
  put_rdf "    <nb6:KZ40>$KZ_40</nb6:KZ40>"

  put_rdf "    <nb6:BK1>$BK_1</nb6:BK1>"
  put_rdf "    <nb6:BK2>$BK_2</nb6:BK2>"
  put_rdf "    <nb6:BK3>$BK_3</nb6:BK3>"
  put_rdf "    <nb6:BK4>$BK_4</nb6:BK4>"
  put_rdf "    <nb6:BK5>$BK_5</nb6:BK5>"
  put_rdf "    <nb6:BK6>$BK_6</nb6:BK6>"
  put_rdf "    <nb6:BK7>$BK_7</nb6:BK7>"
  put_rdf "    <nb6:BK8>$BK_8</nb6:BK8>"
  put_rdf "    <nb6:BK9>$BK_9</nb6:BK9>"
  put_rdf "    <nb6:BK10>$BK_10</nb6:BK10>"

  put_rdf "    <nb6:DTMIN>$DTMIN</nb6:DTMIN>"
  put_rdf "    <nb6:RMIN>$RMIN</nb6:RMIN>"
  put_rdf "    <nb6:ETAU>$ETAU</nb6:ETAU>"
  put_rdf "    <nb6:ECLOSE>$ECLOSE</nb6:ECLOSE>"
  put_rdf "    <nb6:GMIN>$GMIN</nb6:GMIN>"
  put_rdf "    <nb6:GMAX>$GMAX</nb6:GMAX>"
  
  put_rdf "    <nb6:ALPHA>$ALPHA</nb6:ALPHA>"
  put_rdf "    <nb6:BODY1>$BODY1</nb6:BODY1>"
  put_rdf "    <nb6:BODYN>$BODYN</nb6:BODYN>"
  put_rdf "    <nb6:NBIN0>$NBIN0</nb6:NBIN0>"
  put_rdf "    <nb6:ZMET>$ZMET</nb6:ZMET>"
  put_rdf "    <nb6:EPOCH0>$EPOCH0</nb6:EPOCH0>"
  put_rdf "    <nb6:DTPLOT>$DTPLOT</nb6:DTPLOT>"
  
  put_rdf "    <nb6:Q>$Q</nb6:Q>"
  put_rdf "    <nb6:VXROT>$VXROT</nb6:VXROT>"
  put_rdf "    <nb6:VZROT>$VZROT</nb6:VZROT>"
  put_rdf "    <nb6:RSPH2>$RSPH2</nb6:RSPH2>"
  
  put_rdf "    <nb6:NBIN>$NBIN</nb6:NBIN>"
  put_rdf "    <nb6:SEMI0>$SEMI0</nb6:SEMI0>"
  put_rdf "    <nb6:ECC0>$ECC0</nb6:ECC0>"
  put_rdf "    <nb6:RATIO>$RATIO</nb6:RATIO>"
  put_rdf "    <nb6:RANGE>$RANGE</nb6:RANGE>"
  put_rdf "    <nb6:NSKIP>$NSKIP</nb6:NSKIP>"
  put_rdf "    <nb6:IDORM>$IDORM</nb6:IDORM>"

  put_rdf "  </rdf:Description>"
  put_rdf "</rdf:RDF>"
}


function submit_job_rdf()
{
  
  if [ -f $RDF_DOCUMENT ]; then
    if type curl 1>/dev/null; then
      $PUT_CMD
    elif type wget 1>/dev/null; then
      $PUT_CMD2
    else
      echo "WARNING: Cannot submit job statistics. No suitable program found."
      return 1
    fi  
    return $?
  else
    return 1
  fi  
}
