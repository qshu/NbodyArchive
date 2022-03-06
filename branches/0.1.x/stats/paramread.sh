#!/bin/sh

# Script to generate and submit Nbody6++ metadata to
# Stellaris
# Thomas Bruesemeister <tbruese@ari.uni-heidelberg.de>
#
# $Id: paramread.sh 2956 2007-10-04 11:44:54Z tbruese $

####################
# Expected Parameter input file format:
###################

# KSTART, TCOMP, TCRITp, isernb,iserreg                         (nbody6.F)
# N, NFIX, NCRIT, NRAND, NNBOPT, NRUN                           (input.F)
# ETAI, ETAR, RS0, DTADJ, DELTAT, TCRIT, QE, RBAR, ZMBAR        (input.F)
# (KZ(J),J=1,40)                                                (input.F)
# (BK(J),J=1,10)                                                (input.F)
# DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX                         (input.F)
# ALPHA, BODY1, BODYN, NBIN0, ZMET, EPOCH0, DTPLOT              (data.F)
# Q, VXROT, VZROT, RSPH2                                        (scale.F)
# NBIN, SEMI0, ECC0, RATIO,RANGE, NSKIP, IDORM                  (binpop.f)


####################
# VARIABLES
####################

READOUT_LINES=12
CURRENT_LINE=1

LINE_TOKENS[1]=5
LINE_TOKENS[2]=6
LINE_TOKENS[3]=9
LINE_TOKENS[4]=10
LINE_TOKENS[5]=10
LINE_TOKENS[6]=10
LINE_TOKENS[7]=10
LINE_TOKENS[8]=10
LINE_TOKENS[9]=6
LINE_TOKENS[10]=7
LINE_TOKENS[11]=4
LINE_TOKENS[12]=7

####################
# FUNCTIONS
####################

function read_parameter_input_file()
{
  exec 3<$1

  for i in {1..12}; do
    read <&3 line[$i]
  done

get 1 1; KSTART=$val
get 1 2; TCOMP=$val
get 1 3; TCRITp=$val
get 1 4; isernb=$val
get 1 5; iserreg=$val
get 2 1; N=$val
get 2 2; NFIX=$val
get 2 3; NCRIT=$val
get 2 4; NRAND=$val
get 2 5; NNBOPT=$val
get 2 6; NRUN=$val
get 3 1; ETAI=$val
get 3 2; ETAR=$val
get 3 3; RS0=$val
get 3 4; DTADJ=$val
get 3 5; DELTAT=$val
get 3 6; TCRIT=$val
get 3 7; QE=$val
get 3 8; RBAR=$val
get 3 9; ZMBAR=$val
get 4 1; KZ_1=$val
get 4 2; KZ_2=$val
get 4 3; KZ_3=$val
get 4 4; KZ_4=$val
get 4 5; KZ_5=$val
get 4 6; KZ_6=$val
get 4 7; KZ_7=$val
get 4 8; KZ_8=$val
get 4 9; KZ_9=$val
get 4 10; KZ_10=$val
get 5 1; KZ_11=$val
get 5 2; KZ_12=$val
get 5 3; KZ_13=$val
get 5 4; KZ_14=$val
get 5 5; KZ_15=$val
get 5 6; KZ_16=$val
get 5 7; KZ_17=$val
get 5 8; KZ_18=$val
get 5 9; KZ_19=$val
get 5 10;KZ_20=$val
get 6 1; KZ_21=$val
get 6 2; KZ_22=$val
get 6 3; KZ_23=$val
get 6 4; KZ_24=$val
get 6 5; KZ_25=$val
get 6 6; KZ_26=$val
get 6 7; KZ_27=$val
get 6 8; KZ_28=$val
get 6 9; KZ_29=$val
get 6 10; KZ_30=$val
get 7 1; KZ_31=$val
get 7 2; KZ_32=$val
get 7 3; KZ_33=$val
get 7 4; KZ_34=$val
get 7 5; KZ_35=$val
get 7 6; KZ_36=$val
get 7 7; KZ_37=$val
get 7 8; KZ_38=$val
get 7 9; KZ_39=$val
get 7 10; KZ_40=$val
get 8 1; BK_1=$val
get 8 2; BK_2=$val
get 8 3; BK_3=$val
get 8 4; BK_4=$val
get 8 5; BK_5=$val
get 8 6; BK_6=$val
get 8 7; BK_7=$val
get 8 8; BK_8=$val
get 8 9; BK_9=$val
get 8 10; BK_10=$val
get 9 1; DTMIN=$val
get 9 2; RMIN=$val
get 9 3; ETAU=$val
get 9 4; ECLOSE=$val
get 9 5; GMIN=$val
get 9 6; GMAX=$val
get 10 1; ALPHA=$val
get 10 2; BODY1=$val
get 10 3; BODYN=$val
get 10 4; NBIN0=$val
get 10 5; ZMET=$val
get 10 6; EPOCH0=$val
get 10 7; DTPLOT=$val
get 11 1; Q=$val
get 11 2; VXROT=$val
get 11 3; VZROT=$val
get 11 4; RSPH2=$val
get 12 1; NBIN=$val
get 12 2; SEMI0=$val
get 12 3; ECC0=$val
get 12 4; RATIO=$val
get 12 5; RANGE=$val
get 12 6; NSKIP=$val
get 12 7; IDORM=$val


}


function get()
{
  val=$(echo ${line[$1]} | awk "{printf \$$2}")
}

# echo "Number of Particles: $N"

#for i in {1..12}; do
#  echo "line $i: ${line[$i]}"
#done


#cat in1000.comment | while read line; do
#
#  if [ $CURRENT_LINE -le $READOUT_LINES ]; then
#    if [ $CURRENT_LINE -eq 1 ]; then
#      #KSTART=$(echo $line | awk '{printf $1}')
#      KSTART="bla"
#      TESTVAR=$KSTART; echo testvarset
#    fi
#    #echo KSTARTA:$KSTART
#  fi
#  CURRENT_LINE=$((CURRENT_LINE+1))
#  echo "KSTARTB:$KSTART, $$"
#done

# cat in1000.comment | while read line; do
#  
#  if [ $CURRENT_LINE -le $READOUT_LINES ]; then
#    echo "TOKENS $TOKENS"
#    while [ $i -le ${LINE_TOKENS[$CURRENT_LINE]} ]; do
#      printf "Zeile: $CURRENT_LINE, Spalte $i:"
#      echo $line | awk "{printf \$$i}"
#      echo ""
#      i=$((i+1))
#    done
#
#  fi
#  CURRENT_LINE=$((CURRENT_LINE+1))
# done


