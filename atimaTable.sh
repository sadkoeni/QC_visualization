#!/bin/bash
# source $NOTES/JIRA/ATO-432/atimaTable.sh
#-
# miranov@lxir128:~$ /WWW/weick/atima/atimawww-1.41 2 4 2 12 26.998 1.632
#  <I> read ATIMA splines for Mg.
# 1: Z=12, A=24.305000, w=1.000000, pot=156.000000
# /WWW/weick/atima/atimawww-1.41 2 4 2 12 26.998 1.632 | grep "dE/dx at entrance" |sed s_.*:__
#
# Energy loss according Table 1.5 http://www.rmki.kfki.hu/~laszloa/downloads/literature/particle_detection_with_drift_chambers.pdf
# 1.52 MeV/g cm−2
# ATIMA: dE/dx at entrance : 0.001515 MeV/(mg/cm^2)

#export  rndm0=32767.
function myrndm2(){
    echo $(echo $RANDOM/$rndm0 | bc -l)
}

runRandom(){
    export  rndm0=32767.
    # Argon properties https://en.wikipedia.org/wiki/Argon
    ZT=18
    AT=39.948
    thickness=1.661 #mg/cm^2     # https://www.engineeringtoolbox.com/gas-density-d_158.html NTP - Normal Temperature and Pressure - is defined as 20oC (293.15 K, 68oF)
    echo "index/D:etype/D:Z/D:A/D:E/D:EN/D:ZT/D:AT/D:Thickness/D:dEdx/D:specRange/D" >elossAr.txt
    npoints=$1
    # generate random fragment
    for ((index=0;index<npoints; index++)); do
        type=$(echo 2*$RANDOM/$rndm0  | bc)
        if  [ $type -gt 0 ]; then
            E=$(echo 3000*$RANDOM/$rndm0 | bc -l)
            #echo "type>0" $E
        else
            E=$(echo 500*"($RANDOM/$rndm0)"^3 | bc -l)    # power law distribution
            #echo "type==0" $E
        fi
        #E=$(echo $(echo 2^$(echo 20*$RANDOM/$rndm0)/1000. | bc -l))
        Z=$(echo 18*$RANDOM/$rndm0 +1 | bc)
        A=$(echo $(echo 2*$Z+  2*$RANDOM/$rndm0-1  | bc -l)/1| bc)
        #E=$(echo 1000*$RANDOM/$rndm0 | bc -l)
        EN=$(echo $E/$A | bc -l)
        dEdx=$(  /WWW/weick/atima/atimawww-1.41 $Z $A  $EN  $ZT $AT $thickness | grep "dE/dx at entrance" |sed s_.*:__ | sed s_MeV.*__) # MeV/(mg/cm^2)\
        range=$( /WWW/weick/atima/atimawww-1.41 $Z $A $EN  $ZT $AT $thickness | grep "range  " | sed s_.*:__ | sed s_mg.*__)
        echo  $index $type $Z $A $E $EN  $ZT $AT $thickness $dEdx $range
    done |tee  -a elossAr.txt
}

runParallel(){
  wdir=$(pwd)
  for i in `seq 40 80`; do
    echo $i ;
    mkdir $wdir/job$i;
    cd  $wdir/job$i;
    runRandom 600 &
    cd $wdir
  done
}

mergeScan(){
    head -n 1 job1/elossAr.txt >elossAr.txt
    cat job*/elossAr.txt | grep -v index | grep -v inf  >>elossAr.txt
}
