#!/usr/bin/env bash
# source $NOTES/JIRA/ATO-432/code/jetTriggerCounter.sh
dataPrefix=/lustre/nyx/alice/DETdata/triggeredESD/
workingDir=$NOTESData/JIRA/ATO-432

production="/alice/data/2015/LHC15o"
pass="/pass1/"

# make directory structure with flitered data
# assumed production and pass varables are defined before
makeDirsPeriod(){
    cd $workingDir/$production/
    find $dataPrefix/$production -iname "Filter*root" | grep $pass > filtered.list
    cat filtered.list |sed s_"$dataPrefix/"__g |sed s_"$pass.*"_"$pass"_| sort -u  > dir.list
    for adir in `cat dir.list`; do
        echo mkdir -p $workingDir$adir
        mkdir -p $workingDir$adir
        cat filtered.list | grep $adir > $workingDir$adir/filtered.list
    done;
}

makeDirsYear(){
  year=$1
  for a in `ls $dataPrefix/alice/data/$year`; do
    echo $a
    production=/alice/data/$year/$a
    echo $production
    makeDirsPeriod
  done;
  # for a in `ls -d LHC*/`; do echo $a; done

}


batchSetting(){
    wdir=`pwd`
    batchCommand="/usr/bin/sbatch"
    batchFlags="--get-user-env --mem-per-cpu=4096 --time=6:00:00"
    batchFlagsM="--get-user-env --mem-per-cpu=4096 --time=2:00:00"
}

# create commands
makeSubmitCommand(){
    nEvents=100000000
    cat <<EOF > commandHisto.sh
#!/bin/bash
    aliroot -n -b -l <<\EOF 2>&1 | tee commandHisto.log
    .L $NOTES/JIRA/ATO-432/code/jetTriggerCounter.C+
     LoadChains(100000);
    MakeHistogramsMult($nEvents);
    MakeHistogramsV0($nEvents);
    .q
EOF
    chmod a+x commandHisto.sh
}

# for a in `ls -d LHC*/`; do echo $a; find `pwd`/$a -iname "filtered.list" |grep 000| xargs dirname > $a/dir.list;  done

submitJobsPeriod(){
   wdir=`pwd`
   for adir in `cat dir.list | grep 000`; do
       cd $adir;
       rm -f $adir/*log
       rm -f $adir/*root
       cp $wdir/commandHisto.sh .
       cp $NOTES/JIRA/ATO-432/code/trigger.list .
       echo mapOutput.root >performance.list
       $batchCommand $batchFlagsM -o commandHisto.blog commandHisto.sh
       cd $wdir
   done;
}

submitJobsYear(){
    wdirY=`pwd`;
    for a in `ls -d LHC*/`; do
        cd $wdirY/$a
        pwd
        makeSubmitCommand
        submitJobsPeriod;
        cd $wdirY
    done;
}
mergeYear(){
    wdirY=`pwd`;
    for a in `ls -d LHC*/`; do
        cd $wdirY/$a
        pwd
        hadd -f  hisOutput.root */pass1/hisOutput*.root
        cd $wdirY
    done;
}