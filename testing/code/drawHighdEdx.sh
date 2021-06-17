#!/usr/bin/env bash
# source $NOTES/JIRA/ATO-432/code/drawHighdEdx.sh
dataPrefix=/lustre/nyx/alice/DETdata/triggeredESD/
workingDir=$NOTESData/JIRA/ATO-432

production="/alice/data/2015/LHC15o"
pass="/pass1/"

# make directory structure with flitered data
# assumed production and pass varables are defined before
makeDirs(){
    find $dataPrefix/$production -iname "Filter*root" | grep $pass > filtered.list
    cat filtered.list |sed s_"$dataPrefix/"__g |sed s_"$pass.*"_"$pass"_| sort -u  > dir.list
    for adir in `cat dir.list`; do
        echo mkdir -p $workingDir$adir
        mkdir -p $workingDir$adir
        cat filtered.list | grep $adir > $workingDir$adir/filtered.list
    done;
}

batchSetting(){
    wdir=`pwd`
    batchCommand="/usr/bin/sbatch"
    batchFlags="--get-user-env --mem-per-cpu=4096 --time=6:00:00"
    batchFlagsM="--get-user-env --mem-per-cpu=4096 --time=2:00:00"
}

# create commands
makeSubmitCommand(){
    nEvents=10000000
    cat <<EOF > commandHisto.sh
#!/bin/bash
    aliroot -n -b -l <<\EOF 2>&1 | tee commandHisto.log
    .L $NOTES/JIRA/ATO-432/code/drawHighdEdx.C+
     LoadChain("filtered.list",100000);
    MakeHistograms($nEvents)
    MakeMaps();
    LoadSummaryTree();
    MakeNDLocalRegressions();
    .q
EOF
    chmod a+x commandHisto.sh
    cat <<EOF > commandND.sh
#!/bin/bash
    aliroot -n -b -l <<\EOF 2>&1 | tee commandND.log
    .L $NOTES/JIRA/ATO-432/code/drawHighdEdx.C+
    LoadSummaryTree();
    MakeNDLocalRegressions();
    .q
EOF
    chmod a+x commandND.sh
}

# Submit histogramming
submitJobs(){
   wdir=`pwd`
   #find `pwd` -iname "filtered.list" | xargs dirname > dir.list
   for adir in `cat dir.list`; do
       cd $adir;
       rm $adir/*log
       rm $adir/*root
       cp $wdir/commandHisto.sh .
       echo mapOutput.root >performance.list
       $batchCommand $batchFlagsM -o commandHisto.blog commandHisto.sh
       cd $wdir
   done;
}
# Submit histogramming
submitND(){
   wdir=`pwd`
   #find `pwd` -iname "filtered.list" | xargs dirname > dir.list
   for adir in `cat dir.list`; do
       cd $adir;
       cp $wdir/commandND.sh .
       echo mapOutput.root >performance.list
       $batchCommand $batchFlagsM -o commandND.blog commandND.sh
       cd $wdir
   done;
}



mergeHistograms(){
    hadd -f  hisOutput.root $(find -iname  hisOutput.root | grep pass)
}