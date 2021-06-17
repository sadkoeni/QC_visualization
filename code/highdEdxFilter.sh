date > log.log
cat  wn.xml  | sed s_.*turl=\"__ | sed s_\.root.*_\.root_ | grep alien\.*root | head -n 4   > chunk.list

a=0

for f in $(cat chunk.list);
    do  
        var=${f#*LHC18q/000};
        run=${var%/[0-9]*};
        mkdir -p $run;
        cd $run
        mkdir ./$a
        cd $a
        cp ../../../recHLT.C ./;
        cp ../../../tpcHLTClusterLoop.C ./;
        cp ../../../highdEdxFilter.C ./;
        cp ../../../tpcPadRowClustering.C ./;

        root -n -b -l -q "highdEdxFilter.C(\"$f\",$a)" &> root_$a.log

        cp clusterDumpFullHLT.root ../../../clusterDumpFullHLT_$a.root;
        cp selectedGIDs.root ../../../selectedGIDs_$a.root;
        cp eventDumpHLT.root ../../../eventDumpHLT_$a.root;
        cp *.log ../../../;      
        cd ../../../;

        a=$[a+1];

    done

#hadd merged_clusterDumpFullHLT.root */*/*/clusterDumpFullHLT.root       //if this uses too much disk space just move split clusterDump files to base dir and merge them later locally
#hadd merged_selectedGIDs.root */*/*/selectedGIDs_$a.root;

#hadd merged_eventDumpHLT_.root */*/*/eventDumpHLT_*.root


rm */*/*/clusterDump700HLT.root


pwd > log.log
find ./ -exec du -sh {} \ >> log.log;
