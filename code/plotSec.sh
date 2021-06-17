for i in $(find -maxdepth 1 -mindepth 1 -type d);
 do echo $i;
#if [ ! -f ./$i/merged* ]; then 
cp plotSec.C $i

 cd $i

 hadd -f merged_clusterDumpFullHLT.root clust*
 hadd -f selectedGIDs.root selectedGIDs_*.root
 root -n -b -l -q "plotSec.C()"
 cd ..
#fi
done
