for i in {0001..1783}; 
do
if [ ! -f ./$i/clusterDumpFullHLT_0.root ]; then

echo get files nr $i

alien_cp alien:///alice/cern.ch/user/s/selehner/triggeredRaw18q/alice/data/2018/LHC18q/000296619/$i/clusterDumpFullHLT_0.root file:./$i;
alien_cp alien:///alice/cern.ch/user/s/selehner/triggeredRaw18q/alice/data/2018/LHC18q/000296619/$i/clusterDumpFullHLT_1.root file:./$i;
alien_cp alien:///alice/cern.ch/user/s/selehner/triggeredRaw18q/alice/data/2018/LHC18q/000296619/$i/clusterDumpFullHLT_2.root file:./$i;
alien_cp alien:///alice/cern.ch/user/s/selehner/triggeredRaw18q/alice/data/2018/LHC18q/000296619/$i/clusterDumpFullHLT_3.root file:./$i;

alien_cp alien:///alice/cern.ch/user/s/selehner/triggeredRaw18q/alice/data/2018/LHC18q/000296619/$i/selectedGIDs_0.root file:./$i; 
alien_cp alien:///alice/cern.ch/user/s/selehner/triggeredRaw18q/alice/data/2018/LHC18q/000296619/$i/selectedGIDs_1.root file:./$i;
alien_cp alien:///alice/cern.ch/user/s/selehner/triggeredRaw18q/alice/data/2018/LHC18q/000296619/$i/selectedGIDs_2.root file:./$i;
alien_cp alien:///alice/cern.ch/user/s/selehner/triggeredRaw18q/alice/data/2018/LHC18q/000296619/$i/selectedGIDs_3.root file:./$i;

else 
echo files nr $i already existing

fi
done
