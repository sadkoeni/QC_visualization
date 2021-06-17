//High dEdx analysis description by Sebastian Lehner

All the files required to run the high dEdx analysis on TPC raw data (i.e. HLT clusters) can be found in the directory of the JIRA ticket "ATO-432". If there are technical problems submitting jobs or running the scripts in a job Marian Ivanov should be able to help.

submit.jdl needs to be submitted from aliensh. Also the input files listed in "submit.jdl" need to be accessible on aliensh. One needs to replace "PATH" correctly in "submit.jdl".
The script that will be executed is highdEdxFilter.sh, it consecutively calls the scripts that perform the data preparation, analysis and filtering on the HLT cluster data.

recHLT.C reconstructs the HLT cluster information and writes it to a root file.

processTree() in tpcHLTClusterLoop.C analyses all events (which are referred to via their Global Identification number (GID)). It saves a tree with sector and event wise statistics on the charge amplitude in the individual HLT clusters. In addition it writes a tree (to clusterDump700HLT.root) with the full HLT cluster information for HLT clusters that have a high charge, i.e. Qmax>700.

clusterDump700HLT.root serves as the input for tpcPadRowClustering.C. This script analysis for each event and sector the individual high-charge HLT clusters. It searches for clusters of these HLT clusters, which could point to a stopping nuclear fragment (or monopole). It extracts for each event and sector the number of consecutive pad rows that have an HLT cluster above a certain threshold charge. This information is written to tree in locCount700HLT.root. Based on this tree the file selectedGIDs.root is created which contains the event (i.e. GIDs) which are selected by a criterion that is based on the size of the found clusters of high-charge HLT clusters. These criteria are defined in line 202.

processTree(x,kTRUE) in tpcHLTClusterLoop.C will finally read selectedGIDs.root and write the HLT information for all HLT clusters (i.e. without a charge restriction) in the selected events and sectors to clusterDumpFullHLT.root. This is the last output in of the grid job.

Once clusterDumpFullHLT.root files from the grid are downloaded the XY and YZ projections of the selected events and sectors can be plotted with plotSec.sh. Alternatively the clusterDumpFullHLT.root files could serve as input for further automated event selection.


The information in selectedGIDs.root might be interesting to learn how loser or tighter selection criteria would affect the number of selected events.









