# source $NOTES/JIRA/ATO-432/code/monoAnalysis.sh
# make list

alias helpCat='pygmentize -O style=borland,linenos=1 -l bash'

init(){
    # shellcheck disable=SC2154
    export baseDirData="${NOTESData}"/JIRA/ATO-432
    [[ -z "${ALILOG_HOST}" ]] && source "${ALICE_ROOT}"/libexec/alilog4bash.sh
    [[ -z "${PWGPP_runMap}" ]] && source "${ALICE_ROOT}"/libexec/utilities.sh
    source "${ALICE_PHYSICS}"/PWGPP/CalibMacros/AliOCDBtoolkit.sh
    source "${NOTES}"/JIRA/ATO-492/bin/analyzeFiltered.sh

    helpCat<<HELP_USAGE
    List of functions:

    makeList
    makeChunkDirs
    makeRunParallel
HELP_USAGE

}



makeList(){
    heloCat <<HELP_USAGE | cat
    makeList
    Parameters:
       $1 - inputFile
       $2 - number of chunks
    Algorithm:
HELP_USAGE
    #/lustre/alice/users/thellman/NOTESData/alice-tpc-notes/JIRA/ATO-432
    #ls -rS /lustre/alice/users/thellman/NOTESData/alice-tpc-notes/JIRA/ATO-432/*/*/merged_clusterDumpFullHLT.root | tail -n 10 > cluster.list
    find /lustre/alice/DETdata/triggeredRawHighdEdx/alice/cern.ch/user/ -size +40k -iname "clusterDumpFullH*" > clusterAll.list
    echo Number of files $(wc clusterAll.list)
}

makeLists(){
  head -n 1000  clusterAll.list | xargs ls -S | head -n 50 > cluster50.list
}

makeChunkDirs(){
   helpCat <<HELP_USAGE | cat
    makeChunkDir - split input file $1 to $2 chunks
    Parameters:
       $1 - inputFile           - e.g cluster50.list
       $2 - number of chunks    - recommended 3-5
    Algorithm:
HELP_USAGE
  [[ $# -ne 2 ]] &&return
  wdir=$(pwd)
  nChunks=$2
  inputFile=$1
  split --suffix-length=3 -d -l ${nChunks} ${inputFile} Clusters
  for sdir in `ls -d Clusters*`; do
      mkdir -p  ${wdir}/dir${sdir}
      mv $sdir  ${wdir}/dir${sdir}/cluster.list
  done
}

makeRunParallel(){
   helpCat <<HELP_USAGE | cat
    makeRunParallel  $1 jobs
    Parameters:
       $1 - number of jobs    - to run
    Algorithm:
HELP_USAGE
  [[ $# -ne 1 ]] && alilog_error  "makeRunParallel: input parameters not specified" && return
  nJobs=$1
  alilog_info  "makeRunParallel: BEGIN $1"
  ls -d dirClusters0* | head -n ${nJobs}| parallel "echo {}; cd {}; pwd;  root.exe -b -q $NOTES/JIRA/ATO-432/code/monoAnalysis.C+g\(200\) |tee monoAnalysis.log"
  alilog_info  "makeRunParallel: END $1"
}

init



