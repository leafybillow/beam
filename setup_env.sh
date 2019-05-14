
source /share/apps/root-5.34.36-build/bin/thisroot.sh

export ANALYZER=/export/home/slactest/AnalysisTool/analyzer
export EVIO_INCDIR=/export/home/slactest/AnalysisTool/analyzer/evio
export EVIO_LIBDIR=/export/home/slactest/AnalysisTool/analyzer/
export SBS_LIBDIR=/export/home/slactest/AnalysisTool/SBS-offline

export PATH=$PATH:$ANALYZER

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SBS_LIBDIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ANALYZER
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$EVIO_LIBDIR

export BEAM_DATA=/export/home/slactest/raw_data/
export BEAM_ROOTFILES=/export/home/slactest/rootfiles/
export BEAM_DB_PATH=/export/home/slactest/DBfiles/
