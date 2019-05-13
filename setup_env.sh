
source /usr/local/root-5/bin/thisroot.sh
export ANALYZER=/home/yetao/workarea/slac_beam/analyzer
export EVIO_INCDIR=/home/yetao/workarea/slac_beam/analyzer/evio
export EVIO_LIBDIR=/home/yetao/workarea/slac_beam/analyzer
export SBS_LIBDIR=/home/yetao/workarea/slac_beam/SBS-offline/build
export PATH=$PATH:$ANALYZER

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SBS_LIBDIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ANALYZER
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$EVIO_LIBDIR

export BEAM_DATA=/home/yetao/workarea/slac_beam/raw_data/
export BEAM_ROOTFILES=/home/yetao/workarea/slac_beam/rootfiles/
export BEAM_DB_PATH=/home/yetao/workarea/slac_beam/beam/DBfiles/
