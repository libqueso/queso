
# type "source setup_env.sh"

# working as of 5/23/12 on hera

use ic-11.1.073
use openmpi-intel-1.4.3

# S/W paths

export PECOS_DIR=/usr/gapps/pecos/libs/intel-11

export GSL_DIR=$PECOS_DIR/gsl-1.15
###export HDF5_DIR=/usr/local/tools/hdf5-intel-serial-1.8.1
export GRVY_DIR=$PECOS_DIR/grvy-0.32.0
export BOOST_DIR=$PECOS_DIR/boost-1.46.1
###export GLPK_DIR=$PECOS_DIR/glpk-4.45
###export TRILINOS_DIR=$PECOS_DIR/trilinos-9.0.3

export QUESO_DIR=$PECOS_DIR/queso-0.45.1

# Library search path

#export LD_LIBRARY_PATH=$GSL_DIR/lib:$HDF5_DIR/lib:$GRVY_DIR/lib:$BOOST_DIR/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$GLPK_DIR/lib:$TRILINOS_DIR/lib:$QUESO_DIR:/lib:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=$GSL_DIR/lib:$GRVY_DIR/lib:$BOOST_DIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$QUESO_DIR:/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/tools/openmpi-intel-1.4.3/lib:$LD_LIBRARY_PATH
