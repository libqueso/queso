#!/bin/bash

# compile fortran likelihood function
mpif90 -c -g likelihood_function.f90


# compile queso main which calls fortran likelihood

mpic++ -c -g -I./ -I$QUESO_DIR/include -I$BOOST_DIR/include \
    -I$HDF5_DIR/include -DMPICH_IGNORE_CXX_SEEK exStatisticalInverseProblem1_gsl.C 

# link final executable

mpif90 -nofor-main -o queso_app likelihood_function.o exStatisticalInverseProblem1_gsl.o \
    -L$QUESO_DIR/lib -Wl,-rpath,$QUESO_DIR/lib -lqueso

#mpic++ -g -I./ \
#    -I$QUESO_DIR/include \
#    -I$BOOST_DIR/include \
#    -I$HDF5_DIR/include  \
#    -DMPICH_IGNORE_CXX_SEEK \
#    exStatisticalInverseProblem1_gsl.C likelihood_function.o -L$QUESO_DIR/lib -Wl,-rpath,$QUESO_DIR/lib -lqueso -lifcore


#    -I$TRILINOS_DIR/include \
#    -I$GLPK_DIR/include  \
