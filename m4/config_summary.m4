# SYNOPSIS
#
#   Summarizes configuration settings.
#
#   AX_SUMMARIZE_CONFIG([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Outputs a summary of relevant configuration settings.
#
# LAST MODIFICATION
#
#   2009-07-16
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

m4_pattern_allow([BOOST_DIR])

echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version............... : $PACKAGE-$VERSION
echo
echo Debug mode.................... : $enable_debug
echo C++ compiler.................. : $CXX
echo C++ compiler flags............ : $CXXFLAGS
echo ' '
echo GSL_LIBS...................... : $GSL_LIBS
###echo GRVY DIR...................... : $GRVY_PREFIX
echo BOOST_DIR..................... : ${BOOST_DIR}
###echo Boost program options... ..... : $BOOST_PROGRAM_OPTIONS_LDFLAGS $BOOST_PROGRAM_OPTIONS_LIBS
echo ' '
echo Install dir................... : $prefix 
echo Build user.................... : $USER
echo Build host.................... : $BUILD_HOST
echo Configure date................ : $BUILD_DATE
echo Build architecture............ : $BUILD_ARCH
echo Source control revision....... : $BUILD_VERSION
echo
echo Optional Features:

# Optional Features Enabled?

if test "$HAVE_MPI" = "1"; then
  echo '   'Link with MPI.............. : yes
else
  echo '   'Link with MPI.............. : no
fi

if test "$HAVE_GRVY" = "0"; then
  echo '   'Link with GRVY............. : no
else
  echo '   'Link with GRVY............. : yes
fi

if test "$HAVE_GLPK" = "0"; then
  echo '   'Link with GLPK............. : no
else
  echo '   'Link with GLPK............. : yes
fi

if test "$HAVE_HDF5" = "0"; then
  echo '   'Link with HDF5............. : no
else
  echo '   'Link with HDF5............. : yes
fi

if test "$HAVE_TRILINOS" == "0"; then
  echo '   'Link with Trilinos......... : no
else
  echo '   'Link with Trilinos......... : yes
fi

if test "$HAVE_LIBMESH" == "0"; then
  echo '   'Link with libmesh.......... : no
else
  echo '   'Link with libmesh.......... : yes
fi

if test "$HAVE_GCOV_TOOLS" = "0"; then
   echo '   'Enable gcov code coverage.. : no
else
   echo '   'Enable gcov code coverage.. : yes
fi

if test "$HAVE_ANN" == "0"; then
   echo '   'Build internal ANN library. : no
else
   echo '   'Build internal ANN library. : yes
fi

# Paths for optional packages which are enabled

echo 
echo Optional Feature Paths:

if test "$HAVE_GRVY" = "1"; then
   echo '   'GRVY DIR................... : $GRVY_PREFIX
fi

if test "$HAVE_HDF5" = "1"; then
   echo '   'HDF5 DIR................... : $HDF5_PREFIX
fi

if test "$HAVE_GLPK" = "1"; then
   echo '   'GLPK DIR................... : $GLPK_PREFIX
fi

if test "$HAVE_TRILINOS" = "1"; then
   echo '   'Trilinos DIR............... : $TRILINOS_HOME
fi

if test "$HAVE_LIBMESH" == "0"; then
   echo '   'libmesh DIR................ : $LIBMESH_PREFIX
else
   echo '   'libmesh DIR................ : $LIBMESH_PREFIX
fi

echo 		   

echo '-------------------------------------------------------------------------------'

echo
echo Configure complete, now type \'make\' and then \'make install\'.
echo

])
