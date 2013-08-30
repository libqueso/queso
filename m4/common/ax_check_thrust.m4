##### 
#
# SYNOPSIS
#
# AX_CHECK_THRUST( <Minimum Required Version>, # <package-required=yes/no> )
#
# DESCRIPTION
#
# Provides a --with-thrust=DIR option.  Searches --with-thrust, $THRUST_DIR,
# and the standard directories for Thrust compilers and libraries.
#
# Figures out if Thrust API with THRUST Driver/nvcc is available, i.e. existence of:
#   thrust/*.h
#   libthrustrt.so
#   nvcc
#
# Upon success, locations of these are included in 
#   THRUST_CPPFLAGS and 
#   THRUST_LDFLAGS.
# Path to nvcc is included as
#   NVCC_PATH
# in config.h
# 
# Assumes package is optional unless overridden with $2=yes.
#
# The original author is personally using CUDA such that the .cu code
# is generated at runtime, so don't expect any automake magic to exist
# for compile time compilation of .cu files.
#
# LICENCE
#   Public domain
#
# LAST MODIFICATION
#
#   $Id: ax_check_thrust.m4 $
#
# AUTHOR
#   2010-2012 wili <wili@wili.cc> 
#   2013      Roy Stogner <roystgnr@ices.utexas.edu>
#
##### 

AC_DEFUN([AX_CHECK_THRUST], [

# Provide your THRUST path with this		
AC_ARG_WITH(thrust, [AS_HELP_STRING([--with-thrust[=PREFIX]],[Prefix of your THRUST installation])],
  [with_thrust=$withval
    if test "${with_thrust}" = yes; then
      thrust_prefix=$withval
    elif test "${with_thrust}" = no; then
      thrust_prefix=no
    fi
  ], [])

# Setting the prefix to the default if only --with-thrust was given
if test "$thrust_prefix" != "no"; then
	if test "x${THRUST_DIR}" != "x"; then
		thrust_prefix=${THRUST_DIR}
	elif test "$withval" == "yes"; then
		for dir in /usr /usr/local /usr/local/thrust; do
			if test -x "$dir/bin/nvcc"; then
				echo "Found $dir/bin/nvcc"
				thrust_prefix=$dir
			fi
		done
	fi
fi

is_package_required=ifelse([$2], ,no, $2 )

HAVE_THRUST=1

if test "$thrust_prefix" != "no"; then

THRUST_CXX=$CXX

# Checking for nvcc
AC_MSG_CHECKING([nvcc in $thrust_prefix/bin])
if test -x "$thrust_prefix/bin/nvcc"; then
	AC_MSG_RESULT([found])
	AC_DEFINE_UNQUOTED([NVCC_PATH], ["$thrust_prefix/bin/nvcc"], [Path to nvcc binary])
	THRUST_CXX="$thrust_prefix/bin/nvcc"
else
	AC_MSG_RESULT([not found])
	HAVE_THRUST=0
fi

AC_SUBST([THRUST_CXX])

# We need to add the THRUST search directories for header and lib searches

# Saving the current flags
ax_save_CPPFLAGS="${CPPFLAGS}"
ax_save_LDFLAGS="${LDFLAGS}"

THRUST_CPPFLAGS="-I$thrust_prefix/include"
CPPFLAGS="$THRUST_CPPFLAGS $CPPFLAGS"
THRUST_LDFLAGS="-L$thrust_prefix/lib"
LDFLAGS="$THRUST_LDFLAGS $LDFLAGS"

# Thrust headers require C++
AC_LANG_PUSH([C++])

# And the header and the lib
AC_CHECK_HEADER([thrust/host_vector.h], [], 
	[AC_MSG_NOTICE([Couldn't find thrust/host_vector.h])
	 HAVE_THRUST=0
	 THRUST_CPPFLAGS=""],
	[#include <thrust/host_vector.h>])
AC_CHECK_LIB([cudart], [cudaDriverGetVersion], 
	[THRUST_LIBS='-lcudart'], 
	[AC_MSG_NOTICE([Couldn't find libcudart])
	 HAVE_THRUST=0
	 THRUST_LDFLAGS=""])

AC_LANG_POP([C++])

# Returning to the original flags
CPPFLAGS=${ax_save_CPPFLAGS}
LDFLAGS=${ax_save_LDFLAGS}

if test $HAVE_THRUST = 0 -a x$is_package_required != xno; then
	AC_MSG_ERROR([Thrust not found.  Try --with-thrust to set a prefix])
fi

AC_SUBST([THRUST_CPPFLAGS])
AC_SUBST([THRUST_LDFLAGS])
AC_SUBST([THRUST_LIBS])
AC_SUBST([THRUST_PREFIX])
AC_SUBST([HAVE_THRUST])

else
HAVE_THRUST=0
fi

AM_CONDITIONAL(THRUST_ENABLED,test x$HAVE_THRUST = x1)

])-
