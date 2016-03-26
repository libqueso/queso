dnl ----------------------------------------------------------------
dnl Tests for various C++11 features.  These will probably only work
dnl if they are run after the autoconf test that sets -std=c++11.
dnl ----------------------------------------------------------------


AC_DEFUN([QUESO_TEST_CXX11_ISNAN],
  [
    have_cxx11_isnan=no

    # Only run the test if enablecxx11==yes
    if (test "x$enablecxx11" = "xyes"); then
      AC_MSG_CHECKING(for C++11 std::isnan support)
      AC_LANG_PUSH([C++])

      AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                        [[#include<cmath>]],
                        [[bool isnan = std::isnan(1);]])]
      ,[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_ISNAN, 1, [Flag indicating whether compiler supports std::isnan;])
        have_cxx11_isnan=yes
      ],[
        AC_MSG_RESULT(no)
      ])

      AC_LANG_POP([C++])
    fi

    AM_CONDITIONAL(HAVE_CXX11_ISNAN, test x$have_cxx11_isnan == xyes)
  ])
