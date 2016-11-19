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

AC_DEFUN([QUESO_TEST_CXX11_ISFINITE],
  [
    have_cxx11_isfinite=no

    # Only run the test if enablecxx11==yes
    if (test "x$enablecxx11" = "xyes"); then
      AC_MSG_CHECKING(for C++11 std::isfinite support)
      AC_LANG_PUSH([C++])

      AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                        [[#include<cmath>]],
                        [[bool isnan = std::isfinite(1);]])]
      ,[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_ISFINITE, 1, [Flag indicating whether compiler supports std::isfinite;])
        have_cxx11_isfinite=yes
      ],[
        AC_MSG_RESULT(no)
      ])

      AC_LANG_POP([C++])
    fi

    AM_CONDITIONAL(HAVE_CXX11_ISFINITE, test x$have_cxx11_isfinite == xyes)
  ])

AC_DEFUN([QUESO_TEST_CXX11_UNIQUE_PTR],
  [
    have_cxx11_unique_ptr=no

    # Only run the test if enablecxx11==yes
    if (test "x$enablecxx11" = "xyes"); then
      AC_MSG_CHECKING(for C++11 std::unique_ptr support)
      AC_LANG_PUSH([C++])


      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <iostream>
      @%:@include <memory>
      struct Foo
      {
        Foo()      { std::cout << "Foo::Foo\n";  }
        ~Foo()     { std::cout << "Foo::~Foo\n"; }
      };
          ]], [[
      {
        // up now owns a Foo
        std::unique_ptr<Foo> up(new Foo);
      } // Foo deleted when up goes out of scope
      ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_UNIQUE_PTR, 1, [Flag indicating whether compiler supports std::unique_ptr])
        have_cxx11_unique_ptr=yes
      ],[
        AC_MSG_RESULT(no)
      ])

      AC_LANG_POP([C++])
    fi

    AM_CONDITIONAL(HAVE_CXX11_UNIQUE_PTR, test x$have_cxx11_unique_ptr == xyes)
  ])

AC_DEFUN([QUESO_TEST_CXX11_SHARED_PTR],
  [
    have_cxx11_shared_ptr=no

    # Only run the test if enablecxx11==yes
    if (test "x$enablecxx11" = "xyes"); then
      AC_MSG_CHECKING(for C++11 std::shared_ptr support)
      AC_LANG_PUSH([C++])


      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <iostream>
      @%:@include <memory>
      struct Foo
      {
        Foo()      { std::cout << "Foo::Foo\n";  }
        ~Foo()     { std::cout << "Foo::~Foo\n"; }
      };
          ]], [[
      {
        // up now owns a Foo
        std::shared_ptr<Foo> up(new Foo);
      } // Foo deleted when up goes out of scope
      ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_SHARED_PTR, 1, [Flag indicating whether compiler supports std::shared_ptr])
        have_cxx11_shared_ptr=yes
      ],[
        AC_MSG_RESULT(no)
      ])

      AC_LANG_POP([C++])
    fi

    AM_CONDITIONAL(HAVE_CXX11_SHARED_PTR, test x$have_cxx11_shared_ptr == xyes)
  ])
