if GIT_CLONE
BUILT_SOURCES = .license.stamp
endif

check_PROGRAMS = # Append below
check_PROGRAMS += get_set_row_column
check_PROGRAMS += get_min_max_vec
check_PROGRAMS += multiple_rhs_matrix_solve
check_PROGRAMS += power_method
check_PROGRAMS += inverse_power_method

AM_CPPFLAGS =
AM_CPPFLAGS += -I.
AM_CPPFLAGS += -I$(top_builddir)/inc
AM_CPPFLAGS +=  $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(ANN_CFLAGS)

LIBS         = -L$(top_builddir)/src/ -lqueso

if GRVY_ENABLED	
  AM_CPPFLAGS += $(GRVY_CFLAGS)
endif

if TRILINOS_ENABLED
  AM_CPPFLAGS += -I$(TRILINOS_INCLUDE) 
  LIBS += -lteuchoscore -lteuchoscomm -lteuchosnumerics -lteuchosparameterlist -lteuchosremainder -lepetra
endif

if GLPK_ENABLED	
  AM_CPPFLAGS += $(GLPK_CFLAGS)
  LIBS += $(GLPK_LIBS)
endif

if HDF5_ENABLED	
  AM_CPPFLAGS += $(HDF5_CFLAGS) 
  LIBS += $(HDF5_LIBS)
endif

# C Tests

get_set_row_column_SOURCES = get_set_row_column.C get_set_row_column.h
get_min_max_vec_SOURCES = get_min_max_vec.C get_min_max_vec.h
multiple_rhs_matrix_solve_SOURCES = multiple_rhs_matrix_solve.C multiple_rhs_matrix_solve.h 
power_method_SOURCES = power_method.C power_method.h
inverse_power_method_SOURCES = inverse_power_method.C inverse_power_method.h

# Fortran Tests

TESTS =
TESTS += init.sh
TESTS += get_set_row_column
TESTS += get_min_max_vec
TESTS += multiple_rhs_matrix_solve
TESTS += power_method
TESTS += inverse_power_method
TESTS += finalize.sh

EXTRA_DIST = init.sh finalize.sh input

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

if GIT_CLONE
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(QUESO_STAMPED_FILES)
	echo 'updated source license headers' >$@
endif

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
