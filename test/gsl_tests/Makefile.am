if SVN_CHECKOUT
BUILT_SOURCES = .license.stamp
endif

check_PROGRAMS = # Append below
check_PROGRAMS += get_set_row_column \
                  get_min_max_vec \
	          multiple_rhs_matrix_solve \
                  power_method \
                  inverse_power_method

AM_CPPFLAGS  = -I$(top_srcdir)/src -I$(top_builddir)/src
LIBS         = -L$(top_builddir)/src/ -lqueso 

AM_CPPFLAGS  += -I. \
                -I$(top_srcdir)/src/core/inc \
                -I$(top_srcdir)/src/misc/inc \
                -I$(top_srcdir)/src/basic/inc \
                -I$(top_srcdir)/src/stats/inc \
                -I$(top_srcdir)/src/gp/inc \
	              -I$(top_srcdir)/inc \
	                $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(ANN_CFLAGS)

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

TESTS             = init.sh \
                    get_set_row_column \
                    get_min_max_vec \
	            multiple_rhs_matrix_solve \
                    power_method \
                    inverse_power_method \
		    finalize.sh

EXTRA_DIST = init.sh finalize.sh input

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

if SVN_CHECKOUT
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(srcdir)/*.C $(srcdir)/*.h
	echo 'updated source license headers' >$@
endif

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
