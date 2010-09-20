check_PROGRAMS = # Append below
check_PROGRAMS += get_set_row_column \
                  get_min_max_vec \
	          multiple_rhs_matrix_solve \
                  power_method \
                  inverse_power_method

INCLUDES     = -I$(top_srcdir)/src -I$(top_builddir)/src

LIBS         = -L$(top_builddir)/src/.libs -lqueso \
               $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS) ${HDF5_LIBS} -lhdf5

AM_LDFLAGS   = $(BOOST_LDFLAGS) ${GRVY_LIBS} -lgrvy ${GLPK_LIBS} -lglpk ${HDF5_LIBS} -lhdf5

AM_CPPFLAGS  = -I. \
	       -I$(top_srcdir)/src/core/inc \
               -I$(top_srcdir)/src/misc/inc \
               -I$(top_srcdir)/src/basic/inc \
               -I$(top_srcdir)/src/stats/inc \
               -I$(TRILINOS_INCLUDE) \
                $(EPETRA_INCLUDES) $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(GRVY_CFLAGS) $(GLPK_CFLAGS) $(HDF5_CFLAGS)

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

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
