## Process this file with automake to produce Makefile.in

QUESO_CPPFLAGS  = -I$(top_srcdir)/src/core/inc
QUESO_CPPFLAGS += -I$(top_srcdir)/src/misc/inc
QUESO_CPPFLAGS += -I$(top_srcdir)/src/basic/inc
QUESO_CPPFLAGS += -I$(top_srcdir)/src/stats/inc
QUESO_CPPFLAGS += $(BOOST_CPPFLAGS)
QUESO_CPPFLAGS += $(GSL_CFLAGS)
QUESO_CPPFLAGS += $(GRVY_CFLAGS)
QUESO_CPPFLAGS += $(ANN_CFLAGS)

LIBS += $(GSL_LIBS)

if TRILINOS_ENABLED
  QUESO_CPPFLAGS += -I$(TRILINOS_INCLUDE) $(EPETRA_INCLUDES)
  LIBS += -lteuchos -lepetra
endif

if GLPK_ENABLED	
  QUESO_CPPFLAGS += $(GLPK_CFLAGS)
  LIBS += $(GLPK_LIBS)
endif

if HDF5_ENABLED	
  QUESO_CPPFLAGS += $(HDF5_CFLAGS) 
  LIBS += $(HDF5_LIBS)
endif

####################################
# statisticalForwardProblemExample #
####################################
exStatisticalForwardProblem_gsldir = $(prefix)/examples/statisticalForwardProblem

exStatisticalForwardProblem_gsl_DATA = $(top_srcdir)/examples/statisticalForwardProblem/tests/test_2012_11_13/sfp.inp
exStatisticalForwardProblem_gsl_DATA += $(top_srcdir)/examples/statisticalForwardProblem/tests/test_2012_11_13/sfp_plot.m

exStatisticalForwardProblem_gsl_PROGRAMS = exStatisticalForwardProblem_gsl

exStatisticalForwardProblem_gsl_SOURCES = $(top_srcdir)/examples/statisticalForwardProblem/src/exStatisticalForwardProblem_gsl.C

exStatisticalForwardProblem_gsl_HEADERS = $(top_srcdir)/examples/statisticalForwardProblem/src/exStatisticalForwardProblem_appl.h
exStatisticalForwardProblem_gsl_HEADERS += $(top_srcdir)/examples/statisticalForwardProblem/src/exStatisticalForwardProblem_qoi.h

exStatisticalForwardProblem_gsl_LDADD = $(top_builddir)/src/libqueso.la

exStatisticalForwardProblem_gsl_CPPFLAGS = -I$(top_srcdir)/examples/statisticalForwardProblem/src $(QUESO_CPPFLAGS)

# Make sure everything relevant gets added when doing 'make dist'
dist_exStatisticalForwardProblem_gsl_DATA = ${exStatisticalForwardProblem_gsl_DATA} \
																						${exStatisticalForwardProblem_gsl_SOURCES} \
																						${exStatisticalForwardProblem_gsl_HEADERS}

#############################
# statisticalInverseProblem #
#############################
exStatisticalInverseProblem_gsldir = $(prefix)/examples/statisticalInverseProblem

exStatisticalInverseProblem_gsl_DATA = $(top_srcdir)/examples/statisticalInverseProblem/tests/test_2012_11_15/sip.inp
exStatisticalInverseProblem_gsl_DATA += $(top_srcdir)/examples/statisticalInverseProblem/tests/test_2012_11_15/sip_plot.m

exStatisticalInverseProblem_gsl_PROGRAMS = exStatisticalInverseProblem_gsl

exStatisticalInverseProblem_gsl_SOURCES = $(top_srcdir)/examples/statisticalInverseProblem/src/exStatisticalInverseProblem_gsl.C

exStatisticalInverseProblem_gsl_HEADERS = $(top_srcdir)/examples/statisticalInverseProblem/src/exStatisticalInverseProblem_appl.h
exStatisticalInverseProblem_gsl_HEADERS += $(top_srcdir)/examples/statisticalInverseProblem/src/exStatisticalInverseProblem_likelihood.h

exStatisticalInverseProblem_gsl_LDADD = $(top_builddir)/src/libqueso.la

exStatisticalInverseProblem_gsl_CPPFLAGS = -I$(top_srcdir)/examples/statisticalInverseProblem/src $(QUESO_CPPFLAGS)

dist_exStatisticalInverseProblem_gsl_DATA = ${exStatisticalInverseProblem_gsl_DATA} \
																						${exStatisticalInverseProblem_gsl_HEADERS} \
																						${exStatisticalInverseProblem_gsl_SOURCES}

###################
# validationCycle #
###################
exTgaValidationCycle_gsldir = $(prefix)/examples/validationCycle

exTgaValidationCycle_gsl_DATA = $(top_srcdir)/examples/validationCycle/tests/test_2012_11_15/inputData/scenario_100_K_min.dat
exTgaValidationCycle_gsl_DATA += $(top_srcdir)/examples/validationCycle/tests/test_2012_11_15/inputData/scenario_25_K_min.dat
exTgaValidationCycle_gsl_DATA += $(top_srcdir)/examples/validationCycle/tests/test_2012_11_15/inputData/scenario_50_K_min.dat
exTgaValidationCycle_gsl_DATA += $(top_srcdir)/examples/validationCycle/tests/test_2012_11_15/inputData/scenario_5_K_min.dat
exTgaValidationCycle_gsl_DATA += $(top_srcdir)/examples/validationCycle/tests/test_2012_11_15/tgaCycle.inp
exTgaValidationCycle_gsl_DATA += $(top_srcdir)/examples/validationCycle/tests/test_2012_11_15/tgaCycle_2009_03_30.inp
exTgaValidationCycle_gsl_DATA += $(top_srcdir)/examples/validationCycle/tests/test_2012_11_15/tga_cycle_plot.m

exTgaValidationCycle_gsl_PROGRAMS = exTgaValidationCycle_gsl

exTgaValidationCycle_gsl_SOURCES = $(top_srcdir)/examples/validationCycle/src/exTgaValidationCycle_gsl.C

exTgaValidationCycle_gsl_HEADERS = $(top_srcdir)/examples/validationCycle/src/exTgaValidationCycle_appl.h
exTgaValidationCycle_gsl_HEADERS += $(top_srcdir)/examples/validationCycle/src/exTgaValidationCycle_likelihood.h
exTgaValidationCycle_gsl_HEADERS += $(top_srcdir)/examples/validationCycle/src/exTgaValidationCycle_qoi.h

exTgaValidationCycle_gsl_LDADD = $(top_builddir)/src/libqueso.la

exTgaValidationCycle_gsl_CPPFLAGS = -I$(top_srcdir)/examples/validationCycle/src $(QUESO_CPPFLAGS)

dist_exTgaValidationCycle_gsl_DATA = $(exTgaValidationCycle_gsl_DATA) \
																		 ${exTgaValidationCycle_gsl_HEADERS} \
																		 ${exTgaValidationCycle_gsl_SOURCES}

####################
# validationCycle2 #
####################
tga2_gsldir = $(prefix)/examples/validationCycle2

tga2_gsl_DATA = $(top_srcdir)/examples/validationCycle2/tests/test_2012_11_15/inputData/scenario_100_K_min.dat
tga2_gsl_DATA += $(top_srcdir)/examples/validationCycle2/tests/test_2012_11_15/inputData/scenario_25_K_min.dat
tga2_gsl_DATA += $(top_srcdir)/examples/validationCycle2/tests/test_2012_11_15/inputData/scenario_50_K_min.dat
tga2_gsl_DATA += $(top_srcdir)/examples/validationCycle2/tests/test_2012_11_15/inputData/scenario_5_K_min.dat
tga2_gsl_DATA += $(top_srcdir)/examples/validationCycle2/tests/test_2012_11_15/tgaCycle.inp
tga2_gsl_DATA += $(top_srcdir)/examples/validationCycle2/tests/test_2012_11_15/tga_cycle_plot.m

tga2_gsl_PROGRAMS = tga2_gsl

tga2_gsl_SOURCES = $(top_srcdir)/examples/validationCycle2/src/tga2_gsl.C
tga2_gsl_SOURCES += $(top_srcdir)/examples/validationCycle2/src/tga2_appl.C
tga2_gsl_SOURCES += $(top_srcdir)/examples/validationCycle2/src/tga2_likelihood.C
tga2_gsl_SOURCES += $(top_srcdir)/examples/validationCycle2/src/tga2_qoi.C
tga2_gsl_SOURCES += $(top_srcdir)/examples/validationCycle2/src/tga2_func.C

tga2_gsl_HEADERS = $(top_srcdir)/examples/validationCycle2/src/tga2_appl.h
tga2_gsl_HEADERS += $(top_srcdir)/examples/validationCycle2/src/tga2_likelihood.h
tga2_gsl_HEADERS += $(top_srcdir)/examples/validationCycle2/src/tga2_qoi.h
tga2_gsl_HEADERS += $(top_srcdir)/examples/validationCycle2/src/tga2_func.h

tga2_gsl_LDADD = $(top_builddir)/src/libqueso.la

tga2_gsl_CPPFLAGS = -I$(top_srcdir)/examples/validationCycle2/src $(QUESO_CPPFLAGS)

dist_tga2_gsl_DATA = $(tga2_gsl_DATA) \
										 ${tga2_gsl_HEADERS} \
										 ${tga2_gsl_SOURCES}

#####################
# infoTheoryProblem #
#####################
exInfoTheory_gsldir = $(prefix)/examples/infoTheoryProblem

exInfoTheory_gsl_DATA = $(top_srcdir)/examples/infoTheoryProblem/tests/test_2011_02_03/infTh.inp

exInfoTheory_gsl_PROGRAMS = exInfoTheory_gsl

exInfoTheory_gsl_SOURCES = $(top_srcdir)/examples/infoTheoryProblem/src/exInfoTheory_gsl.C

exInfoTheory_gsl_HEADERS = $(top_srcdir)/examples/infoTheoryProblem/src/exInfoTheory_appl.h

exInfoTheory_gsl_LDADD = $(top_builddir)/src/libqueso.la

exInfoTheory_gsl_CPPFLAGS = -I$(top_srcdir)/examples/infoTheoryProblem/src $(QUESO_CPPFLAGS)

dist_exInfoTheory_gsl_DATA = ${exInfoTheory_gsl_DATA} \
														 ${exInfoTheory_gsl_HEADERS} \
														 ${exInfoTheory_gsl_SOURCES}


####################
# gravity #
####################
gravity_gsldir = $(prefix)/examples/gravity

gravity_gsl_DATA = $(top_srcdir)/examples/gravity/tests/test_2013_01_22/gravity_inv_fwd.inp
gravity_gsl_DATA += $(top_srcdir)/examples/gravity/tests/test_2013_01_22/gravity_plots_ip.m
gravity_gsl_DATA += $(top_srcdir)/examples/gravity/tests/test_2013_01_22/gravity_plots_fp.m
gravity_gsl_DATA += $(top_srcdir)/examples/gravity/tests/test_2013_01_22/autocorr_local.m
gravity_gsl_DATA += $(top_srcdir)/examples/gravity/src/Makefile_example_violeta

gravity_gsl_PROGRAMS = gravity_gsl

gravity_gsl_SOURCES = $(top_srcdir)/examples/gravity/src/gravity_main.C
gravity_gsl_SOURCES += $(top_srcdir)/examples/gravity/src/gravity_compute.C
gravity_gsl_SOURCES += $(top_srcdir)/examples/gravity/src/gravity_likelihood.C
gravity_gsl_SOURCES += $(top_srcdir)/examples/gravity/src/gravity_qoi.C

gravity_gsl_HEADERS = $(top_srcdir)/examples/gravity/src/gravity_compute.h
gravity_gsl_HEADERS += $(top_srcdir)/examples/gravity/src/gravity_likelihood.h
gravity_gsl_HEADERS += $(top_srcdir)/examples/gravity/src/gravity_qoi.h

gravity_gsl_LDADD = $(top_builddir)/src/libqueso.la

gravity_gsl_CPPFLAGS = -I$(top_srcdir)/examples/gravity/src $(QUESO_CPPFLAGS)

dist_gravity_gsl_DATA = $(gravity_gsl_DATA) \
										 ${gravity_gsl_HEADERS} \
										 ${gravity_gsl_SOURCES}
										 
if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

MAINTAINERCLEANFILES = Makefile.in
