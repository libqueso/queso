## Process this file with automake to produce Makefile.in

AM_CPPFLAGS  = $(QUESO_CPPFLAGS)
AM_CPPFLAGS += -I$(top_builddir)/src/core/inc
AM_CPPFLAGS += -I$(top_builddir)/inc
AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(GSL_CFLAGS)
AM_CPPFLAGS += $(GRVY_CFLAGS)

if HAVE_ANN
AM_CPPFLAGS += $(ANN_CFLAGS)
endif

LIBS += $(GSL_LIBS)

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

if LIBMESH_SLEPC_ENABLED
  AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)
  LIBS += $(LIBMESH_LIBS)

##################################
# inifinite dimensional examples #
##################################
operatordir = $(prefix)/examples/infinite_dim/operator
operator_DATA = $(top_srcdir)/examples/infinite_dim/operator/mesh.e
operator_PROGRAMS = operator
operator_SOURCES = infinite_dim/operator/operator.C
operator_LDADD = $(top_builddir)/src/libqueso.la
operator_CPPFLAGS = $(AM_CPPFLAGS)
dist_operator_DATA = $(operator_SOURCES) $(operator_DATA)

gaussian_random_fielddir = $(prefix)/examples/infinite_dim/gaussian_fields
gaussian_random_field_DATA = $(top_srcdir)/examples/infinite_dim/gaussian_fields/mesh.e
gaussian_random_field_PROGRAMS = gaussian_random_field
gaussian_random_field_SOURCES = infinite_dim/gaussian_fields/gaussian_random_field.C
gaussian_random_field_LDADD = $(top_builddir)/src/libqueso.la
gaussian_random_field_CPPFLAGS = $(AM_CPPFLAGS)
dist_gaussian_random_field_DATA = $(gaussian_random_field_SOURCES) $(gaussian_random_field_DATA)

gaussian_random_field_genmeshdir = $(prefix)/examples/infinite_dim/gaussian_fields
gaussian_random_field_genmesh_DATA = $(top_srcdir)/examples/infinite_dim/gaussian_fields/mesh.e
gaussian_random_field_genmesh_PROGRAMS = gaussian_random_field_genmesh
gaussian_random_field_genmesh_SOURCES = infinite_dim/gaussian_fields/gaussian_random_field_genmesh.C
gaussian_random_field_genmesh_LDADD = $(top_builddir)/src/libqueso.la
gaussian_random_field_genmesh_CPPFLAGS = $(AM_CPPFLAGS)
dist_gaussian_random_field_genmesh_DATA = $(gaussian_random_field_genmesh_SOURCES) $(gaussian_random_field_genmesh_DATA)

invprobdir = $(prefix)/examples/infinite_dim
invprob_DATA = infinite_dim/inverse_options
invprob_PROGRAMS = inverse_problem
inverse_problem_SOURCES = infinite_dim/inverse_problem.C
inverse_problem_LDADD = $(top_builddir)/src/libqueso.la
inverse_problem_CPPFLAGS = $(AM_CPPFLAGS)
dist_invprob_DATA = $(inverse_problem_SOURCES) $(invprob_DATA)

# Parallel inverse problem example
parinvprobdir = $(prefix)/examples/infinite_dim
parinvprob_DATA = infinite_dim/parallel_inverse_options
parinvprob_PROGRAMS = parallel_inverse_problem
parallel_inverse_problem_SOURCES = infinite_dim/parallel_inverse_problem.C
parallel_inverse_problem_LDADD = $(top_builddir)/src/libqueso.la
parallel_inverse_problem_CPPFLAGS = $(AM_CPPFLAGS)
dist_parinvprob_DATA = $(parallel_inverse_problem_SOURCES) $(parinvprob_DATA)
endif

#############################
# simpleStatisticalInverseProblem #
###################################
exSimpleStatisticalInverseProblem_gsldir = $(prefix)/examples/simpleStatisticalInverseProblem

exSimpleStatisticalInverseProblem_gsl_DATA = $(top_srcdir)/examples/simpleStatisticalInverseProblem/tests/test_2013_08_26/example.inp
exSimpleStatisticalInverseProblem_gsl_DATA += $(top_srcdir)/examples/simpleStatisticalInverseProblem/tests/test_2013_08_26/simple_ip_plots.m
exSimpleStatisticalInverseProblem_gsl_DATA += $(top_srcdir)/examples/simpleStatisticalInverseProblem/tests/test_2013_08_26/histyy.m

exSimpleStatisticalInverseProblem_gsl_PROGRAMS = exSimpleStatisticalInverseProblem_gsl

exSimpleStatisticalInverseProblem_gsl_SOURCES = simpleStatisticalInverseProblem/src/example_main.C
exSimpleStatisticalInverseProblem_gsl_SOURCES += simpleStatisticalInverseProblem/src/example_likelihood.C
exSimpleStatisticalInverseProblem_gsl_SOURCES += simpleStatisticalInverseProblem/src/example_compute.C

exSimpleStatisticalInverseProblem_gsl_HEADERS = simpleStatisticalInverseProblem/src/example_likelihood.h
exSimpleStatisticalInverseProblem_gsl_HEADERS += simpleStatisticalInverseProblem/src/example_compute.h

exSimpleStatisticalInverseProblem_gsl_LDADD = $(top_builddir)/src/libqueso.la

exSimpleStatisticalInverseProblem_gsl_CPPFLAGS = -I$(top_srcdir)/examples/simpleStatisticalInverseProblem/src $(AM_CPPFLAGS)

# Make sure everything relevant gets added when doing 'make dist'
dist_exSimpleStatisticalInverseProblem_gsl_DATA =
dist_exSimpleStatisticalInverseProblem_gsl_DATA += ${exSimpleStatisticalInverseProblem_gsl_DATA}
dist_exSimpleStatisticalInverseProblem_gsl_DATA += ${exSimpleStatisticalInverseProblem_gsl_HEADERS}
dist_exSimpleStatisticalInverseProblem_gsl_DATA += ${exSimpleStatisticalInverseProblem_gsl_SOURCES}

###################################
# simpleStatisticalForwardProblem #
###################################
exSimpleStatisticalForwardProblem_gsldir = $(prefix)/examples/simpleStatisticalForwardProblem

exSimpleStatisticalForwardProblem_gsl_DATA = $(top_srcdir)/examples/simpleStatisticalForwardProblem/tests/test_2013_08_27/simple_sfp_example.inp
exSimpleStatisticalForwardProblem_gsl_DATA += $(top_srcdir)/examples/simpleStatisticalForwardProblem/tests/test_2013_08_27/simple_fp_plots.m

exSimpleStatisticalForwardProblem_gsl_PROGRAMS = exSimpleStatisticalForwardProblem_gsl

exSimpleStatisticalForwardProblem_gsl_SOURCES = simpleStatisticalForwardProblem/src/simple_sfp_example_main.C
exSimpleStatisticalForwardProblem_gsl_SOURCES += simpleStatisticalForwardProblem/src/simple_sfp_example_qoi.C
exSimpleStatisticalForwardProblem_gsl_SOURCES += simpleStatisticalForwardProblem/src/simple_sfp_example_compute.C

exSimpleStatisticalForwardProblem_gsl_HEADERS = simpleStatisticalForwardProblem/src/simple_sfp_example_qoi.h
exSimpleStatisticalForwardProblem_gsl_HEADERS += simpleStatisticalForwardProblem/src/simple_sfp_example_compute.h

exSimpleStatisticalForwardProblem_gsl_LDADD = $(top_builddir)/src/libqueso.la

exSimpleStatisticalForwardProblem_gsl_CPPFLAGS = -I$(top_srcdir)/examples/simpleStatisticalForwardProblem/src $(AM_CPPFLAGS)

dist_exSimpleStatisticalForwardProblem_gsl_DATA = ${exSimpleStatisticalForwardProblem_gsl_DATA}
dist_exSimpleStatisticalForwardProblem_gsl_DATA += ${exSimpleStatisticalForwardProblem_gsl_HEADERS}
dist_exSimpleStatisticalForwardProblem_gsl_DATA += ${exSimpleStatisticalForwardProblem_gsl_SOURCES}

####################################
# statisticalForwardProblemExample #
####################################
# Kemelli, 8/26/13: removing this example from the distro
# exStatisticalForwardProblem_gsldir = $(prefix)/examples/statisticalForwardProblem
#
# exStatisticalForwardProblem_gsl_DATA = $(top_srcdir)/examples/statisticalForwardProblem/tests/test_2012_11_13/sfp.inp
# exStatisticalForwardProblem_gsl_DATA += $(top_srcdir)/examples/statisticalForwardProblem/tests/test_2012_11_13/sfp_plot.m
#
# exStatisticalForwardProblem_gsl_PROGRAMS = exStatisticalForwardProblem_gsl
#
# exStatisticalForwardProblem_gsl_SOURCES = statisticalForwardProblem/src/exStatisticalForwardProblem_gsl.C
#
# exStatisticalForwardProblem_gsl_HEADERS = statisticalForwardProblem/src/exStatisticalForwardProblem_appl.h
# exStatisticalForwardProblem_gsl_HEADERS += statisticalForwardProblem/src/exStatisticalForwardProblem_qoi.h
#
# exStatisticalForwardProblem_gsl_LDADD = $(top_builddir)/src/libqueso.la
#
# exStatisticalForwardProblem_gsl_CPPFLAGS = -I$(top_srcdir)/examples/statisticalForwardProblem/src $(AM_CPPFLAGS)
#
#dist_exStatisticalForwardProblem_gsl_DATA = ${exStatisticalForwardProblem_gsl_DATA}
#
#${exStatisticalForwardProblem_gsl_SOURCES} \
#																						#${exStatisticalForwardProblem_gsl_HEADERS}

#############################
# statisticalInverseProblem #
#############################
# Kemelli, 8/26/13: removing this example from the distro
# exStatisticalInverseProblem_gsldir = $(prefix)/examples/statisticalInverseProblem
#
# exStatisticalInverseProblem_gsl_DATA = $(top_srcdir)/examples/statisticalInverseProblem/tests/test_2012_11_15/sip.inp
# exStatisticalInverseProblem_gsl_DATA += $(top_srcdir)/examples/statisticalInverseProblem/tests/test_2012_11_15/sip_plot.m
#
# exStatisticalInverseProblem_gsl_PROGRAMS = exStatisticalInverseProblem_gsl
#
# exStatisticalInverseProblem_gsl_SOURCES = statisticalInverseProblem/src/exStatisticalInverseProblem_gsl.C
#
# exStatisticalInverseProblem_gsl_HEADERS = statisticalInverseProblem/src/exStatisticalInverseProblem_appl.h
# exStatisticalInverseProblem_gsl_HEADERS += statisticalInverseProblem/src/exStatisticalInverseProblem_likelihood.h
#
# exStatisticalInverseProblem_gsl_LDADD = $(top_builddir)/src/libqueso.la
#
# exStatisticalInverseProblem_gsl_CPPFLAGS = -I$(top_srcdir)/examples/statisticalInverseProblem/src $(AM_CPPFLAGS)
#
#dist_exStatisticalInverseProblem_gsl_DATA = ${exStatisticalInverseProblem_gsl_DATA} \
#																						#${exStatisticalInverseProblem_gsl_HEADERS} \
#																						#${exStatisticalInverseProblem_gsl_SOURCES}

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

exTgaValidationCycle_gsl_SOURCES = validationCycle/src/exTgaValidationCycle_gsl.C

exTgaValidationCycle_gsl_HEADERS = validationCycle/src/exTgaValidationCycle_appl.h
exTgaValidationCycle_gsl_HEADERS += validationCycle/src/exTgaValidationCycle_likelihood.h
exTgaValidationCycle_gsl_HEADERS += validationCycle/src/exTgaValidationCycle_qoi.h

exTgaValidationCycle_gsl_LDADD = $(top_builddir)/src/libqueso.la

exTgaValidationCycle_gsl_CPPFLAGS = -I$(top_srcdir)/examples/validationCycle/src $(AM_CPPFLAGS)

dist_exTgaValidationCycle_gsl_DATA =
dist_exTgaValidationCycle_gsl_DATA += $(exTgaValidationCycle_gsl_DATA)
dist_exTgaValidationCycle_gsl_DATA +=	${exTgaValidationCycle_gsl_HEADERS}
dist_exTgaValidationCycle_gsl_DATA +=	${exTgaValidationCycle_gsl_SOURCES}

####################
# validationCycle2 #
####################
tga2_gsldir = $(prefix)/examples/validationCycle2

tga2_gslinputdatadir = $(prefix)/examples/validationCycle2/inputData

tga2_gslinputdata_DATA  = $(top_srcdir)/examples/validationCycle2/tests/test_2012_11_15/inputData/scenario_100_K_min.dat
tga2_gslinputdata_DATA += $(top_srcdir)/examples/validationCycle2/tests/test_2012_11_15/inputData/scenario_25_K_min.dat
tga2_gslinputdata_DATA += $(top_srcdir)/examples/validationCycle2/tests/test_2012_11_15/inputData/scenario_50_K_min.dat
tga2_gslinputdata_DATA += $(top_srcdir)/examples/validationCycle2/tests/test_2012_11_15/inputData/scenario_5_K_min.dat
tga2_gsl_DATA  = $(top_srcdir)/examples/validationCycle2/tests/test_2012_11_15/tgaCycle.inp
tga2_gsl_DATA += $(top_srcdir)/examples/validationCycle2/tests/test_2012_11_15/tga_cycle_plot.m

tga2_gsl_PROGRAMS = tga2_gsl

tga2_gsl_SOURCES = validationCycle2/src/tga2_gsl.C
tga2_gsl_SOURCES += validationCycle2/src/tga2_appl.C
tga2_gsl_SOURCES += validationCycle2/src/tga2_likelihood.C
tga2_gsl_SOURCES += validationCycle2/src/tga2_qoi.C
tga2_gsl_SOURCES += validationCycle2/src/tga2_func.C

tga2_gsl_HEADERS = validationCycle2/src/tga2_appl.h
tga2_gsl_HEADERS += validationCycle2/src/tga2_likelihood.h
tga2_gsl_HEADERS += validationCycle2/src/tga2_qoi.h
tga2_gsl_HEADERS += validationCycle2/src/tga2_func.h

tga2_gsl_LDADD = $(top_builddir)/src/libqueso.la

tga2_gsl_CPPFLAGS = -I$(top_srcdir)/examples/validationCycle2/src $(AM_CPPFLAGS)

dist_tga2_gsl_DATA =
dist_tga2_gsl_DATA += $(tga2_gsl_DATA)
dist_tga2_gsl_DATA +=	${tga2_gslinputdata_DATA}
dist_tga2_gsl_DATA +=	${tga2_gsl_HEADERS}
dist_tga2_gsl_DATA +=	${tga2_gsl_SOURCES}

#####################
# infoTheoryProblem #
#####################
exInfoTheory_gsldir = $(prefix)/examples/infoTheoryProblem

exInfoTheory_gsl_DATA = $(top_srcdir)/examples/infoTheoryProblem/tests/test_2011_02_03/infTh.inp

exInfoTheory_gsl_PROGRAMS = exInfoTheory_gsl

exInfoTheory_gsl_SOURCES = infoTheoryProblem/src/exInfoTheory_gsl.C

exInfoTheory_gsl_HEADERS = infoTheoryProblem/src/exInfoTheory_appl.h

exInfoTheory_gsl_LDADD = $(top_builddir)/src/libqueso.la

exInfoTheory_gsl_CPPFLAGS = -I$(top_srcdir)/examples/infoTheoryProblem/src $(AM_CPPFLAGS)

dist_exInfoTheory_gsl_DATA =
dist_exInfoTheory_gsl_DATA += ${exInfoTheory_gsl_DATA}
dist_exInfoTheory_gsl_DATA +=	${exInfoTheory_gsl_HEADERS}
dist_exInfoTheory_gsl_DATA +=	${exInfoTheory_gsl_SOURCES}


####################
# gravity #
####################
gravity_gsldir = $(prefix)/examples/gravity

gravity_gsl_DATA = $(top_srcdir)/examples/gravity/tests/test_2013_01_22/gravity_inv_fwd.inp
gravity_gsl_DATA += $(top_srcdir)/examples/gravity/tests/test_2013_01_22/gravity_plots_ip.m
gravity_gsl_DATA += $(top_srcdir)/examples/gravity/tests/test_2013_01_22/gravity_plots_fp.m
gravity_gsl_DATA += $(top_srcdir)/examples/gravity/tests/test_2013_01_22/autocorr_local.m

gravity_gsl_PROGRAMS = gravity_gsl

gravity_gsl_SOURCES = gravity/src/gravity_main.C
gravity_gsl_SOURCES += gravity/src/gravity_compute.C
gravity_gsl_SOURCES += gravity/src/gravity_likelihood.C
gravity_gsl_SOURCES += gravity/src/gravity_qoi.C

gravity_gsl_HEADERS = gravity/src/gravity_compute.h
gravity_gsl_HEADERS += gravity/src/gravity_likelihood.h
gravity_gsl_HEADERS += gravity/src/gravity_qoi.h

gravity_gsl_LDADD = $(top_builddir)/src/libqueso.la

gravity_gsl_CPPFLAGS = -I$(top_srcdir)/examples/gravity/src $(AM_CPPFLAGS)

dist_gravity_gsl_DATA =
dist_gravity_gsl_DATA += $(gravity_gsl_DATA)
dist_gravity_gsl_DATA += ${gravity_gsl_HEADERS}
dist_gravity_gsl_DATA += ${gravity_gsl_SOURCES}

####################
# line2D #
####################
line2D_gsldir = $(prefix)/examples/line2D_with_Sensitivity

line2D_gsl_DATA = $(top_srcdir)/examples/line2D_with_Sensitivity/slope_inv_fwd.inp
line2D_gsl_DATA += $(top_srcdir)/examples/line2D_with_Sensitivity/README
line2D_gsl_DATA += $(top_srcdir)/examples/line2D_with_Sensitivity/gen_plots.py
line2D_gsl_DATA += $(top_srcdir)/examples/line2D_with_Sensitivity/gen_truth.py

line2D_gsl_PROGRAMS = line2D_gsl

line2D_gsl_SOURCES = line2D_with_Sensitivity/sensitivity_c.C
line2D_gsl_SOURCES += line2D_with_Sensitivity/sensitivity_m.C
line2D_gsl_SOURCES += line2D_with_Sensitivity/sensitivity_mc.C
line2D_gsl_SOURCES += line2D_with_Sensitivity/slope_compute.C
line2D_gsl_SOURCES += line2D_with_Sensitivity/slope_likelihood.C
line2D_gsl_SOURCES += line2D_with_Sensitivity/slope_main.C
line2D_gsl_SOURCES += line2D_with_Sensitivity/slope_qoi.C

line2D_gsl_HEADERS = line2D_with_Sensitivity/sensitivity_c.h
line2D_gsl_HEADERS += line2D_with_Sensitivity/sensitivity_mc.h
line2D_gsl_HEADERS += line2D_with_Sensitivity/sensitivity_m.h
line2D_gsl_HEADERS += line2D_with_Sensitivity/slope_compute.h
line2D_gsl_HEADERS += line2D_with_Sensitivity/slope_likelihood.h
line2D_gsl_HEADERS += line2D_with_Sensitivity/slope_qoi.h

line2D_gsl_LDADD = $(top_builddir)/src/libqueso.la

line2D_gsl_CPPFLAGS = -I$(top_srcdir)/examples/line2D_with_Sensitivity $(AM_CPPFLAGS)

dist_line2D_gsl_DATA =
dist_line2D_gsl_DATA += $(line2D_gsl_DATA)
dist_line2D_gsl_DATA += ${line2D_gsl_HEADERS}
dist_line2D_gsl_DATA += ${line2D_gsl_SOURCES}

####################
# bimodal          #
####################
bimodal_gsldir = $(prefix)/examples/bimodal

bimodal_gsl_DATA = $(top_srcdir)/examples/bimodal/tests/test_2013_12_02/bimodal_1chain.inp
bimodal_gsl_DATA += $(top_srcdir)/examples/bimodal/matlab/plot_all.m
bimodal_gsl_DATA += $(top_srcdir)/examples/bimodal/matlab/plot_likelihood.m
bimodal_gsl_DATA += $(top_srcdir)/examples/bimodal/matlab/plot_likelihood_normalized_taus.m
bimodal_gsl_DATA += $(top_srcdir)/examples/bimodal/matlab/plot_likelihood_unnormalized_taus.m

bimodal_gsl_PROGRAMS = bimodal_gsl

bimodal_gsl_SOURCES = bimodal/src/bimodal_main.C
bimodal_gsl_SOURCES += bimodal/src/bimodal_compute.C
bimodal_gsl_SOURCES += bimodal/src/bimodal_likelihood.C

bimodal_gsl_HEADERS = bimodal/src/bimodal_compute.h
bimodal_gsl_HEADERS += bimodal/src/bimodal_likelihood.h

bimodal_gsl_LDADD = $(top_builddir)/src/libqueso.la

bimodal_gsl_CPPFLAGS = -I$(top_srcdir)/examples/bimodal/src $(AM_CPPFLAGS)

dist_bimodal_gsl_DATA =
dist_bimodal_gsl_DATA += $(bimodal_gsl_DATA)
dist_bimodal_gsl_DATA += ${bimodal_gsl_HEADERS}
dist_bimodal_gsl_DATA += ${bimodal_gsl_SOURCES}

####################
# hysteretic          #
####################
hysteretic_gsldir = $(prefix)/examples/hysteretic

hysteretic_gsl_DATA = $(top_srcdir)/examples/hysteretic/tests/test_2013_12_11/example.inp
hysteretic_gsl_DATA += $(top_srcdir)/examples/hysteretic/tests/test_2013_12_11/an.txt
hysteretic_gsl_DATA += $(top_srcdir)/examples/hysteretic/tests/test_2013_12_11/measured_data1_1.txt
hysteretic_gsl_DATA += $(top_srcdir)/examples/hysteretic/tests/test_2013_12_11/measured_data1_2.txt
hysteretic_gsl_DATA += $(top_srcdir)/examples/hysteretic/tests/test_2013_12_11/measured_data1_3.txt
hysteretic_gsl_DATA += $(top_srcdir)/examples/hysteretic/tests/test_2013_12_11/measured_data1_4.txt
hysteretic_gsl_DATA += $(top_srcdir)/examples/hysteretic/matlab/plot_all.m
hysteretic_gsl_DATA += $(top_srcdir)/examples/hysteretic/matlab/cpp_gen.m

hysteretic_gsl_PROGRAMS = hysteretic_gsl

hysteretic_gsl_SOURCES = hysteretic/src/example_main.C
hysteretic_gsl_SOURCES += hysteretic/src/example_compute.C
hysteretic_gsl_SOURCES += hysteretic/src/example_likelihood.C
hysteretic_gsl_SOURCES += hysteretic/src/example_hyst.C

hysteretic_gsl_HEADERS = hysteretic/src/example_compute.h
hysteretic_gsl_HEADERS += hysteretic/src/example_likelihood.h
hysteretic_gsl_HEADERS += hysteretic/src/example_hyst.h

hysteretic_gsl_LDADD = $(top_builddir)/src/libqueso.la

hysteretic_gsl_CPPFLAGS = -I$(top_srcdir)/examples/hysteretic/src $(AM_CPPFLAGS)

dist_hysteretic_gsl_DATA =
dist_hysteretic_gsl_DATA += $(hysteretic_gsl_DATA)
dist_hysteretic_gsl_DATA += ${hysteretic_gsl_HEADERS}
dist_hysteretic_gsl_DATA += ${hysteretic_gsl_SOURCES}

##################
# GPMSA examples #
##################
gpmsa_scalardir = $(prefix)/examples/gp/scalar

gpmsa_scalar_PROGRAMS = gpmsa_scalar
gpmsa_scalar_SOURCES = gp/scalar/gpmsa_scalar.C
gpmsa_scalar_LDADD = $(top_builddir)/src/libqueso.la
gpmsa_scalar_CPPFLAGS = -I$(top_srcdir)/examples/gp/scalar $(AM_CPPFLAGS)

dist_gpmsa_scalar_DATA =
dist_gpmsa_scalar_DATA += ${gpmsa_scalar_SOURCES}
dist_gpmsa_scalar_DATA += gp/scalar/dakota_pstudy.dat
dist_gpmsa_scalar_DATA += gp/scalar/gpmsa_input.txt
dist_gpmsa_scalar_DATA += gp/scalar/ctf_dat.txt

gpmsa_pseudovectordir = $(prefix)/examples/gp/pseudovector

gpmsa_pseudovector_PROGRAMS = gpmsa_pseudovector
gpmsa_pseudovector_SOURCES = gp/pseudovector/gpmsa_pseudovector.C
gpmsa_pseudovector_LDADD = $(top_builddir)/src/libqueso.la
gpmsa_pseudovector_CPPFLAGS = -I$(top_srcdir)/examples/gp/pseudovector $(AM_CPPFLAGS)

dist_gpmsa_pseudovector_DATA =
dist_gpmsa_pseudovector_DATA += ${gpmsa_pseudovector_SOURCES}
dist_gpmsa_pseudovector_DATA += gp/scalar/dakota_pstudy.dat
dist_gpmsa_pseudovector_DATA += gp/pseudovector/gpmsa_input.txt
dist_gpmsa_pseudovector_DATA += gp/scalar/ctf_dat.txt

gpmsa_vectordir = $(prefix)/examples/gp/vector

gpmsa_vector_PROGRAMS = gpmsa_vector
gpmsa_vector_SOURCES = gp/vector/gpmsa_vector.C
gpmsa_vector_LDADD = $(top_builddir)/src/libqueso.la
gpmsa_vector_CPPFLAGS = -I$(top_srcdir)/examples/gp/vector $(AM_CPPFLAGS)

dist_gpmsa_vector_DATA =
dist_gpmsa_vector_DATA += ${gpmsa_vector_SOURCES}
dist_gpmsa_vector_DATA += gp/scalar/dakota_pstudy.dat
dist_gpmsa_vector_DATA += gp/vector/gpmsa_input.txt
dist_gpmsa_vector_DATA += gp/scalar/ctf_dat.txt


#######################################
# Canned Gaussian likelihood examples #
#######################################
scalarCovariancedir = $(prefix)/examples/gaussian_likelihoods
scalarCovariance_PROGRAMS = scalarCovariance
scalarCovariance_SOURCES = gaussian_likelihoods/scalarCovariance.C
scalarCovariance_LDADD = $(top_builddir)/src/libqueso.la
scalarCovariance_CPPFLAGS = -I$(top_srcdir)/examples/gaussian_likelihoods/scalarCovariance $(AM_CPPFLAGS)

diagonalCovariancedir = $(prefix)/examples/gaussian_likelihoods
diagonalCovariance_PROGRAMS = diagonalCovariance
diagonalCovariance_SOURCES = gaussian_likelihoods/diagonalCovariance.C
diagonalCovariance_LDADD = $(top_builddir)/src/libqueso.la
diagonalCovariance_CPPFLAGS = -I$(top_srcdir)/examples/gaussian_likelihoods/diagonalCovariance $(AM_CPPFLAGS)

fullCovariancedir = $(prefix)/examples/gaussian_likelihoods
fullCovariance_PROGRAMS = fullCovariance
fullCovariance_SOURCES = gaussian_likelihoods/fullCovariance.C
fullCovariance_LDADD = $(top_builddir)/src/libqueso.la
fullCovariance_CPPFLAGS = -I$(top_srcdir)/examples/gaussian_likelihoods/fullCovariance $(AM_CPPFLAGS)

blockDiagonalCovariancedir = $(prefix)/examples/gaussian_likelihoods
blockDiagonalCovariance_PROGRAMS = blockDiagonalCovariance
blockDiagonalCovariance_SOURCES = gaussian_likelihoods/blockDiagonalCovariance.C
blockDiagonalCovariance_LDADD = $(top_builddir)/src/libqueso.la
blockDiagonalCovariance_CPPFLAGS = -I$(top_srcdir)/examples/gaussian_likelihoods/blockDiagonalCovariance $(AM_CPPFLAGS)

fullCovarianceRandomCoefficientdir = $(prefix)/examples/gaussian_likelihoods
fullCovarianceRandomCoefficient_PROGRAMS = fullCovarianceRandomCoefficient
fullCovarianceRandomCoefficient_SOURCES = gaussian_likelihoods/fullCovarianceRandomCoefficient.C
fullCovarianceRandomCoefficient_LDADD = $(top_builddir)/src/libqueso.la
fullCovarianceRandomCoefficient_CPPFLAGS = -I$(top_srcdir)/examples/gaussian_likelihoods/fullCovarianceRandomCoefficient $(AM_CPPFLAGS)

blockDiagonalCovarianceRandomCoefficientsdir = $(prefix)/examples/gaussian_likelihoods
blockDiagonalCovarianceRandomCoefficients_PROGRAMS = blockDiagonalCovarianceRandomCoefficients
blockDiagonalCovarianceRandomCoefficients_SOURCES = gaussian_likelihoods/blockDiagonalCovarianceRandomCoefficients.C
blockDiagonalCovarianceRandomCoefficients_LDADD = $(top_builddir)/src/libqueso.la
blockDiagonalCovarianceRandomCoefficients_CPPFLAGS = -I$(top_srcdir)/examples/gaussian_likelihoods/blockDiagonalCovarianceRandomCoefficients $(AM_CPPFLAGS)

####################################
# Interpolation Surrogate examples #
####################################
interpsurrgdir = $(prefix)/examples/interpolation_surrogate

interpsurrg_PROGRAMS = 4d_interp 4d_interp_read

4d_interp_SOURCES = interpolation_surrogate/4d_interp.C
4d_interp_LDADD = $(top_builddir)/src/libqueso.la
4d_interp_CPPFLAGS = $(AM_CPPFLAGS)

4d_interp_read_SOURCES = interpolation_surrogate/4d_interp_read.C
4d_interp_read_LDADD = $(top_builddir)/src/libqueso.la
4d_interp_read_CPPFLAGS = $(AM_CPPFLAGS)

dist_interpsurrg_DATA =
dist_interpsurrg_DATA += ${4d_interp_SOURCES}
dist_interpsurrg_DATA += ${4d_interp_read_SOURCES}
dist_interpsurrg_DATA += interpolation_surrogate/input.in
dist_interpsurrg_DATA += interpolation_surrogate/4d_interp_data_coarse.dat

#####################
# SharedPtr example #
#####################
pointersdir = $(prefix)/examples/pointers

pointers_PROGRAMS = example_SharedPtr

example_SharedPtr_SOURCES = pointers/example_SharedPtr.C
example_SharedPtr_LDADD = $(top_builddir)/src/libqueso.la
example_SharedPtr_CPPFLAGS = $(AM_CPPFLAGS)

dist_pointers_DATA =
dist_pointers_DATA += ${example_SharedPtr_SOURCES}

####################################
# Custom Transition Kernel Example #
####################################
custom_tkdir = $(prefix)/examples/custom_tk

custom_tk_PROGRAMS = hessian_eg

hessian_eg_SOURCES =
hessian_eg_SOURCES += custom_tk/hessian_eg.C
hessian_eg_SOURCES += custom_tk/hessian_eg_main.C
custom_tk_HEADERS = custom_tk/hessian_eg.h
hessian_eg_LDADD = $(top_builddir)/src/libqueso.la
hessian_eg_CPPFLAGS = $(AM_CPPFLAGS)

dist_custom_tk_DATA =
dist_custom_tk_DATA += ${hessian_eg_SOURCES}
dist_custom_tk_DATA += ${hessian_eg_HEADERS}


if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

runalldir = $(prefix)/examples
dist_runall_SCRIPTS = run_all.sh

MAINTAINERCLEANFILES = Makefile.in
