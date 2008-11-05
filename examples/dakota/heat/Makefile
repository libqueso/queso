ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

ifeq ($(BOOST_INCLUDE_PATH_FOR_MCMC_TOOL),)
override BOOST_INCLUDE_PATH_FOR_MCMC_TOOL = .
endif
ifeq ($(GSL_INCLUDE_PATH_FOR_MCMC_TOOL),)
override GSL_INCLUDE_PATH_FOR_MCMC_TOOL = .
endif
ifeq ($(BOOST_LIB_PATH_FOR_MCMC_TOOL),)
override BOOST_LIB_PATH_FOR_MCMC_TOOL = .
endif
ifeq ($(GSL_LIB_PATH_FOR_MCMC_TOOL),)
override GSL_LIB_PATH_FOR_MCMC_TOOL = .
endif

ifeq ($(UQBT),gsl)

DEFINES = \
	-D__EXAMPLE_USES_GSL__
INC_PATHS = \
	-I. \
	-I../inc \
	-I../../../libs/basic/inc \
	-I../../../libs/fp/inc \
	-I$(MPI_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(TRILINOS_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(BOOST_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(GSL_INCLUDE_PATH_FOR_MCMC_TOOL)
LIB_PATHS = \
	-L../../../libs/basic/lib \
	-L../../../libs/fp/lib \
	-L$(MPI_LIB_PATH_FOR_MCMC_TOOL) \
	-L$(TRILINOS_LIB_PATH_FOR_MCMC_TOOL) \
	-L$(BOOST_LIB_PATH_FOR_MCMC_TOOL) \
	-L$(GSL_LIB_PATH_FOR_MCMC_TOOL)
LIBS = \
	-luq_fp_any \
	-luq_basic_any \
	-luq_basic_$(UQBT) \
	-lepetra \
	-lteuchos \
	-lmpi_cxx \
	-lmpi \
	-lgsl \
	-lgslcblas \
	-lblas \
	-lboost_program_options-$(COMPILER_SPEC_IN_BOOST_LIB_NAME)-mt-1_35

else

DEFINES = \
	-D__APPL_USES_TRILINOS__ -D__UQ_USES_TRILINOS__
INC_PATHS = \
	-I. \
	-I../inc \
	-I../../../libs/basic/inc \
	-I../../../libs/fp/inc \
	-I$(MPI_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(TRILINOS_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(BOOST_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(GSL_INCLUDE_PATH_FOR_MCMC_TOOL)

LIB_PATHS = \
	-L../../../libs/basic/lib \
	-L../../../libs/fp/lib \
	-L$(MPI_LIB_PATH_FOR_MCMC_TOOL) \
	-L$(TRILINOS_LIB_PATH_FOR_MCMC_TOOL) \
	-L$(BOOST_LIB_PATH_FOR_MCMC_TOOL) \
	-L$(GSL_LIB_PATH_FOR_MCMC_TOOL)

LIBS = \
	-luq_fp_any \
	-luq_basic_$(UQBT) \
	-luq_basic_any \
	-lepetra \
	-lteuchos \
	-lmpi_cxx \
	-lmpi \
	-lgsl \
	-lgslcblas \
	-lblas \
	-lboost_program_options-$(COMPILER_SPEC_IN_BOOST_LIB_NAME)-mt-1_35

endif


CXX = g++
CXXFLAGS = -Wall -c

all:	uqHeatEx

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_lib clean_exec

clean_lib:
	cd ../../../libs/basic/src; make clean_all UQBT=$(UQBT)
	cd ../../../libs/fp/src; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~
	rm -f *.o
	rm -f uqHeatEx_$(UQBT)

uqHeatEx: libs uqHeatEx.o 
	$(CXX) $(LIB_PATHS) uqHeatEx.o -o uqHeatEx_$(UQBT) $(LIBS)

libs:
	cd ../../../libs/basic/src; make libuq_basic_all UQBT=$(UQBT)
	cd ../../../libs/fp/src; make libuq_fp_all UQBT=$(UQBT)

.SUFFIXES: .o .C

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<
