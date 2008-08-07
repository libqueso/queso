ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

ifeq ($(UQBT),gsl)

DEFINES = \
	-D__APPL_USES_GSL__ -D__UQ_USES_GSL__
INC_PATHS = \
	-I. \
	-I../inc \
	-I../../../libs/basic/inc \
	-I../../../libs/mcmc/inc \
	-I$(BOOST_INCLUDE_PATH_FOR_PECOS_TOOLKIT) \
	-I$(GSL_INCLUDE_PATH_FOR_PECOS_TOOLKIT)
LIB_PATHS = \
	-L../../../libs/basic/lib \
	-L../../../libs/mcmc/lib \
	-L$(BOOST_LIB_PATH_FOR_PECOS_TOOLKIT) \
	-L$(GSL_LIB_PATH_FOR_PECOS_TOOLKIT)
LIBS = \
	-luq_basic_$(UQBT) \
	-luq_basic_any \
	-luq_mcmc_$(UQBT) \
	-luq_mcmc_any \
	-lgsl \
	-lgslcblas \
	-lboost_program_options-$(COMPILER_SPEC_IN_BOOST_LIB_NAME)-mt-1_35

else

DEFINES = \
	-D__APPL_USES_TRILINOS__ -D__UQ_USES_TRILINOS__
INC_PATHS = \
	-I. \
	-I../../../inc \
	-I$(MPI_INCLUDE_PATH_FOR_PECOS_TOOLKIT) \
	-I$(TRILINOS_INCLUDE_PATH_FOR_PECOS_TOOLKIT) \
	-I$(BOOST_INCLUDE_PATH_FOR_PECOS_TOOLKIT) \
	-I$(GSL_INCLUDE_PATH_FOR_PECOS_TOOLKIT)

LIB_PATHS = \
	-L../../../lib \
	-L$(MPI_LIB_PATH_FOR_PECOS_TOOLKIT) \
	-L$(TRILINOS_LIB_PATH_FOR_PECOS_TOOLKIT) \
	-L$(BOOST_LIB_PATH_FOR_PECOS_TOOLKIT) \
	-L$(GSL_LIB_PATH_FOR_PECOS_TOOLKIT)

LIBS = \
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

all:	uqNormalEx

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_lib clean_exec

clean_lib:
	cd ../../../libs/basic/src; make clean_all UQBT=$(UQBT)
	cd ../../../libs/mcmc/src; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~
	rm -f *.o
	rm -f uqNormalEx_$(UQBT)

uqNormalEx: libs uqNormalEx_$(UQBT).o
	$(CXX) $(LIB_PATHS) uqNormalEx_$(UQBT).o -o uqNormalEx_$(UQBT) $(LIBS)

libs:
	cd ../../../libs/basic/src; make libuq_basic_all UQBT=$(UQBT)
	cd ../../../libs/mcmc/src; make libuq_mcmc_all UQBT=$(UQBT)

.SUFFIXES: .o .C

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<
