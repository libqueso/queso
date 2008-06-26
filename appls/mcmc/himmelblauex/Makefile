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
	-I../../../libs/basic/inc \
	-I../../../libs/mcmc/inc \
	-I/h2/prudenci/Installations/Boost_1_35_0/include/boost-1_35
LIB_PATHS = \
	-L../../../libs/basic/lib \
	-L../../../libs/mcmc/lib \
	-L/h2/prudenci/Installations/Boost_1_35_0/lib
LIBS = \
	-luq_basic_$(UQBT) \
	-luq_basic_any \
	-luq_mcmc_$(UQBT) \
	-luq_mcmc_any \
	-lgsl \
	-lgslcblas \
	-lboost_program_options-gcc41-mt-1_35

else

DEFINES = \
	-D__APPL_USES_TRILINOS__ -D__UQ_USES_TRILINOS__
INC_PATHS = \
	-I. \
	-I../../../inc \
	-I/opt/openmpi-gnu/gnu/1.2.4/include \
	-I/h2/prudenci/Installations/Trilinos_8_0_7/include \
	-I/h2/prudenci/Installations/Boost_1_35_0/include/boost-1_35

LIB_PATHS = \
	-L../../../lib \
	-L/opt/openmpi-gnu/gnu/1.2.4/lib \
	-L/h2/prudenci/Installations/Trilinos_8_0_7/lib
	-L/h2/prudenci/Installations/Boost_1_35_0/lib

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
	-lboost_program_options-gcc41-mt-1_35

endif

CXX = g++
CXXFLAGS = -Wall -c

all:	uqHimmelEx

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_lib clean_exec

clean_lib:
	cd ../../../libs/basic/src; make clean_all UQBT=$(UQBT)
	cd ../../../libs/mcmc/src; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~
	rm -f *.o
	rm -f uqHimmelEx_$(UQBT)

uqHimmelEx: libs uqHimmelEx_$(UQBT).o
	$(CXX) $(LIB_PATHS) uqHimmelEx_$(UQBT).o -o uqHimmelEx_$(UQBT) $(LIBS)

libs:
	cd ../../../libs/basic/src; make libuq_basic_all UQBT=$(UQBT)
	cd ../../../libs/mcmc/src; make libuq_mcmc_all UQBT=$(UQBT)

.SUFFIXES: .o .C

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<
