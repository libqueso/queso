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

ifeq ($(UQBT),gsl)

DEFINES = \
	-D__UQ_USES_GSL__
INC_PATHS = \
	-I../../basic/inc \
	-I../inc \
	-I$(BOOST_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(GSL_INCLUDE_PATH_FOR_MCMC_TOOL)

else

DEFINES = \
	-D__UQ_USES_TRILINOS__
INC_PATHS = \
	-I../../basic/inc \
	-I../inc \
	-I$(MPI_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(TRILINOS_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(BOOST_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(GSL_INCLUDE_PATH_FOR_MCMC_TOOL)

endif

AR = ar
ARFLAGS = -rv
CXX = g++
CXXFLAGS = -Wall -c

OBJS_ANY = uqParameter.o \
	uqObservable.o

OBJS_gsl = uqGslParamSpace.o \
	uqGslObservableSpace.o \
	uqGslFinDimLinearSpace.o

OBJS_tri = uqTrilinosParamSpace.o \
	uqTrilinosObservableSpace.o \
	uqTrilinosFinDimLinearSpace.o

all:	libuq_mcmc_all

echo_num_proc:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_num_proc clean_obj
	rm -f ../lib/libuq_mcmc_any.a
	rm -f ../lib/libuq_mcmc_$(UQBT).a

clean_obj:
	rm -f ../inc/*~
	rm -rf ../inc/html
	rm -f *~
	rm -f *.o

all_docs:
	doxygen uqLibsMcmcSrc.dox

libuq_mcmc_all: libuq_mcmc_any libuq_mcmc_specific

libuq_mcmc_any: echo_num_proc $(OBJS_ANY)
	$(AR) $(ARFLAGS) ../lib/libuq_mcmc_any.a $(OBJS_ANY)

libuq_mcmc_specific: echo_num_proc $(OBJS_$(UQBT))
	$(AR) $(ARFLAGS) ../lib/libuq_mcmc_$(UQBT).a $(OBJS_$(UQBT))

.SUFFIXES: .o .C

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<
