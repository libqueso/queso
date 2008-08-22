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
ifeq ($(MPI_INCLUDE_PATH_FOR_MCMC_TOOL),)
override MPI_INCLUDE_PATH_FOR_MCMC_TOOL = .
endif
ifeq ($(TRILINOS_INCLUDE_PATH_FOR_MCMC_TOOL),)
override TRILINOS_INCLUDE_PATH_FOR_MCMC_TOOL = .
endif

ifeq ($(UQBT),gsl)

DEFINES = \
	-D__UQ_USES_GSL__
INC_PATHS = \
	-I../inc \
	-I$(MPI_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(TRILINOS_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(BOOST_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(GSL_INCLUDE_PATH_FOR_MCMC_TOOL)

else

DEFINES = \
	-D__UQ_USES_TRILINOS__
INC_PATHS = \
	-I../inc \
	-I$(MPI_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(TRILINOS_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(BOOST_INCLUDE_PATH_FOR_MCMC_TOOL) \
	-I$(GSL_INCLUDE_PATH_FOR_MCMC_TOOL)

endif

AR = ar
ARFLAGS = -rv
CXX = g++
CXXFLAGS = -Wall -c -g

OBJS_ANY = uqEnvironment.o \
	uqDefines.o \
	uqVector.o \
	uqMatrix.o \
	uqMiscellaneous.o

OBJS_gsl = uqGslVector.o \
	uqGslMatrix.o \
	uqGslOdeSolver.o

OBJS_tri = uqTrilinosVector.o \
	uqTrilinosMatrix.o

all:	libuq_basic_all

all_docs:
	doxygen uqLibsBasicSrc.dox

echo_num_proc:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_num_proc clean_obj
	rm -f ../lib/libuq_basic_any.a
	rm -f ../lib/libuq_basic_$(UQBT).a

clean_obj:
	rm -f ../inc/*~
	rm -f *~
	rm -f *.o

libuq_basic_all: libuq_basic_any libuq_basic_specific

libuq_basic_any: echo_num_proc $(OBJS_ANY)
	$(AR) $(ARFLAGS) ../lib/libuq_basic_any.a $(OBJS_ANY)

libuq_basic_specific: echo_num_proc $(OBJS_$(UQBT))
	$(AR) $(ARFLAGS) ../lib/libuq_basic_$(UQBT).a $(OBJS_$(UQBT))

.SUFFIXES: .o .C

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<
