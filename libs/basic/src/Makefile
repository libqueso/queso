include ../../../Make.common

ifeq ($(UQBT),gsl)

DEFINES = \
	-D__UQ_USES_GSL__
INC_PATHS = \
	-I../inc \
	-I$(MPI_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(TRILINOS_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(BOOST_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(GSL_INCLUDE_PATH_FOR_QUESO_TOOL)

else

DEFINES = \
	-D__UQ_USES_TRILINOS__
INC_PATHS = \
	-I../inc \
	-I$(MPI_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(TRILINOS_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(BOOST_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(GSL_INCLUDE_PATH_FOR_QUESO_TOOL)

endif

AR = ar
ARFLAGS = -rv
CXX = g++
CXXFLAGS += -Wall -c -g

OBJS_ANY = uqEnvironment.o \
	uqDefines.o \
	uqVector.o \
	uqMatrix.o \
	uqMiscellaneous.o \
	uqRealFft.o \
	uqComplexFft.o \
	uqChainStatisticalOptions.o

OBJS_gsl = uqGslVector.o \
	uqGslMatrix.o \
	uqGslOdeSolver.o \
	uqGslVectorSpace.o

OBJS_tri = uqTrilinosVector.o \
	uqTrilinosMatrix.o

all:	libuq_basic_all

all_docs:
	doxygen uqLibsBasicSrc.dox

echo_num_proc:
	@echo 'UQBT = ' $(UQBT)

clean_lib: echo_num_proc clean_obj
	rm -f ../lib/libuq_basic_any.a
	rm -f ../lib/libuq_basic_$(UQBT).a

clean_exec:

clean_obj:
	rm -f ../inc/*~
	rm -f *~
	rm -f *.o

libuq_basic_all: libuq_basic_any libuq_basic_specific

libuq_basic_any: echo_num_proc $(OBJS_ANY)
	$(AR) $(ARFLAGS) ../lib/libuq_basic_any.a $(OBJS_ANY)

libuq_basic_specific: echo_num_proc $(OBJS_$(UQBT))
	$(AR) $(ARFLAGS) ../lib/libuq_basic_$(UQBT).a $(OBJS_$(UQBT))

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<
