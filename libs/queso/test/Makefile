ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

ifeq ($(BOOST_INCLUDE_PATH_FOR_QUESO_TOOL),)
override BOOST_INCLUDE_PATH_FOR_QUESO_TOOL = .
endif
ifeq ($(GSL_INCLUDE_PATH_FOR_QUESO_TOOL),)
override GSL_INCLUDE_PATH_FOR_QUESO_TOOL = .
endif

ifeq ($(UQBT),gsl)

DEFINES = \
	-D__UQ_USES_GSL__

else

DEFINES = \
	-D__UQ_USES_TRILINOS__

endif

INC_PATHS = \
	-I../../basic/inc \
	-I../inc \
	-I$(MPI_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(TRILINOS_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(BOOST_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(GSL_INCLUDE_PATH_FOR_QUESO_TOOL)

LIB_PATHS = \
	-L../../basic/lib \
	-L$(MPI_LIB_PATH_FOR_QUESO_TOOL) \
	-L$(TRILINOS_LIB_PATH_FOR_QUESO_TOOL) \
	-L$(BOOST_LIB_PATH_FOR_QUESO_TOOL) \
	-L$(GSL_LIB_PATH_FOR_QUESO_TOOL)

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
	-lboost_program_options-$(COMPILER_SPEC_IN_BOOST_LIB_NAME)-mt-1_35 \
	-l$(BOOST_UNIT_TEST_LIB_FOR_QUESO_TOOL)


AR = ar
ARFLAGS = -rv
CXX = g++
CXXFLAGS = -Wall -c -g

OBJS_ANY =

OBJS_gsl = uqGslVectorSpace.o

OBJS_tri = uqTrilinosParamSpace.o \
	uqTrilinosFinDimLinearSpace.o

UNIT_TESTS = uqTestExample \
	uqTestVectorPdf_gsl

all: check

echo_num_proc:
	@echo 'UQBT = ' $(UQBT)

clean: clean_all

clean_all: echo_num_proc clean_obj clean_docs 
	rm -f ../lib/libuq_queso_any.a
	rm -f ../lib/libuq_queso_$(UQBT).a
	rm -f $(UNIT_TESTS)

clean_obj:
	rm -f ../inc/*~
	rm -rf ../inc/html
	rm -f *~
	rm -f *.o

clean_docs:
	rm -rf html

.SUFFIXES: .o .C

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<

# Unit testing
#
# For new unit test, add executable name to UNIT_TESTS variable above, and
# create target and rule for building exectuable.  Use 'make check' to
# build and run all unit tests.

check: $(UNIT_TESTS)
	@for mytest in $^; do echo "Running ./$${mytest}..."; ./$${mytest}; done

uqTestExample: uqTestExample.o
	$(CXX) -L$(BOOST_LIB_PATH_FOR_QUESO_TOOL) $< -o $@ -l$(BOOST_UNIT_TEST_LIB_FOR_QUESO_TOOL)


uqTestVectorPdf_gsl: uqTestVectorPdf_gsl.o libs_basic
	$(CXX) $(LIB_PATHS) $< -o $@ $(LIBS)

uqTestVectorPdf_gsl.o: uqTestVectorPdf_gsl.C ../inc/uqVectorPdf.h
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<

libs_basic:
	cd ../../basic/src; make libuq_basic_all UQBT=$(UQBT)