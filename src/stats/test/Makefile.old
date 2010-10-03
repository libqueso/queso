include ../../../Make.common

INC_PATHS = \
	-I../../basic/inc \
	-I../inc \
	$(QUESO_TOOL_INCLUDE)

LIBS = \
	-L../../basic/lib \
	-luq_basic_$(UQBT) \
	-luq_basic_any \
	$(QUESO_TOOL_LIBS) \
	-l$(BOOST_UNIT_TEST_LIB_FOR_QUESO_TOOL)


OBJS_ANY =

OBJS_gsl = uqGslVectorSpace.o

OBJS_tri = uqTrilinosParamSpace.o \
	uqTrilinosFinDimLinearSpace.o

UNIT_TESTS = uqTestExample \
	uqTestVectorPdf_gsl \
	uqTestVectorRealizer_gsl \
	uqTestVectorRV_gsl

all: check

echo_num_proc:
	@echo 'UQBT = ' $(UQBT)

clean: clean_all

clean_lib: echo_num_proc clean_obj clean_docs 
	rm -f ../lib/libuq_queso_any.a
	rm -f ../lib/libuq_queso_$(UQBT).a
	rm -f $(UNIT_TESTS)

clean_obj:
	rm -f ../inc/*~
	rm -rf ../inc/html
	rm -f *~
	rm -f *.o

clean_exec:

clean_docs:
	rm -rf html

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
	$(CXX) $< -o $@ $(LIBS)

uqTestVectorPdf_gsl.o: uqTestVectorPdf_gsl.C ../inc/uqVectorPdf.h
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<

uqTestVectorRealizer_gsl: uqTestVectorRealizer_gsl.o libs_basic
	$(CXX) $< -o $@ $(LIBS)

uqTestVectorRealizer_gsl.o: uqTestVectorRealizer_gsl.C ../inc/uqVectorRealizer.h
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<

uqTestVectorRV_gsl: uqTestVectorRV_gsl.o libs_basic
	$(CXX) $< -o $@ $(LIBS)

uqTestVectorRV_gsl.o: uqTestVectorRV_gsl.C ../inc/uqVectorRV.h
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<

libs_basic:
	cd ../../basic/src; make libuq_basic_all UQBT=$(UQBT)
