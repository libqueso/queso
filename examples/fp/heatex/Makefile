override UQBT = gsl

DEFINES = \
	-D__EXAMPLE_USES_GSL__
INC_PATHS = \
	-I. \
	-I../../../libs/basic/inc \
	-I../../../libs/fp/inc \
	-I../inc \
	-I$(BOOST_INCLUDE_PATH_FOR_PECOS_TOOLKIT) \
	-I$(GSL_INCLUDE_PATH_FOR_PECOS_TOOLKIT)
LIB_PATHS = \
	-L../../../libs/basic/lib \
	-L../../../libs/fp/lib \
	-L../lib \
	-L$(BOOST_LIB_PATH_FOR_PECOS_TOOLKIT) \
	-L$(GSL_LIB_PATH_FOR_PECOS_TOOLKIT)
LIBS = \
	-luq_fp_any \
	-luq_basic_$(UQBT) \
	-luq_basic_any \
	-luq_fp_any \
	-lgsl \
	-lgslcblas \
	-lboost_program_options-$(COMPILER_SPEC_IN_BOOST_LIB_NAME)-mt-1_35

CXX = g++
CXXFLAGS = -Wall -c

all:	uq1DFem

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_lib clean_exec

clean_lib:
	cd ../../../libs/basic/src; make clean_all UQBT=$(UQBT)
	cd ../../../libs/fp/src; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~
	rm -f *.o
	rm -f uq1DFem_$(UQBT)

uq1DFem: libs uq1DFem.o 
	$(CXX) $(LIB_PATHS) uq1DFem.o -o uq1DFem_$(UQBT) $(LIBS)

libs:
	cd ../../../libs/basic/src; make libuq_basic_all UQBT=$(UQBT)
	cd ../../../libs/fp/src; make libuq_fp_all UQBT=$(UQBT)

.SUFFIXES: .o .C

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<
