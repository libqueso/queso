include ../../../Make.common

ifeq ($(UQBT),gsl)

DEFINES = \
	-D__UQ_USES_GSL__
INC_PATHS = \
	-I../../basic/inc \
	-I../inc \
	-I$(MPI_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(TRILINOS_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(BOOST_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(GSL_INCLUDE_PATH_FOR_QUESO_TOOL)

else

DEFINES = \
	-D__UQ_USES_TRILINOS__
INC_PATHS = \
	-I../../basic/inc \
	-I../inc \
	-I$(MPI_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(TRILINOS_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(BOOST_INCLUDE_PATH_FOR_QUESO_TOOL) \
	-I$(GSL_INCLUDE_PATH_FOR_QUESO_TOOL)

endif

OBJS_ANY =

OBJS_gsl = 

OBJS_tri = uqTrilinosParamSpace.o \
	uqTrilinosFinDimLinearSpace.o

all:	libuq_queso_all

echo_num_proc:
	@echo 'UQBT = ' $(UQBT)

clean_lib: echo_num_proc clean_obj clean_docs 
	rm -f ../lib/libuq_queso_any.a
	rm -f ../lib/libuq_queso_$(UQBT).a

clean_exec:

clean_obj:
	rm -f ../inc/*~
	rm -rf ../inc/html
	rm -f *~
	rm -f *.o

clean_docs:
	rm -rf html

all_docs:
	doxygen uqLibsMcmcSrc.dox

libuq_queso_all: libuq_queso_any libuq_queso_specific

libuq_queso_any: echo_num_proc $(OBJS_ANY)
	$(AR) $(ARFLAGS) ../lib/libuq_queso_any.a $(OBJS_ANY)

libuq_queso_specific: echo_num_proc $(OBJS_$(UQBT))
	$(AR) $(ARFLAGS) ../lib/libuq_queso_$(UQBT).a $(OBJS_$(UQBT))

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<

