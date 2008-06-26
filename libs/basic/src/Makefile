ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

ifeq ($(UQBT),gsl)

DEFINES = \
	-D__UQ_USES_GSL__
INC_PATHS = \
	-I../inc \
	-I/h2/prudenci/Installations/Boost_1_35_0/include/boost-1_35

else

DEFINES = \
	-D__UQ_USES_TRILINOS__
INC_PATHS = \
	-I../inc \
	-I/h2/prudenci/Installations/Boost_1_35_0/include/boost-1_35 \
	-I/opt/openmpi-gnu/gnu/1.2.4/include \
	-I/h2/prudenci/Installations/Trilinos_8_0_7/include

endif

AR = ar
ARFLAGS = -rv
CXX = g++
CXXFLAGS = -Wall -c

OBJS_ANY = uqEnvironment.o \
	uqVector.o \
	uqMatrix.o \
	uqMiscellaneous.o

OBJS_gsl = uqGslVector.o \
	uqGslMatrix.o \
	uqGslOdeSolver.o

OBJS_tri = uqTrilinosVector.o \
	uqTrilinosMatrix.o \
	uqTrilinosParamSpace.o

all:	libuq_basic_all

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
