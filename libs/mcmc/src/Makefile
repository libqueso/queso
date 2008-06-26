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
	-I../../basic/inc \
	-I../inc \
	-I/h2/prudenci/Installations/Boost_1_35_0/include/boost-1_35

else

DEFINES = \
	-D__UQ_USES_TRILINOS__
INC_PATHS = \
	-I../../basic/inc \
	-I../inc \
	-I/h2/prudenci/Installations/Boost_1_35_0/include/boost-1_35 \
	-I/opt/openmpi-gnu/gnu/1.2.4/include \
	-I/h2/prudenci/Installations/Trilinos_8_0_7/include

endif

AR = ar
ARFLAGS = -rv
CXX = g++
CXXFLAGS = -Wall -c

OBJS_ANY = uqParameter.o

OBJS_gsl = uqGslParamSpace.o \
	uqGslFinDimLinearSpace.o

OBJS_tri = uqTrilinosVector.o \
	uqTrilinosMatrix.o \
	uqTrilinosParamSpace.o

all:	libuq_mcmc_all

echo_num_proc:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_num_proc clean_obj
	rm -f ../lib/libuq_mcmc_any.a
	rm -f ../lib/libuq_mcmc_$(UQBT).a

clean_obj:
	rm -f ../inc/*~
	rm -f *~
	rm -f *.o

libuq_mcmc_all: libuq_mcmc_any libuq_mcmc_specific

libuq_mcmc_any: echo_num_proc $(OBJS_ANY)
	$(AR) $(ARFLAGS) ../lib/libuq_mcmc_any.a $(OBJS_ANY)

libuq_mcmc_specific: echo_num_proc $(OBJS_$(UQBT))
	$(AR) $(ARFLAGS) ../lib/libuq_mcmc_$(UQBT).a $(OBJS_$(UQBT))

.SUFFIXES: .o .C

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<
