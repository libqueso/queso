override UQBT = gsl

INC_PATHS = \
	-I../inc \
	-I../../basic/inc

AR = ar
ARFLAGS = -rv
CXX = g++
CXXFLAGS = -Wall -c

OBJS_ANY = uq1DScalarElement.o \
	uq1DNode.o \
	uq1DLocalDof.o \
	uq1DProblem.o \
	uq1DFemDakotaHandler.o

all:	libuq_fp_all

echo_num_proc:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_num_proc clean_obj
	rm -f ../lib/libuq_fp_any.a

clean_obj:
	rm -f ../inc/*~
	rm -f *~
	rm -f *.o

libuq_fp_all: libuq_fp_any libuq_fp_specific

libuq_fp_any: echo_num_proc $(OBJS_ANY)
	$(AR) $(ARFLAGS) ../lib/libuq_fp_any.a $(OBJS_ANY)

libuq_fp_specific: echo_num_proc $(OBJS_$(UQBT))

.SUFFIXES: .o .C

%.o: %.C
	$(CXX) $(INC_PATHS) $(CXXFLAGS) $<
