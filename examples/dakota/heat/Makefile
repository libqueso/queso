include ../../../Make.common

INC_PATHS = \
	-I. \
	-I../inc \
	-I../../../libs/basic/inc \
	-I../../../libs/fp/inc \
	$(QUESO_TOOL_INCLUDE)

LIBS = \
	-L../../../libs/basic/lib \
	-L../../../libs/fp/lib \
	-luq_fp_any \
	-luq_basic_any \
	-luq_basic_$(UQBT) \
	$(QUESO_TOOL_LIBS)

all:	uqHeatEx

clean_lib:
	cd ../../../libs/basic/src; make clean_all UQBT=$(UQBT)
	cd ../../../libs/fp/src; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~
	rm -f *.o
	rm -f uqHeatEx_$(UQBT)

uqHeatEx: libs uqHeatEx.o 
	$(CXX) uqHeatEx.o -o uqHeatEx_$(UQBT) $(LIBS)

libs:
	cd ../../../libs/basic/src; make libuq_basic_all UQBT=$(UQBT)
	cd ../../../libs/fp/src; make libuq_fp_all UQBT=$(UQBT)

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<
