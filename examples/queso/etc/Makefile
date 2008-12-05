include ../../../Make.common

INC_PATHS = \
	-I. \
	-I../inc \
	-I../../../libs/basic/inc \
	-I../../../libs/queso/inc \
	$(QUESO_TOOL_INCLUDE)

LIBS = \
	-L../../../libs/basic/lib \
	-L../../../libs/queso/lib \
	-luq_basic_$(UQBT) \
	-luq_basic_any \
	-luq_queso_$(UQBT) \
	-luq_queso_any \
	$(QUESO_TOOL_LIBS)

all:	uqEtcEx

clean_lib:
	cd ../../../libs/basic/src; make clean_all UQBT=$(UQBT)
	cd ../../../libs/queso/src; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~
	rm -f *.o
	rm -f uqEtcEx_$(UQBT)

uqEtcEx: libs uqEtcEx_$(UQBT).o
	$(CXX) uqEtcEx_$(UQBT).o -o uqEtcEx_$(UQBT) $(LIBS)

libs:
	cd ../../../libs/basic/src; make libuq_basic_all UQBT=$(UQBT)
	cd ../../../libs/queso/src; make libuq_queso_all UQBT=$(UQBT)

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<
