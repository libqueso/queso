include ../../../../Make.common

INC_PATHS = \
	-I. \
	-I../../inc \
	-I../../../../libs/basic/inc \
	-I../../../../libs/queso/inc \
	$(QUESO_TOOL_INCLUDE)

LIBS = \
	-L../../../../libs/basic/lib \
	-L../../../../libs/queso/lib \
	-luq_queso_$(UQBT) \
	-luq_queso_any \
	-luq_basic_$(UQBT) \
	-luq_basic_any \
	$(QUESO_TOOL_LIBS)

all:	uqTgaEx

clean_lib:
	cd ../../../../libs/basic/src; make clean_all UQBT=$(UQBT)
	cd ../../../../libs/queso/src; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~
	rm -f *.o
	rm -f uqTgaEx_$(UQBT)

uqTgaEx: libs uqTgaEx4_$(UQBT).o
	$(CXX) uqTgaEx4_$(UQBT).o -o uqTgaEx4_$(UQBT) $(LIBS)

libs:
	cd ../../../../libs/basic/src; make libuq_basic_all UQBT=$(UQBT)
	cd ../../../../libs/queso/src; make libuq_queso_all UQBT=$(UQBT)

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<
