include ../../../../Make.common

INC_PATHS = \
	-I. \
	-I../../inc \
	-I../../../../libs/basic/inc \
	-I../../../../libs/queso/inc \
	-I$(BURGERS_INCLUDE_PATH_FOR_QUESO_TOOL) \
	$(QUESO_TOOL_INCLUDE)

LIBS = \
	-L../../../../libs/basic/lib \
	-L../../../../libs/queso/lib \
	-luq_queso_$(UQBT) \
	-luq_queso_any \
	-luq_basic_$(UQBT) \
	-luq_basic_any \
	-L$(BURGERS_LIB_PATH_FOR_QUESO_TOOL) \
	-llegendrePoly \
	$(QUESO_TOOL_LIBS)

all:	Burgers_Cal_Variance

clean_lib:
	cd ../../../../libs/basic/src; make clean_all UQBT=$(UQBT)
	cd ../../../../libs/queso/src; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~
	rm -f *.o
	rm -f uqBurgers_Cal_Variance_$(UQBT)

Burgers_Cal_Variance: libs uqBurgers_Cal_Variance_$(UQBT).o
	$(CXX) uqBurgers_Cal_Variance_$(UQBT).o -o uqBurgers_Cal_Variance_$(UQBT) $(LIBS)

libs:
	cd ../../../../libs/basic/src; make libuq_basic_all UQBT=$(UQBT)
	cd ../../../../libs/queso/src; make libuq_queso_all UQBT=$(UQBT)

%.o: %.C
	$(CXX) $(DEFINES) $(INC_PATHS) $(CXXFLAGS) $<
