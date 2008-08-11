ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_heat

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_heat clean_exec

clean_heat:
	cd heat; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~

all_heat:
	cd heat; make UQBT=$(UQBT)
