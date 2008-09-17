ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_ex3 all_ex4

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_ex3 clean_ex4

clean_ex3:
	cd ex3; make clean_all UQBT=$(UQBT)

clean_ex4:
	cd ex4; make clean_all UQBT=$(UQBT)

all_ex3:
	cd ex3; make UQBT=$(UQBT)

all_ex4:
	cd ex4; make UQBT=$(UQBT)

