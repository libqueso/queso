ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_ex4

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_ex4

clean_ex4:
	cd ex4; make clean_all UQBT=$(UQBT)

all_ex4:
	cd ex4; make UQBT=$(UQBT)

