ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_src

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_src

clean_src:
	cd src; make clean_all UQBT=$(UQBT)

all_src:
	cd src; make UQBT=$(UQBT)

