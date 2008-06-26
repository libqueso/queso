ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_normal all_himmel all_algae all_template

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_normal clean_himmel clean_algae clean_template clean_exec

clean_normal:
	cd normalex; make clean_all UQBT=$(UQBT)

clean_himmel:
	cd himmelblauex; make clean_all UQBT=$(UQBT)

clean_algae:
	cd algaeex; make clean_all UQBT=$(UQBT)

clean_template:
	cd templateex; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~

all_normal:
	cd normalex; make UQBT=$(UQBT)

all_himmel:
	cd himmelblauex; make UQBT=$(UQBT)

all_algae:
	cd algaeex; make UQBT=$(UQBT)

all_template:
	cd templateex; make UQBT=$(UQBT)
