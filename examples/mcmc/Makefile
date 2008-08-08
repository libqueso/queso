ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_normal all_chem all_algae all_etc

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_normal clean_chem clean_algae clean_etc clean_exec

clean_normal:
	cd normal; make clean_all UQBT=$(UQBT)

clean_chem:
	cd chemicalReactions; make clean_all UQBT=$(UQBT)

clean_algae:
	cd algae; make clean_all UQBT=$(UQBT)

clean_etc:
	cd etc; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~

all_normal:
	cd normal; make UQBT=$(UQBT)

all_chem:
	cd chemicalReactions; make UQBT=$(UQBT)

all_algae:
	cd algae; make UQBT=$(UQBT)

all_etc:
	cd etc; make UQBT=$(UQBT)
