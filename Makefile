ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_libs all_examples

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_libs clean_examples clean_docs clean_exec

clean_libs:
	cd libs; make clean_all UQBT=$(UQBT)

clean_examples:
	cd examples; make clean_all UQBT=$(UQBT)

clean_docs:
	rm -rf html
	cd libs; make clean_docs
	cd examples; make clean_docs

clean_exec:
	rm -f *~

all_docs:
	doxygen uq.dox
	cd libs; make all_docs
#	cd examples; make all_docs

all_libs:
	cd libs; make UQBT=$(UQBT)

all_examples:
	cd examples; make UQBT=$(UQBT)
