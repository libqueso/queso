ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_fp

echo_fp_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_fp_type clean_fp clean_docs clean_exec

clean_fp:
	cd src; make clean_all UQBT=$(UQBT)

clean_docs:
	rm -rf html

clean_exec:
	rm -f *~

all_docs:
	doxygen uqLibsFp.dox
	cd src; make all_docs
	cd inc; make all_docs

all_fp:
	cd src; make UQBT=$(UQBT)
