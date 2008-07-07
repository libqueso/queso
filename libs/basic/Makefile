ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_basic

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_basic clean_docs clean_exec

clean_basic:
	cd src; make clean_all UQBT=$(UQBT)

clean_docs:
	rm -rf html

clean_exec:
	rm -f *~

all_docs:
	doxygen uqLibsBasic.dox
	cd src; make all_docs
	cd inc; make all_docs

all_basic:
	cd src; make UQBT=$(UQBT)
