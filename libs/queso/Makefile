ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_queso

echo_queso_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_queso_type clean_queso clean_docs clean_exec

clean_queso:
	cd src; make clean_all UQBT=$(UQBT)
	cd test; make clean_all UQBT=$(UQBT)

clean_docs:
	rm -rf html

clean_exec:
	rm -f *~

all_docs:
	doxygen uqLibsQueso.dox
	cd src; make all_docs
	cd inc; make all_docs

all_queso:
	cd src; make UQBT=$(UQBT)

check:
	cd test; make check UQBT=$(UQBT)