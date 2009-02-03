ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_queso

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_queso clean_docs clean_exec

clean_queso:
	cd queso; make clean_all UQBT=$(UQBT)

clean_docs:
	rm -rf html

clean_exec:
	rm -f *~

all_docs:
	doxygen uqAppls.dox
	cd queso; make all_docs

all_queso:
	cd queso; make UQBT=$(UQBT)
