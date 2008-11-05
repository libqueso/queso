ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_basic all_queso all_fp

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_basic clean_queso clean_fp clean_docs clean_exec

clean_basic:
	cd basic; make clean_all UQBT=$(UQBT)

clean_queso:
	cd queso; make clean_all UQBT=$(UQBT)

clean_fp:
	cd fp; make clean_all UQBT=$(UQBT)

clean_docs:
	rm -rf html

clean_exec:
	rm -f *~

all_docs:
	doxygen uqLibs.dox
	cd basic; make all_docs
	cd queso; make all_docs
	cd fp; make all_docs

all_basic:
	cd basic; make UQBT=$(UQBT)

all_queso:
	cd queso; make UQBT=$(UQBT)

all_fp:
	cd fp; make UQBT=$(UQBT)
