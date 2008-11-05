ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_queso all_dakota

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_queso clean_dakota clean_docs clean_exec

clean_queso:
	cd queso; make clean_all UQBT=$(UQBT)

clean_dakota:
	cd dakota; make clean_all UQBT=$(UQBT)

clean_docs:
	rm -rf html

clean_exec:
	rm -f *~

all_docs:
	doxygen uqAppls.dox
	cd queso; make all_docs
	cd dakota; make all_docs

all_queso:
	cd queso; make UQBT=$(UQBT)

all_dakota:
	cd dakota; make UQBT=$(UQBT)
