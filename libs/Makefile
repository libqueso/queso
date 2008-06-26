ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_basic all_mcmc all_fp

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_basic clean_mcmc clean_fp clean_exec

clean_basic:
	cd basic/src; make clean_all UQBT=$(UQBT)

clean_mcmc:
	cd mcmc/src; make clean_all UQBT=$(UQBT)

clean_fp:
	cd fp/src; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~

all_basic:
	cd basic/src; make UQBT=$(UQBT)

all_mcmc:
	cd mcmc/src; make UQBT=$(UQBT)

all_fp:
	cd fp/src; make UQBT=$(UQBT)
