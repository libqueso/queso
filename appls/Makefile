ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_mcmc all_fp

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_mcmc clean_fp clean_exec

clean_mcmc:
	cd mcmc; make clean_all UQBT=$(UQBT)

clean_fp:
	cd fp; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~

all_mcmc:
	cd mcmc; make UQBT=$(UQBT)

all_fp:
	cd fp; make UQBT=$(UQBT)
