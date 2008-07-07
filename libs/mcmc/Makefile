ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_mcmc

echo_mcmc_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_mcmc_type clean_mcmc clean_docs clean_exec

clean_mcmc:
	cd src; make clean_all UQBT=$(UQBT)

clean_docs:
	rm -rf html

clean_exec:
	rm -f *~

all_docs:
	doxygen uqLibsMcmc.dox
	cd src; make all_docs
	cd inc; make all_docs

all_mcmc:
	cd src; make UQBT=$(UQBT)
