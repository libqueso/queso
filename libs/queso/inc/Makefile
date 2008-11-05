ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

all:	all_docs

echo_mcmc_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_mcmc_type clean_docs clean_exec

clean_docs:
	rm -rf html

clean_exec:
	rm -f *~

all_docs:
	doxygen uqLibsMcmcInc.dox
