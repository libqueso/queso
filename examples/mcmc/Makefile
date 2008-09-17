ifeq ($(UQBT),)
override UQBT = gsl
else ifeq ($(UQBT),gsl)
else
override UQBT = tri
endif

override LOWERCASE_EXAMPLE_NAME=$(n1)
override MIXED_EXAMPLE_NAME=$(n2)
override CAPITAL_EXAMPLE_NAME=$(n3)

#all:	all_normal all_chem all_algae all_tga all_etc
all:	all_tga

echo_basic_type:
	@echo 'UQBT = ' $(UQBT)

clean_all: echo_basic_type clean_normal clean_chem clean_algae clean_tga clean_etc clean_exec

clean_normal:
	cd normal; make clean_all UQBT=$(UQBT)

clean_chem:
	cd chemicalReactions; make clean_all UQBT=$(UQBT)

clean_algae:
	cd algae; make clean_all UQBT=$(UQBT)

clean_tga:
	cd tga; make clean_all UQBT=$(UQBT)

clean_etc:
	cd etc; make clean_all UQBT=$(UQBT)

clean_exec:
	rm -f *~

all_normal:
	cd normal; make UQBT=$(UQBT)

all_chem:
	cd chemicalReactions; make UQBT=$(UQBT)

all_algae:
	cd algae; make UQBT=$(UQBT)

all_tga:
	cd tga; make UQBT=$(UQBT)

all_etc:
	cd etc; make UQBT=$(UQBT)

create_example:
	@echo 'NAME OF EXAMPLE =' $(LOWERCASE_EXAMPLE_NAME)
	mkdir $(LOWERCASE_EXAMPLE_NAME)
	cp etc/Makefile      $(LOWERCASE_EXAMPLE_NAME)
	cp etc/uqEtcEx.h     $(LOWERCASE_EXAMPLE_NAME)/uq$(MIXED_EXAMPLE_NAME)Ex.h
	cp etc/uqEtcEx_gsl.C $(LOWERCASE_EXAMPLE_NAME)/uq$(MIXED_EXAMPLE_NAME)Ex_gsl.C
	cp etc/uqEtcEx_tri.C $(LOWERCASE_EXAMPLE_NAME)/uq$(MIXED_EXAMPLE_NAME)Ex_tri.C
	cp etc/uqEtcEx.inp   $(LOWERCASE_EXAMPLE_NAME)/uq$(MIXED_EXAMPLE_NAME)Ex.inp
	cp etc/calib.par     $(LOWERCASE_EXAMPLE_NAME)/calib.par
	cp etc/calib.obs     $(LOWERCASE_EXAMPLE_NAME)/calib.obs
	cd $(LOWERCASE_EXAMPLE_NAME); sed 's/uqEtc/uq$(MIXED_EXAMPLE_NAME)/g' Makefile > Makefile2 && mv Makefile2 Makefile
	cd $(LOWERCASE_EXAMPLE_NAME); sed 's/uqEtc/uq$(MIXED_EXAMPLE_NAME)/g' uq$(MIXED_EXAMPLE_NAME)Ex.h > uq$(MIXED_EXAMPLE_NAME)Ex2.h && mv uq$(MIXED_EXAMPLE_NAME)Ex2.h uq$(MIXED_EXAMPLE_NAME)Ex.h
	cd $(LOWERCASE_EXAMPLE_NAME); sed 's/etc/$(LOWERCASE_EXAMPLE_NAME)/g' uq$(MIXED_EXAMPLE_NAME)Ex.h > uq$(MIXED_EXAMPLE_NAME)Ex2.h && mv uq$(MIXED_EXAMPLE_NAME)Ex2.h uq$(MIXED_EXAMPLE_NAME)Ex.h
	cd $(LOWERCASE_EXAMPLE_NAME); sed 's/ETC/$(CAPITAL_EXAMPLE_NAME)/g'   uq$(MIXED_EXAMPLE_NAME)Ex.h > uq$(MIXED_EXAMPLE_NAME)Ex2.h && mv uq$(MIXED_EXAMPLE_NAME)Ex2.h uq$(MIXED_EXAMPLE_NAME)Ex.h
	cd $(LOWERCASE_EXAMPLE_NAME); sed 's/etc/$(LOWER_EXAMPLE_NAME)/g'     uq$(MIXED_EXAMPLE_NAME)Ex_gsl.C > uq$(MIXED_EXAMPLE_NAME)Ex2_gsl.C && mv uq$(MIXED_EXAMPLE_NAME)Ex2_gsl.C uq$(MIXED_EXAMPLE_NAME)Ex_gsl.C
	cd $(LOWERCASE_EXAMPLE_NAME); sed 's/Etc/$(MIXED_EXAMPLE_NAME)/g'     uq$(MIXED_EXAMPLE_NAME)Ex_gsl.C > uq$(MIXED_EXAMPLE_NAME)Ex2_gsl.C && mv uq$(MIXED_EXAMPLE_NAME)Ex2_gsl.C uq$(MIXED_EXAMPLE_NAME)Ex_gsl.C
	cd $(LOWERCASE_EXAMPLE_NAME); sed 's/etc/$(LOWER_EXAMPLE_NAME)/g'     uq$(MIXED_EXAMPLE_NAME)Ex_tri.C > uq$(MIXED_EXAMPLE_NAME)Ex2_tri.C && mv uq$(MIXED_EXAMPLE_NAME)Ex2_tri.C uq$(MIXED_EXAMPLE_NAME)Ex_tri.C
	cd $(LOWERCASE_EXAMPLE_NAME); sed 's/Etc/$(MIXED_EXAMPLE_NAME)/g'     uq$(MIXED_EXAMPLE_NAME)Ex_tri.C > uq$(MIXED_EXAMPLE_NAME)Ex2_tri.C && mv uq$(MIXED_EXAMPLE_NAME)Ex2_tri.C uq$(MIXED_EXAMPLE_NAME)Ex_tri.C
	cd $(LOWERCASE_EXAMPLE_NAME); sed 's/Etc/$(MIXED_EXAMPLE_NAME)/g'     uq$(MIXED_EXAMPLE_NAME)Ex.inp   > uq$(MIXED_EXAMPLE_NAME)Ex2.inp && mv uq$(MIXED_EXAMPLE_NAME)Ex2.inp uq$(MIXED_EXAMPLE_NAME)Ex.inp
