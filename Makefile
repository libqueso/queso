
default: all
	@true # override later % rule

dirlist := $(wildcard */)

%:
	@echo make $@
	@for dir in $(dirlist); do \
	  if [ -f "$${dir}"/Makefile ]; then \
	    $(MAKE) -C "$${dir}" $@ || exit 1; \
	  fi; \
	done
