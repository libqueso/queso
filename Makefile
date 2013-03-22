
default: all
	@true # override later % rule

dirlist := $(wildcard */)

clean:
	@echo make clean
	@for dir in $(dirlist); do \
	  if [ -f "$${dir}"/Makefile ]; then \
	    if [ -f "$${dir}"/bootstrap ]; then \
	      $(MAKE) -C "$${dir}" maintainer-clean || touch "$${dir}"/maintainer_clean_failed; \
            else \
	      $(MAKE) -C "$${dir}" clean || touch "$${dir}"/clean_failed; \
	    fi; \
	  fi; \
	done

%:
	@echo make $@
	@for dir in $(dirlist); do \
	  if [ -f "$${dir}"/bootstrap ]; then \
            echo "Running bootstrap in $${dir}"; \
	    (cd "$${dir}" && ./bootstrap || touch bootstrap_failed); \
	  fi; \
	  if [ -f "$${dir}"/configure ]; then \
            echo "Running configure in $${dir}"; \
	    (cd "$${dir}" && ./configure || touch configure_failed); \
	  fi; \
	  if [ -f "$${dir}"/Makefile ]; then \
	    $(MAKE) -C "$${dir}" $@ || touch "$${dir}"/make_$@_failed; \
	  fi; \
	done
