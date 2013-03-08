
default: all
	@true # override later % rule

dirlist := $(wildcard */)

clean:
	@echo make clean
	@for dir in $(dirlist); do \
	  if [ -f "$${dir}"/Makefile ]; then \
	    $(MAKE) -C "$${dir}" clean || exit 1; \
	  fi; \
	done

%:
	@echo make $@
	@for dir in $(dirlist); do \
	  if [ -f "$${dir}"/bootstrap ]; then \
            echo "Running bootstrap in $${dir}"; \
	    (cd "$${dir}" && ./bootstrap) || exit 1; \
	  fi; \
	  if [ -f "$${dir}"/configure ]; then \
            echo "Running configure in $${dir}"; \
	    (cd "$${dir}" && ./configure) || exit 1; \
	  fi; \
	  if [ -f "$${dir}"/Makefile ]; then \
	    $(MAKE) -C "$${dir}" $@ || exit 1; \
	  fi; \
	done
