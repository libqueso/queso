
default: all
	@true # override later % rule

dirlist := $(wildcard */)

clean:
	@echo make clean
	@for dir in $(dirlist); do \
	  if [ -f "$${dir}"/Makefile ]; then \
	    if [ -f "$${dir}"/bootstrap ]; then \
	      if [ -f "$${dir}"/maintainer_clean_failed ]; then \
                echo "Skipping broken maintainer-clean in $${dir}"; \
              else \
                echo "Running maintainer-clean in $${dir}"; \
	        $(MAKE) -C "$${dir}" maintainer-clean || touch "$${dir}"/maintainer_clean_failed; \
	      fi; \
            else \
	      if [ -f "$${dir}"/clean_failed ]; then \
                echo "Skipping broken clean in $${dir}"; \
              else \
                echo "Running clean in $${dir}"; \
	        $(MAKE) -C "$${dir}" clean || touch "$${dir}"/clean_failed; \
	      fi; \
	    fi; \
	  fi; \
	done

%:
	@echo make $@
	@for dir in $(dirlist); do \
	  if [ -f "$${dir}"/bootstrap ]; then \
	    if [ -f "$${dir}"/bootstrap_failed ]; then \
              echo "Skipping broken bootstrap in $${dir}"; \
            else \
              echo "Running bootstrap in $${dir}"; \
	      (cd "$${dir}" && ./bootstrap || touch bootstrap_failed); \
	    fi; \
	  fi; \
	  if [ -f "$${dir}"/configure ]; then \
	    if [ -f "$${dir}"/configure_failed ]; then \
              echo "Skipping broken configure in $${dir}"; \
            else \
               echo "Running configure in $${dir}"; \
	       (cd "$${dir}" && ./configure || touch configure_failed); \
	    fi; \
	  fi; \
	  if [ -f "$${dir}"/Makefile ]; then \
	    if [ -f "$${dir}"/configure_failed ]; then \
              echo "Skipping broken make $@ in $${dir}"; \
            else \
              echo "Running make $@ in $${dir}"; \
	      $(MAKE) -C "$${dir}" $@ || touch "$${dir}"/make_$@_failed; \
	    fi; \
	  fi; \
	done
