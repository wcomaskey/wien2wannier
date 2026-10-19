### wien2wannier/Makefile
###
###    wien2wannier main Makefile
###
### Copyright 2013-2016 Elias Assmann

-include make.sys

SHELL=/bin/bash

version = v3.0

VERSION = $(shell git describe 2>/dev/null || echo $(version))

SIMPLE      := SRC_trig doc test
REALCOMPLEX := SRC_w2w SRC_wplot

SUBDIRS := $(SIMPLE) $(REALCOMPLEX)

.PHONY: all clean distclean $(SUBDIRS) install install-cli

all: SRC_w2w SRC_wplot SRC_trig

$(REALCOMPLEX):
	$(MAKE) -C $@ real
	$(MAKE) -C $@ complex

$(SIMPLE):
	$(MAKE) -C $@

test: target-dir = $(WIENROOT_TEST)
test: install

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
		rm -f $$dir/Makefile.orig; \
	done

distclean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir distclean; \
	done

%/Makefile.orig: %/Makefile
	perl -pe 's/^\#.orig\#//' $^ >$@


### ‘install’ target
###
### Installs what version 3.0 adds to or changes in a WIEN2k
### installation: the programs w2w and w2wc, the standalone
### w2w_hybrid (and w2w_hybrid_acc when it was built), and the driver.
### It installs no script: WIEN2k's own `x' script and helper scripts
### stay as they are.
target-dir ?= $(WIENROOT)
cli = SRC/wien2wannier
programs = SRC_w2w/w2w SRC_w2w/w2wc
hybrid   = $(wildcard SRC_w2w/w2w_hybrid SRC_w2w/w2w_hybrid_acc)

install-cli:
	@test -n "$(target-dir)" || { echo "target-dir or WIENROOT is required" >&2; exit 2; }
	install -d "$(target-dir)"
	install -m755 "$(cli)" "$(target-dir)/wien2wannier"

install: SRC_w2w install-cli
	for f in $(programs) $(hybrid); do \
		install -m755 "$$f" "$(target-dir)/$$(basename $$f)" || exit 1; \
	done
