# Makefile for volta.  Nothing has to be built in advance: the `volta` script
# compiles a binary for each parameter set on first use.  See README.md.
#
#   make check                  run the test suite (about a minute)
#   make install [PREFIX=...]   install the script and the sources (default /usr/local)
#   make uninstall [PREFIX=...]
#   make clean-cache            remove the compiled binaries kept by volta

PREFIX  ?= /usr/local
BINDIR   = $(PREFIX)/bin
DATADIR  = $(PREFIX)/share/volta

.PHONY: all check install uninstall clean-cache

all:
	@echo "volta needs no build step: run ./volta ORBITS [options], e.g. ./volta 5,5"
	@echo "  make check                  run the test suite"
	@echo "  make install [PREFIX=...]   install into $(PREFIX)"

check:
	@sh tests/check.sh

install:
	install -d $(DESTDIR)$(BINDIR) $(DESTDIR)$(DATADIR)/src $(DESTDIR)$(DATADIR)/doc
	install -m 755 volta $(DESTDIR)$(BINDIR)/volta
	install -m 644 src/volta.cpp src/bitset.h $(DESTDIR)$(DATADIR)/src
	install -m 644 README.md LICENSE doc/DESIGN.md $(DESTDIR)$(DATADIR)/doc

uninstall:
	rm -f $(DESTDIR)$(BINDIR)/volta
	rm -f $(DESTDIR)$(DATADIR)/src/volta.cpp $(DESTDIR)$(DATADIR)/src/bitset.h
	rm -f $(DESTDIR)$(DATADIR)/doc/README.md $(DESTDIR)$(DATADIR)/doc/LICENSE $(DESTDIR)$(DATADIR)/doc/DESIGN.md
	-rmdir $(DESTDIR)$(DATADIR)/src $(DESTDIR)$(DATADIR)/doc $(DESTDIR)$(DATADIR)

clean-cache:
	./volta --clear-cache
