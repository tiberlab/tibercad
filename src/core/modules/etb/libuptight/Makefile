# $Id$

topdir := ../../../../..

include $(topdir)/Make.defs
include Make.defs


ARCH = $(shell uname -m)

CXXFLAGS = -fPIC -I$(MPIINCL) -I . 

UPTREV=HEAD
UPTREPO=https://afrodite.eln.uniroma2.it/svn/uptight/branches/tibercad/src

target := libuptight$(shared_libext)


libuptight$(shared_libext) :  .$(ARCH)/UptWrapper$(obj-suffix) libuptight_$(ARCH).a
	@echo "Linking libuptight$(shared_libext) ..."
	@$(CXX) -shared -o $@ $< \
		-Wl,--version-script,export.map \
		-Wl,-Bsymbolic -Wl,-Bstatic \
		-L. libuptight_$(ARCH).a \
		-Wl,-Bdynamic
	
.$(ARCH)/UptWrapper$(obj-suffix) : libuptight_$(ARCH).a 

libuptight_$(ARCH).a: checkout
	@echo "Updating libuptight ..."
	@cd src; svn -r $(UPTREV) up lib_uptight; svn -r $(UPTREV) up api; svn -r $(UPTREV) up make.common; \
	echo "Compiling libuptight ..."; \
	$(MAKE) -C api libuptight; cd ..; 
	@mv src/api/libuptight_$(ARCH).a .	
	@cp src/api/fortran.h .
	@cp src/api/uptight.h . 
	
checkout:
	@echo "Checkout libuptight ..."
	@if [ ! -d src/api ]; then svn -r $(UPTREV) co $(UPTREPO) --depth empty; fi
	

	

include $(topdir)/Make.rules
