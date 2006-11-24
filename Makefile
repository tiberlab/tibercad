# $Id$
#

topdir := .
dirs := src examples

include Make.common


.PHONY: all doc clean



all:
	@$(MAKE) -C src all
	@$(CXX) -shared -fPIC $(CXXFLAGS) \
		$(TIBER_OBJECTS) -o lib/libtibercad.so $(LIBS) $(LDFLAGS)
	

clean:
	@for i in $(dirs);	\
	do			\
	  $(MAKE) -C $$i clean;	\
	done

distclean: clean
	@$(MAKE) -C examples distclean
	@$(MAKE) -C doc clean
	@rm -f lib/*


examples: all
	@$(MAKE) -C examples


doc:
	@$(MAKE) -C doc



tags:
	@echo $(tibermodules)
	@rm h.tags C.tags
	@touch h.tags C.tags
	@for i in $(tibermodules); \
	do \
	  etags -a --members src/$$i/*.C -o C.tags ; \
	  etags -a --members include/$$i/*.h -o h.tags ; \
	done 

