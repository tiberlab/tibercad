topdir=..

# the subversion revision number
SVNREVISION = $(shell svn info | awk '/Revision/ {print $$2}')

# the full TiberCAD version number
TIBERVERSION = $(shell echo "`cat $(topdir)/VERSION`-$(SVNREVISION)")


SOURCES=$(wildcard *.tex) version.tex

SUBDIRS=figures

RECURSIVEPDF=$(addsuffix .makepdf,$(SUBDIRS))
RECURSIVEPS=$(addsuffix .makeps,$(SUBDIRS))

MANUAL=tiber_manual

.PHONY : clean


all: $(MANUAL).dvi

manual: pdf 

pdf: $(MANUAL).pdf

ps: $(MANUAL).ps


$(RECURSIVEPDF): $(SUBDIRS)
	$(MAKE) -C $(basename $@) pdf

$(RECURSIVEPS): $(SUBDIRS)
	$(MAKE) -C $(basename $@) ps


%.dvi: %.tex ${SOURCES} clean-local
	latex $< 
	bibtex $(basename $<)
	latex $< 
	latex $< 

%.ps:  %.dvi $(RECURSIVEPS)
	dvips $< -o $@ 

%.pdf: %.tex ${SOURCES} $(RECURSIVEPDF) clean-local
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<

.PHONY: version.tex
version.tex:
	@echo $(TIBERVERSION) > $@

clean:  clean-local
	@for subdir in ${SUBDIRS} ; do	\
		$(MAKE) -C $$subdir $@; \
	done

clean-local:
	-rm -f *.ind *.idx *.aux *.lof *.log *.lot *.dvi *.out \
	      *.toc *.ilg *.bbl *.blg *~

distclean: clean
	-rm -f *.pdf *.ps

	
