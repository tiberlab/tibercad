SHELL=/bin/sh

SOURCES=$(wildcard *.tex)

SUBDIRS=

RECURSIVEPDF=$(addsuffix .makepdf,$(SUBDIRS))
RECURSIVEPS=$(addsuffix .makeps,$(SUBDIRS))

MANUAL=tiber_manual

.PHONY : clean


all: $(MANUAL).dvi
#all: $(MANUAL).pdf 

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

clean:  clean-local
	@for subdir in ${SUBDIRS} ; do	\
		$(MAKE) -C $$subdir $@; \
	done

clean-local:
	rm -f *.ind *.idx *.aux *.lof *.log *.lot *.dvi *.out \
	      *.toc *.ilg *.bbl *.blg *~

distclean: clean
	rm -f *.pdf *.ps

	
