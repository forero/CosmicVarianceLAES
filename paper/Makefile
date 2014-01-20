LATEX = pdflatex
MAIN  = ms


$(MAIN).pdf: 	$(MAIN).tex  references.bib
#	$(LATEX) $(MAIN).tex 
#	bibtex $(MAIN)
	$(LATEX) $(MAIN).tex 
	$(LATEX) $(MAIN).tex 

clean:
	rm -f $(MAIN)~
	rm -f $(MAIN).aux
	rm -f $(MAIN).bbl
	rm -f $(MAIN).ps
	rm -f $(MAIN).pdf
	rm -f $(MAIN).dvi
	rm -f $(MAIN).blg
	rm -f $(MAIN).log
	rm -f *~
	rm -f $(MAIN).tex.flc

submit:
	tar -cvzf draft.tgz Fig* ms.tex references.bib ms.pdf

view: 
	open $(MAIN).pdf &

