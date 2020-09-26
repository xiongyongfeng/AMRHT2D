# make latex para compilar
# aryuivos necessarios els-manuscript-v# and mybibfile
# pacote necessario: texlive-publishers

ARQUIVO=manuscript-amrht2d-v1A
ARQUIVO-BIB=mybibfile

latex: $(ARQUIVO).tex $(ARQUIVO-BIB).bib
	pdflatex $(ARQUIVO)
	bibtex $(ARQUIVO)
	pdflatex $(ARQUIVO)
	pdflatex $(ARQUIVO)

clean:
	rm $(ARQUIVO).pdf $(ARQUIVO).log $(ARQUIVO).aux $(ARQUIVO).bbl $(ARQUIVO).blg $(ARQUIVO).out $(ARQUIVO).spl 
