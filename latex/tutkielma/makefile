name=tutkielma
bibfile=bibliography.bib

all: $(name).pdf
	
simple: 
	pdflatex $(name).tex
	echo ???
	
clean:
	rm -f *.aux *.blg *.out *.bbl *.log *.nlo *.nls *.toc *.ilg $(name).pdf
	echo !!!
	
$(name).pdf: $(name).tex $(bibfile)
	pdflatex $(name).tex
	bibtex $(name).aux
	makeindex $(name).nlo -s nomencl.ist -o $(name).nls
	pdflatex $(name).tex
	pdflatex $(name).tex

