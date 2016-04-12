filename=main

pdf: dvi
	pdflatex ${filename}

dvi:
	latex ${filename}
	bibtex ${filename}
	latex ${filename}
	bibtex ${filename}
	latex ${filename}

clean:
	rm ${filename}.{dvi,ps,aux,log,bbl,blg,pdf}
