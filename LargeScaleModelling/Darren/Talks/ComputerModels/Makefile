# Makefile


lhd.pdf: lhd.R
	Rscript lhd.R
	mv Rplots.pdf lhd.pdf

sim1.pdf: gp.R
	Rscript gp.R

talk.pdf: talk.tex sim1.pdf lhd.pdf
	pdflatex talk
	pdflatex talk

view: talk.pdf
	acroread talk.pdf &

edit:
	gedit talk.tex gp.R &


clean:
	rm -f *~ core *.pdf *.ps *.dvi *.log *.aux *.out *.flc

# eof
