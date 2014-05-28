# Makefile

talk.pdf: talk.tex
	pdflatex talk
	pdflatex talk

view: talk.pdf
	acroread talk.pdf &

edit:
	gedit talk.tex &


clean:
	rm -f *~ core *.pdf *.ps *.dvi *.log *.aux *.out *.flc

# eof
