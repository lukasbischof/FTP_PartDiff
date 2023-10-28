Summary.pdf : Summary.tex
	tectonic --print Summary.tex

.PHONY: clean
clean:
	rm -f Summary.pdf
