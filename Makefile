build/Summary.pdf : Summary.tex
	@mkdir -p build
	tectonic --print --outdir build Summary.tex

.PHONY: clean
clean:
	rm -rf build
