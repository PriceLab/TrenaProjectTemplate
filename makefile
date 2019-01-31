all:  docs install

docs:
	R -e "devtools::document()"
build:
	(cd ..; R CMD build --no-build-vignettes trenaProjectTemplate)

install:
	(cd ..; R CMD INSTALL trenaProjectTemplate)

check:
	(cd ..; R CMD check `ls -t trenaProjectTemplate) | head -1`)

test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

