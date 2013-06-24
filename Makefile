# fix folders 
all:
	mv _static static
	mv _sources sources
	sed -i 's/_static/static/g' *html
	sed -i 's/_sources/sources/g' *html
