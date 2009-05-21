include Makefile.def

install:
	touch src/Makefile.DEPEND
	cp -f src/BATL_size_orig.f90 src/BATL_size.f90

test:
	cd src; make test

clean:
	cd src; make clean

allclean:
	touch src/Makefile.DEPEND
	cd src; make distclean
