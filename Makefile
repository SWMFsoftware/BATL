include Makefile.def

install:
	touch src/Makefile.DEPEND

test:
	cd src; make test

clean:
	cd src; make clean

allclean:
	touch src/Makefile.DEPEND
	cd src; make distclean
