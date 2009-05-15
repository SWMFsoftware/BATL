include Makefile.def

install:
	touch src/Makefile.DEPEND

test:
	cd src; make test

allclean:
	touch src/Makefile.DEPEND
	cd src; make distclean
