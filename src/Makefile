
SHELL =/bin/sh

include ../../Makefile.def

include ../../Makefile.conf

include Makefile.DEPEND

install:
	touch Makefile.DEPEND

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

OBJECTS = \
	test_bamr.o \
	BAMR_size.o \
	BAMR_tree.o \
	BAMR_grid.o

EXE:
	make DEPEND
	make test_bamr.exe

test_bamr.exe: ${OBJECTS}
	${LINK.f90} -o test_bamr.exe ${OBJECTS} ${Lflag2}

test: EXE
	@echo './test_bamr.exe'
	@./test_bamr.exe > test_bamr.log
	-@diff test_bamr.log test_bamr.ref > test_bamr.diff
	ls -l test_bamr.diff
	@echo test_bamr done

distclean: clean
	rm -f *.log *.diff *.rst Makefile.DEPEND

