
SHELL =/bin/sh

include ../Makefile.def

include ../Makefile.conf

include Makefile.DEPEND

install:
	touch Makefile.DEPEND

OBJECTS = \
	main.o \
	BATL_size.o \
	BATL_tree.o \
	BATL_grid.o \
	BATL_geometry.o 

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

EXE:
	make DEPEND
	make BATL.exe

BATL.exe: ${OBJECTS}
	${LINK.f90} -o BATL.exe ${OBJECTS} ${Lflag2}

test: EXE
	BATL.exe > test.log
	-@diff test.log test.ref > test.diff
	ls -l test.diff

distclean: clean
	rm -f *.log *.diff *.rst BATL_size.f90 Makefile.DEPEND

