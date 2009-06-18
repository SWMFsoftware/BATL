
SHELL =/bin/sh

include ../Makefile.def

include ../Makefile.conf

include Makefile.DEPEND

install:
	touch Makefile.DEPEND

OBJECTS = \
	BATL_lib.o  \
	BATL_mpi.o  \
	BATL_size.o \
	BATL_tree.o \
	BATL_grid.o \
	BATL_geometry.o \
	BATL_pass_cell.o

ALLOBJECTS = \
	${OBJECTS} \
	main.o \
	advect_main.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${ALLOBJECTS}

BATL:
	make DEPEND
	make BATL.exe

BATL.exe: main.o ${OBJECTS}
	${LINK.f90} -o BATL.exe main.o ${OBJECTS} \
		-L${LIBDIR} -lSHARE -lTIMING ${Lflag1}

ADVECT:
	make DEPEND
	make ADVECT.exe

ADVECT.exe: advect_main.o ${OBJECTS}
	${LINK.f90} -o ADVECT.exe advect_main.o ${OBJECTS} \
		-L${LIBDIR} -lSHARE -lTIMING ${Lflag1}

distclean: clean
	rm -f *.log *.diff *.rst BATL_size.f90 Makefile.DEPEND

