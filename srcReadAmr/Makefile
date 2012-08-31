#^CFG COPYRIGHT UM

include ../Makefile.def

SEARCH_EXTRA = -I../src

include ../Makefile.conf

OBJECTS= ModReadAmr.o

MY_LIB = ${LIBDIR}/libREADAMR.a

EXE = ${BINDIR}/READAMR.exe

LIB:
	make ${MY_LIB}

LIBBATL = ../src/libBATL.a

LIB: ${LIBBATL} ${OBJECTS}
	@cp -f ${LIBBATL} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

EXE:    ${EXE}

${EXE}: LIB read_amr.o
	${LINK.f90} -o ${EXE} read_amr.o \
		-L${LIBDIR} -lREADAMR -lTIMING -lSHARE ${LFLAG2}

distclean: clean
