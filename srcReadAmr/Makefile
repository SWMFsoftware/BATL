#^CFG COPYRIGHT UM

include ../Makefile.def

SEARCH_EXTRA = -I../src

include ../Makefile.conf

LIBBATL = ../src/libBATL.a

OBJECTS = ModReadAmr.o wrap_read_amr.o

${OBJECTS}: ${LIBBATL}

read_amr.o: ${OBJECTS}

MY_LIB = ${LIBDIR}/libREADAMR.a

LIB: ${MY_LIB}

${MY_LIB}: ${LIBBATL} ${OBJECTS}
	cp -f ${LIBBATL} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

LIBSO: ${LIBBATL} ${OBJECTS}
	${LINK.f90} -o ${LIBDIR}/libREADAMR.so ${OBJECTS} \
		-L ../src -l BATL -L ${LIBDIR} -l SHARE -l TIMING \
		-shared -fPIC
MY_EXE = ${BINDIR}/READAMR.exe

EXE:    ${MY_EXE}

${MY_EXE}: LIB read_amr.o
	${LINK.f90} -o ${MY_EXE} read_amr.o \
		-L${LIBDIR} -lREADAMR -lTIMING -lSHARE ${Lflag1}

clean: cleanfiles

distclean: clean
