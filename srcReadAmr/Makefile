#^CFG COPYRIGHT UM

include ../Makefile.def

SEARCH_EXTRA = -I../src

include ../Makefile.conf

LIBBATL = ../src/libBATL.a

OBJECTS = ModReadAmr.o wrap_read_amr.o

${OBJECTS}: ${LIBBATL}

wrap_read_amr.o: ModReadAmr.o
read_amr.o: ${OBJECTS}

MY_LIB = ${LIBDIR}/libREADAMR.a

LIB: ${MY_LIB}

${MY_LIB}: ${LIBBATL} ${OBJECTS}
	cp -f ${LIBBATL} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

LIBSO: ${LIBBATL} ${OBJECTS}
	${LINK.f90} -o ${LIBDIR}/libWRAPAMR.so ${OBJECTS} \
		-L ../src -lBATL -L${LIBDIR} -lSHARE -lTIMING ${MPILIB} \
		-shared -fPIC

MY_EXE = ${BINDIR}/READAMR.exe

EXE:    ${MY_EXE}

${MY_EXE}: LIB read_amr.o
	${LINK.f90} -o ${MY_EXE} read_amr.o \
		-L${LIBDIR} -lREADAMR -lTIMING -lSHARE ${LflagMpi}

MY_EXE_C = ${BINDIR}/READAMR_C.exe

EXE_C:  ${MY_EXE_C}

# Link with the NOMPI library. 
${MY_EXE_C}: LIB read_amr_c.o
	${LINK.f90} ${NOMAINFLAG} -o ${MY_EXE_C} read_amr_c.o \
		-L${LIBDIR} -lREADAMR -lTIMING -lSHARE ${MPILIB} ${Lflag}

clean: cleanfiles

distclean: clean
	rm -f readamr_py.ref
