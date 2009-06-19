include Makefile.def

install:
	touch src/Makefile.DEPEND
	cp -f src/BATL_size_orig.f90 src/BATL_size.f90

bin:
	mkdir -p bin

run:
	mkdir -p run/plots
	cd run; \
		ln -s ${BINDIR}/ADVECT.exe .; \
		ln -s ${BINDIR}/BATL.exe .; \
		ln -s ${BINDIR}/PostIDL.exe .; \
		ln -s ${BINDIR}/pIDL .

BATL:	bin run
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make BATL

ADVECT: bin run
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make ADVECT

test:	test11 test21 test22 test31 test32 test33
	ls -l *.diff

MPIRUN = mpirun -np 2

test11:
	Config.pl -g=1,1,8
	make BATL
	cd run; ${MPIRUN} BATL.exe > test11.ref
	-@(${SCRIPTDIR}/DiffNum.pl -b \
		run/test11.ref src/test11.ref > test11.diff)
	ls -l test11.diff

test21:
	Config.pl -g=2,1,8,4
	make BATL
	cd run; ${MPIRUN} BATL.exe > test21.ref
	-@(${SCRIPTDIR}/DiffNum.pl -b \
		run/test21.ref src/test21.ref > test21.diff)
	ls -l test21.diff

test22:
	Config.pl -g=2,2,8,4
	make BATL
	cd run; ${MPIRUN} BATL.exe > test22.ref
	-@(${SCRIPTDIR}/DiffNum.pl -b \
		run/test22.ref src/test22.ref > test22.diff)
	ls -l test22.diff

test31:
	Config.pl -g=3,1,8,4,2
	make BATL
	cd run; ${MPIRUN} BATL.exe > test31.ref
	-@(${SCRIPTDIR}/DiffNum.pl -b \
		run/test31.ref src/test31.ref > test31.diff)
	ls -l test31.diff

test32:
	Config.pl -g=3,2,8,4,2
	make BATL
	cd run; ${MPIRUN} BATL.exe > test32.ref
	-@(${SCRIPTDIR}/DiffNum.pl -b \
		run/test32.ref src/test32.ref > test32.diff)
	ls -l test32.diff

test33:
	Config.pl -g=3,3,8,6,4
	make BATL
	cd run; ${MPIRUN} BATL.exe > test33.ref
	-@(${SCRIPTDIR}/DiffNum.pl -b \
		run/test33.ref src/test33.ref > test33.diff)
	ls -l test33.diff

test_advect22: 
	Config.pl -g=2,2,80,40
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect22.log
	cd run; ${MPIRUN} ADVECT.exe > runlog
	make test_advect22_check

test_advect22_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect22.log src/advect22.ref > advect22.diff)
	ls -l advect22.diff

clean:
	cd share; make clean
	cd src; make clean

allclean:
	touch src/Makefile.DEPEND
	cd src; make distclean
	rm -rf run bin *.diff

