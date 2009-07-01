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
	ls -l test??.diff

test_advect:	test_advect11 test_advect21 test_advect22 test_advect31 test_advect32 test_advect33
	ls -l advect??.diff

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

test_advect11: 
	Config.pl -g=1,1,4
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect11.log
	cd run; ${MPIRUN} ADVECT.exe > runlog
	make test_advect11_check

test_advect21:
	Config.pl -g=2,1,4,4
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect21.log
	cd run; ${MPIRUN} ADVECT.exe > runlog
	make test_advect21_check

test_advect22: 
	Config.pl -g=2,2,4,4
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect22.log
	cd run; ${MPIRUN} ADVECT.exe > runlog
	make test_advect22_check

test_advect31: 
	Config.pl -g=3,1,4,4,4
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect31.log
	cd run; ${MPIRUN} ADVECT.exe > runlog
	make test_advect31_check

test_advect32:
	Config.pl -g=3,2,4,4,4
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect32.log
	cd run; ${MPIRUN} ADVECT.exe > runlog
	make test_advect32_check

test_advect33: 
	Config.pl -g=3,3,4,4,4
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect33.log
	cd run; ${MPIRUN} ADVECT.exe > runlog
	make test_advect33_check

test_advect11_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect11.log src/advect11.ref > advect11.diff)
	ls -l advect11.diff

test_advect21_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect21.log src/advect21.ref > advect21.diff)
	ls -l advect21.diff

test_advect22_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect22.log src/advect22.ref > advect22.diff)
	ls -l advect22.diff

test_advect31_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect31.log src/advect31.ref > advect31.diff)
	ls -l advect31.diff

test_advect32_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect32.log src/advect32.ref > advect32.diff)
	ls -l advect32.diff

test_advect33_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect33.log src/advect33.ref > advect33.diff)
	ls -l advect33.diff

clean:
	cd share; make clean
	cd src; make clean

allclean:
	touch src/Makefile.DEPEND
	cd src; make distclean
	rm -rf run bin *.diff

