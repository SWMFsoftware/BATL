include Makefile.def

install:
	touch src/Makefile.DEPEND
	cp -f src/BATL_size_orig.f90 src/BATL_size.f90

BATL:
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make BATL

ADVECT:
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make ADVECT

test:	test11 test21 test22 test31 test32 test33
	ls -l src/*.diff

MPIRUN = mpirun -np 2

test11:
	Config.pl -g=1,1,8
	make BATL
	cd src; ${MPIRUN} BATL.exe > test11.log
	-@(cd src; diff test11.log test11.ref > test11.diff)
	cd src; ls -l test11.diff

test21:
	Config.pl -g=2,1,8,4
	make BATL
	cd src; ${MPIRUN} BATL.exe > test21.log
	-@(cd src; diff test21.log test21.ref > test21.diff)
	cd src; ls -l test21.diff

test22:
	Config.pl -g=2,2,8,4
	make BATL
	cd src; ${MPIRUN} BATL.exe > test22.log
	-@(cd src; diff test22.log test22.ref > test22.diff)
	cd src; ls -l test22.diff

test31:
	Config.pl -g=3,1,8,4,2
	make BATL
	cd src; ${MPIRUN} BATL.exe > test31.log
	-@(cd src; diff test31.log test31.ref > test31.diff)
	cd src; ls -l test31.diff

test32:
	Config.pl -g=3,2,8,4,2
	make BATL
	cd src; ${MPIRUN} BATL.exe > test32.log
	-@(cd src; diff test32.log test32.ref > test32.diff)
	cd src; ls -l test32.diff

test33:
	Config.pl -g=3,3,8,4,2
	make BATL
	cd src; ${MPIRUN} BATL.exe > test33.log
	-@(cd src; diff test33.log test33.ref > test33.diff)
	cd src; ls -l test33.diff

rundir:
	mkdir -p run/plots
	cd run; ln -s ../src/ADVECT.exe .;

test_advect22: 
	Config.pl -g=2,2,80,40
	make ADVECT
	rm -rf run
	make rundir
	cd run; ${MPIRUN} ADVECT.exe > runlog
	make test_advect22_check

test_advect22_check:
	-@(${SCRIPTDIR}/DiffNum.pl \
		-r=1.e-8 run/advect22.log src/advect22.ref > advect22.diff)
	ls -l advect22.diff

clean:
	cd share; make clean
	cd src; make clean

allclean:
	touch src/Makefile.DEPEND
	cd src; make distclean
	rm -rf run *.diff

