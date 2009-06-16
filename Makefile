include Makefile.def

install:
	touch src/Makefile.DEPEND
	cp -f src/BATL_size_orig.f90 src/BATL_size.f90

BATL:
	cd ${SHAREDIR}; make LIB
	cd src; make BATL

ADVECT:
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make ADVECT

test:	test1 test2 test3
	ls -l src/*.diff

MPIRUN = mpirun -np 2

test1:
	Config.pl -g=3,1,8,4,2
	make BATL
	cd src; ${MPIRUN} BATL.exe > test1.log
	-@(cd src; diff test1.log test1.ref > test1.diff)
	cd src; ls -l test1.diff

test2:
	Config.pl -g=3,2,8,4,2
	make BATL
	cd src; ${MPIRUN} BATL.exe > test2.log
	-@(cd src; diff test2.log test2.ref > test2.diff)
	cd src; ls -l test2.diff

test3:
	Config.pl -g=3,3,8,4,2
	make BATL
	cd src; ${MPIRUN} BATL.exe > test3.log
	-@(cd src; diff test3.log test3.ref > test3.diff)
	cd src; ls -l test3.diff

test_advect2:
	Config.pl -g=2,2,80,40
	make ADVECT
	cd src; ${MPIRUN} ADVECT.exe > advect2.log
	make test_advect2_check

test_advect2_check:
	-@(cd src; ${SCRIPTDIR}/DiffNum.pl \
		-r=1.e-8 advect2.log advect2.ref > advect2.diff)
	ls -l src/advect2.diff

clean:
	cd share; make clean
	cd src; make clean

allclean:
	touch src/Makefile.DEPEND
	cd src; make distclean
