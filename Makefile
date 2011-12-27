include Makefile.def

install:
	touch src/Makefile.DEPEND
	cp -f src/BATL_size_orig.f90 src/BATL_size.f90

run:
	mkdir -p run/plots
	cd run; \
		cp ../input/PARAM.in.cart PARAM.in; \
		ln -s ${BINDIR}/ADVECT.exe .; \
		ln -s ${BINDIR}/BATL.exe .; \
		ln -s ${BINDIR}/PostIDL.exe .; \
		ln -s ${BINDIR}/pIDL .

BATL:	run
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make BATL

ADVECT: run
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make ADVECT

NOMPI:
	cd util/NOMPI/src; make LIB

test:	test_unit test_advect

test_unit: test11 test12 test21 test22 test31 test32 test33
	ls -l test??.diff

test_advect:	test_advect11 test_advect21 test_advect22 \
		test_advect12 \
		test_advect31 test_advect32 test_advect33 \
		test_advect22_rz advect22_cyl
	ls -l advect??*.diff

MPIRUN = mpirun -np 2

test11:
	Config.pl -g=8,1,1 -r=2,2,2
	make BATL
	cd run; ${MPIRUN} BATL.exe > test11.ref
	-@(${SCRIPTDIR}/DiffNum.pl -t \
		run/test11.ref output/test11.ref > test11.diff)
	ls -l test11.diff

test12:
	Config.pl -g=8,4,1 -r=1,2,1
	make BATL
	cd run; ${MPIRUN} BATL.exe > test12.ref
	-@(${SCRIPTDIR}/DiffNum.pl -t \
		run/test12.ref output/test12.ref > test12.diff)
	ls -l test12.diff

test21:
	Config.pl -g=8,4,1 -r=2,1,1
	make BATL
	cd run; ${MPIRUN} BATL.exe > test21.ref
	-@(${SCRIPTDIR}/DiffNum.pl -t \
		run/test21.ref output/test21.ref > test21.diff)
	ls -l test21.diff

test22:
	Config.pl -g=8,4,1 -r=2,2,2
	make BATL
	cd run; ${MPIRUN} BATL.exe > test22.ref
	-@(${SCRIPTDIR}/DiffNum.pl -t \
		run/test22.ref output/test22.ref > test22.diff)
	ls -l test22.diff
	sleep 1

test31:
	Config.pl -g=8,4,2 -r=2,1,1
	make BATL
	cd run; ${MPIRUN} BATL.exe > test31.ref
	-@(${SCRIPTDIR}/DiffNum.pl -t \
		run/test31.ref output/test31.ref > test31.diff)
	ls -l test31.diff
	sleep 1

test32:
	Config.pl -g=8,4,2 -r=2,2,1
	make BATL
	cd run; ${MPIRUN} BATL.exe > test32.ref
	-@(${SCRIPTDIR}/DiffNum.pl -t \
		run/test32.ref output/test32.ref > test32.diff)
	ls -l test32.diff
	sleep 1

test33:
	Config.pl -g=8,6,4 -r=2,2,2
	make BATL
	cd run; ${MPIRUN} BATL.exe > test33.ref
	-@(${SCRIPTDIR}/DiffNum.pl -t \
		run/test33.ref output/test33.ref > test33.diff)
	ls -l test33.diff
	sleep 1

test_advect11: 
	Config.pl -g=4,1,1 -r=2,2,2
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect11.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	cd run; ${MPIRUN} ADVECT.exe > runlog; mv advect.log advect11.log
	make test_advect11_check

test_advect12:
	Config.pl -g=4,4,1 -r=1,2,1
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect12.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	cd run; ${MPIRUN} ADVECT.exe > runlog; mv advect.log advect12.log
	make test_advect12_check

test_advect21:
	Config.pl -g=4,4,1 -r=2,1,1
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect21.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	cd run; ${MPIRUN} ADVECT.exe > runlog; mv advect.log advect21.log
	make test_advect21_check

test_advect22: 
	Config.pl -g=4,4,1 -r=2,2,2
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect22.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	cd run; ${MPIRUN} ADVECT.exe > runlog; mv advect.log advect22.log
	make test_advect22_check

test_advect22_rz: 
	Config.pl -g=4,4,1 -r=2,2,2
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect22_rz.log
	rm -f input/PARAM.in; cp input/PARAM.in.rz run/PARAM.in
	cd run; ${MPIRUN} ADVECT.exe > runlog; mv advect.log advect22_rz.log
	make test_advect22_rz_check

test_advect22_cyl: 
	Config.pl -g=4,4,1 -r=2,2,2
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect22_cyl.log
	rm -f input/PARAM.in; cp input/PARAM.in.cyl run/PARAM.in
	cd run; ${MPIRUN} ADVECT.exe > runlog; mv advect.log advect22_cyl.log
	make test_advect22_cyl_check

test_advect31: 
	Config.pl -g=4,4,4 -r=2,1,1
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect31.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	cd run; ${MPIRUN} ADVECT.exe > runlog; mv advect.log advect31.log
	make test_advect31_check

test_advect32:
	Config.pl -g=4,4,4 -r=2,2,1
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect32.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	cd run; ${MPIRUN} ADVECT.exe > runlog; mv advect.log advect32.log
	make test_advect32_check

test_advect33: 
	Config.pl -g=4,4,4 -r=2,2,2
	make ADVECT
	rm -rf run/plots/* run/runlog run/advect33.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	cd run; ${MPIRUN} ADVECT.exe > runlog; mv advect.log advect33.log
	make test_advect33_check

test_advect11_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect11.log output/advect11.log > advect11.diff)
	ls -l advect11.diff

test_advect12_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect12.log output/advect12.log > advect12.diff)
	ls -l advect12.diff

test_advect21_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect21.log output/advect21.log > advect21.diff)
	ls -l advect21.diff

test_advect22_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect22.log output/advect22.log > advect22.diff)
	ls -l advect22.diff

test_advect22_rz_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=5.e-6 -a=1.e-12 \
		run/advect22_rz.log output/advect22_rz.log > advect22_rz.diff)
	ls -l advect22_rz.diff

test_advect22_cyl_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=5.e-6 -a=1.e-12 \
		run/advect22_cyl.log output/advect22_cyl.log > advect22_cyl.diff)
	ls -l advect22_cyl.diff

test_advect31_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect31.log output/advect31.log > advect31.diff)
	ls -l advect31.diff

test_advect32_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect32.log output/advect32.log > advect32.diff)
	ls -l advect32.diff

test_advect33_check:
	-@(${SCRIPTDIR}/DiffNum.pl -r=1.e-8 -a=1.e-12 \
		run/advect33.log output/advect33.log > advect33.diff)
	ls -l advect33.diff

clean:
	cd share; make clean
	cd src; make clean

allclean:
	touch src/Makefile.DEPEND
	cd src; make distclean
	rm -rf run bin/*.exe *.diff

