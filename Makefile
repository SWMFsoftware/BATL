include Makefile.def

# BATL tests are designed for 2 cores
BATLRUN = ${PARALLEL} ${NPFLAG} 2

help:
	@echo "Available targets:"
	@echo "help                  - this help message"
	@echo "run                   - create run directory"
	@echo "BATL                  - create bin/BATL.exe for unit tests"
	@echo "ADVECT                - create bin/ADVECT.exe for advection tests"
	@echo "GAME                  - create bin/GAME.exe for game-of-life"
	@echo "READAMRLIB            - create lib/libREADAMR.a"
	@echo "WRAPAMRLIB            - create lib/libWRAPAMR.so"
	@echo "READAMR               - create bin/READAMR.exe for READAMR tests"
	@echo "READAMR_C             - create bin/READAMR_C.exe to test C wrapper"
	@echo "PIDL                  - create bin/PostIDL.exe for post processing"
	@echo "NOMPI                 - create lib/NOMPI.a for serial execution"
	@echo "clean                 - remove object files"
	@echo "test                  - run full test suite (optional MPIRUN=...)"
	@echo "test_unit             - run unit tests with BATL.exe"
	@echo "test_advect           - run advection tests with ADVECT.exe"
	@echo "test_readamr          - run READAMR tests with READAMR.exe"
	@echo "test11                - run 1D unit test"
	@echo "test11 BLESS=YES      - update test11 reference solution"
	@echo "test12                - 2D unit test with AMR in 2nd dimension"
	@echo "test21                - 2D unit test with AMR in 1st dimension"
	@echo "test22                - 2D unit test with AMR in all dimensions"
	@echo "test31                - 3D unit test with AMR in 1st dimension"
	@echo "test32                - 3D unit test with AMR in 1&2 dimensions"
	@echo "test33                - 3D unit test with AMR in all dimensions"
	@echo "test_advect11         - 1D advection test"
	@echo "test_advect12         - 2D advection test with AMR in 2nd dimension"
	@echo "test_advect21         - 2D advection test with AMR in 1st dimension"
	@echo "test_advect22         - 2D advection test with AMR in all dimensions"
	@echo "test_advect22_rot     - 2D advection test with rotated Cartesian grid"
	@echo "test_advect22_rz      - 2D advection test with R-Z geometry"
	@echo "test_advect22_cyl     - 2D advection test with cylindrical geometry"
	@echo "test_advect22_round   - 2D advection test with round-cube geometry"
	@echo "test_advect31         - 3D advection test with AMR in 1st dimension"
	@echo "test_advect32         - 3D advection test with AMR in 1&2 dimension"
	@echo "test_advect33         - 3D advection test with AMR in all dimensions"
	@echo "test_advect33_sph     - 3D advection test with spherical geometry"
	@echo "test_advect33_rlonlat - 3D advection test with R-Lon-Lat geometry"
	@echo "test_advect33_round   - 3D advection test with round-cube geometry"
	@echo "test_advect33_cubed   - 3D advection test with cubed-sphere geometry"
	@echo "test_readamr_1d       - 1D READAMR test"
	@echo "test_readamr_2d       - 2D READAMR test"
	@echo "test_readamr_3d       - 3D READAMR test"
	@echo "test_readamr_sph      - 3D READAMR test with logarithmic spherical grid"
	@echo "test_readamr_sphgenr  - 3D READAMR test with stretched spherical grid"
	@echo "test_readamr_c        - READAMR test with C wrapper"
	@echo "test_readamr_py       - READAMR test with Python wrapper"

install:
	mkdir -p lib bin
	cp -f src/BATL_size_orig.f90 src/BATL_size.f90

run:
	mkdir -p run/plots
	cd run; \
		cp ../input/PARAM.in.cart PARAM.in; \
		ln -s ${BINDIR}/ADVECT.exe .; \
		ln -s ${BINDIR}/BATL.exe .; \
		ln -s ${BINDIR}/GAME.exe .; \
		ln -s ${BINDIR}/READAMR.exe .; \
		ln -s ${BINDIR}/READAMR_C.exe .; \
		ln -s ${BINDIR}/PostIDL.exe .; \
		ln -s ${SCRIPTDIR}/pIDL .; \
		ln -s ${DIR}/data .

.NOTPARALLEL: BATL ADVECT GAME READAMRLIB WRAPAMRLIB READAMR READAMR_C


BATL:	run
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make BATL

BATL_dynamic:	run
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make LIBSO; make BATL

ADVECT: run
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make ADVECT

GAME: run
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make GAME

READAMRLIB:
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make LIB
	cd srcReadAmr; make LIB

WRAPAMRLIB:
	perl -pi -e 's/^(CFLAG =.*)/$$1 -fPIC/ unless /fPIC/' Makefile.conf
	make clean
	-(cd util/NOMPI/src; make LIB)
	-(cd ${SHAREDIR}; make LIB)
	-(cd ${TIMINGDIR}; make LIB)
	-(cd src; make LIB)
	-(cd srcReadAmr; make LIBSO)
	perl -pi -e 's/^(CFLAG =.*)( \-fPIC)+/$$1/' Makefile.conf

READAMR: run
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make LIB
	cd srcReadAmr; make EXE

READAMR_C: run
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make LIB
	cd srcReadAmr; make EXE_C

PIDL:
	cd ${SHAREDIR}; make PIDL

NOMPI:
	cd util/NOMPI/src; make LIB

.NOTPARALLEL: test test_unit test_readamr

# Setting BLESS=YES or Y will copy the solution into the reference solution
BLESS=NO

DIFFNUM = ${SCRIPTDIR}/DiffNum.pl -BLESS=${BLESS}

test:	
	make test_unit
	make test_advect
	make test_readamr

test_unit:
	@rm -f test??.diff
	-@(${MAKE} test11)
	-@(${MAKE} test12)
	-@(${MAKE} test21)
	-@(${MAKE} test22)
	-@(${MAKE} test31)
	-@(${MAKE} test31)
	-@(${MAKE} test32)
	-@(${MAKE} test33)
	ls -l test??.diff

test11:
	./Config.pl -default -g=10,1,1 -r=2,2,2 -ng=5
	-@(${MAKE} BATL)
	-(cd run; ${BATLRUN} ./BATL.exe > test11.ref)
	-@(${DIFFNUM} -t -a=1e-15 \
		run/test11.ref output/test11.ref > test11.diff)
	ls -l test11.diff

test12:
	./Config.pl -default -g=8,6,1 -r=1,2,1 -ng=3
	-@(${MAKE} BATL)
	-(cd run; ${BATLRUN} ./BATL.exe > test12.ref)
	-@(${DIFFNUM} -t -a=1e-15 \
		run/test12.ref output/test12.ref > test12.diff)
	ls -l test12.diff

test21:
	./Config.pl -default -g=8,4,1 -r=2,1,1 -ng=4
	-@(${MAKE} BATL)
	-(cd run; ${BATLRUN} ./BATL.exe > test21.ref)
	-@(${DIFFNUM} -t -a=1e-15 \
		run/test21.ref output/test21.ref > test21.diff)
	ls -l test21.diff

test22:
	./Config.pl -default -g=8,6,1 -r=2,2,2 -ng=3
	-@(${MAKE} BATL)
	-(cd run; ${BATLRUN} ./BATL.exe > test22.ref)
	-@(${DIFFNUM} -t -a=1e-15 \
		run/test22.ref output/test22.ref > test22.diff)
	ls -l test22.diff
	sleep 1

test31:
	./Config.pl -default -g=8,4,2 -r=2,1,1 -ng=2
	-@(${MAKE} BATL)
	-(cd run; ${BATLRUN} ./BATL.exe > test31.ref)
	-@(${DIFFNUM} -t -a=1e-15 \
		run/test31.ref output/test31.ref > test31.diff)
	ls -l test31.diff
	sleep 1

test32:
	./Config.pl -default -g=8,4,2 -r=2,2,1 -ng=2
	-@(${MAKE} BATL)
	-(cd run; ${BATLRUN} ./BATL.exe > test32.ref)
	-@(${DIFFNUM} -t -a=1e-15 \
		run/test32.ref output/test32.ref > test32.diff)
	ls -l test32.diff
	sleep 1

test33:
	./Config.pl -default -g=10,8,6 -r=2,2,2 -ng=3
	-@(${MAKE} BATL)
	-(cd run; ${BATLRUN} ./BATL.exe > test33.ref)
	-@(${DIFFNUM} -t -a=1e-10 \
		run/test33.ref output/test33.ref > test33.diff)
	ls -l test33.diff
	sleep 1

test_advect:
	@rm -f test_advect*.diff
	-@(${MAKE} test_advect11)
	-@(${MAKE} test_advect21)
	-@(${MAKE} test_advect22)
	-@(${MAKE} test_advect12)
	-@(${MAKE} test_advect31)
	-@(${MAKE} test_advect32)
	-@(${MAKE} test_advect33)
	-@(${MAKE} test_advect22_rz)
	-@(${MAKE} test_advect22_cyl)
	-@(${MAKE} test_advect22_round)
	-@(${MAKE} test_advect33_sph)
	-@(${MAKE} test_advect33_rlonlat)
	-@(${MAKE} test_advect33_round)
	-@(${MAKE} test_advect33_cubed)
	ls -l advect??*.diff

test_advect11:
	./Config.pl -default -g=4,1,1 -r=2,2,2 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect11.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect11.log)
	${MAKE} test_advect11_check

test_advect12:
	./Config.pl -default -g=4,4,1 -r=1,2,1 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect12.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect12.log)
	${MAKE} test_advect12_check

test_advect21:
	./Config.pl -default -g=4,4,1 -r=2,1,1 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect21.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect21.log)
	${MAKE} test_advect21_check

test_advect22: 
	./Config.pl -default -g=4,4,1 -r=2,2,2 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect22.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect22.log)
	${MAKE} test_advect22_check

test_advect22_rot:
	./Config.pl -default -g=4,4,1 -r=2,2,2 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect22.log
	rm -f input/PARAM.in; cp input/PARAM.in.rotcart run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect22.log)
	${MAKE} test_advect22_rot_check

test_advect22_rz: 
	./Config.pl -default -g=4,4,1 -r=2,2,2 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect22_rz.log
	rm -f input/PARAM.in; cp input/PARAM.in.rz run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect22_rz.log)
	${MAKE} test_advect22_rz_check

test_advect22_cyl: 
	./Config.pl -default -g=4,4,1 -r=2,2,2 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect22_cyl.log
	rm -f input/PARAM.in; cp input/PARAM.in.cyl run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect22_cyl.log)
	${MAKE} test_advect22_cyl_check

test_advect22_round: 
	./Config.pl -default -g=4,4,1 -r=2,2,2 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect22_round.log
	rm -f input/PARAM.in; cp input/PARAM.in.round run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect22_round.log)
	${MAKE} test_advect22_round_check

test_advect33_sph: 
	./Config.pl -default -g=4,4,4 -r=2,2,2 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect33_sph.log
	rm -f input/PARAM.in; cp input/PARAM.in.sph run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect33_sph.log)
	${MAKE} test_advect33_sph_check

test_advect33_rlonlat:
	./Config.pl -default -g=4,4,4 -r=2,2,2 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect33_rlonlat.log
	rm -f input/PARAM.in; cp input/PARAM.in.rlonlat run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect33_rlonlat.log)
	${MAKE} test_advect33_rlonlat_check

test_advect33_round: 
	./Config.pl -default -g=4,4,4 -r=2,2,2 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect33_round.log
	rm -f input/PARAM.in; cp input/PARAM.in.round3d run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect33_round.log)
	${MAKE} test_advect33_round_check

test_advect33_cubed: 
	./Config.pl -default -g=4,4,4 -r=2,2,2 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/cubed33_cubed.log
	rm -f input/PARAM.in; cp input/PARAM.in.cubed3d run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect33_cubed.log)
	${MAKE} test_advect33_cubed_check

test_advect31: 
	./Config.pl -default -g=4,4,4 -r=2,1,1 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect31.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect31.log)
	${MAKE} test_advect31_check

test_advect32:
	./Config.pl -default -g=4,4,4 -r=2,2,1 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect32.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect32.log)
	${MAKE} test_advect32_check

test_advect33: 
	./Config.pl -default -g=4,4,4 -r=2,2,2 -ng=2
	-@(${MAKE} ADVECT PIDL)
	rm -rf run/plots/* run/runlog run/advect33.log
	rm -f input/PARAM.in; cp input/PARAM.in.cart run/PARAM.in
	-(cd run; ${MPIRUN} ./ADVECT.exe > runlog; ./pIDL -M; \
		   mv advect.log advect33.log)
	${MAKE} test_advect33_check

test_advect11_check:
	-@(${DIFFNUM} -r=1e-8 -a=1e-12 \
		run/advect11.log output/advect11.log > advect11.diff)
	ls -l advect11.diff

test_advect12_check:
	-@(${DIFFNUM} -r=1e-8 -a=1e-11 \
		run/advect12.log output/advect12.log > advect12.diff)
	ls -l advect12.diff

test_advect21_check:
	-@(${DIFFNUM} -r=1e-8 -a=1e-12 \
		run/advect21.log output/advect21.log > advect21.diff)
	ls -l advect21.diff

test_advect22_check:
	-@(${DIFFNUM} -r=1e-8 -a=1e-12 \
		run/advect22.log output/advect22.log > advect22.diff)
	ls -l advect22.diff

test_advect22_rot_check:
	-@(${DIFFNUM} -r=1e-8 -a=1e-12 \
		run/advect22.log output/advect22.log > advect22_rot.diff)
	ls -l advect22_rot.diff

test_advect22_rz_check:
	-@(${DIFFNUM} -r=1e-5 -a=1e-12 \
		run/advect22_rz.log output/advect22_rz.log > advect22_rz.diff)
	ls -l advect22_rz.diff

test_advect22_cyl_check:
	-@(${DIFFNUM} -r=5e-6 -a=1e-12 \
		run/advect22_cyl.log output/advect22_cyl.log \
		> advect22_cyl.diff)
	ls -l advect22_cyl.diff

test_advect22_round_check:
	-@(${DIFFNUM} -r=5e-6 -a=1e-12 \
		run/advect22_round.log output/advect22_round.log \
		> advect22_round.diff)
	ls -l advect22_round.diff

test_advect31_check:
	-@(${DIFFNUM} -r=1e-8 -a=1e-12 \
		run/advect31.log output/advect31.log > advect31.diff)
	ls -l advect31.diff

test_advect32_check:
	-@(${DIFFNUM} -r=1e-8 -a=1e-12 \
		run/advect32.log output/advect32.log > advect32.diff)
	ls -l advect32.diff

test_advect33_check:
	-@(${DIFFNUM} -r=1e-8 -a=1e-12 \
		run/advect33.log output/advect33.log > advect33.diff)
	ls -l advect33.diff

test_advect33_sph_check:
	-@(${DIFFNUM} -r=5e-6 -a=1e-12 \
		run/advect33_sph.log output/advect33_sph.log \
		> advect33_sph.diff)
	ls -l advect33_sph.diff

test_advect33_rlonlat_check:
	-@(${DIFFNUM} -r=5e-6 -a=1e-12 \
		run/advect33_rlonlat.log output/advect33_rlonlat.log \
		> advect33_rlonlat.diff)
	ls -l advect33_rlonlat.diff

test_advect33_round_check:
	-@(${DIFFNUM} -r=5e-6 -a=1e-12 \
		run/advect33_round.log output/advect33_round.log \
		> advect33_round.diff)
	ls -l advect33_round.diff

test_advect33_cubed_check:
	-@(${DIFFNUM} -r=5e-6 -a=1e-12 \
		run/advect33_cubed.log output/advect33_cubed.log \
		> advect33_cubed.diff)
	ls -l advect33_cubed.diff

test_readamr: 
	rm -f readamr_*.diff
	-@(make -j1 test_readamr_1d)
	-@(make -j1 test_readamr_2d)
	-@(make -j1 test_readamr_3d)
	-@(make -j1 test_readamr_sph)
	-@(make -j1 test_readamr_sphgenr)
	-@(make -j1 test_readamr_c)
	if [ ${COMPILER} == gfortran ]; then make -j1 test_readamr_py; fi
	./Config.pl -mpi
	ls -l readamr_*.diff

test_readamr_1d:
	echo "test_readamr_1d_compile" > readamr_1d.diff
	./Config.pl -default -g=4,1,1 -r=2,1,1 -ng=0
	@(make READAMR)
	echo "test_readamr_1d_run" >> readamr_1d.diff
	(cd run; ${MPIRUN} ./READAMR.exe > readamr_1d.ref)
	-@(${DIFFNUM} -t \
		run/readamr_1d.ref output/readamr_1d.ref > readamr_1d.diff)
	ls -l readamr_1d.diff

test_readamr_2d:
	echo "test_readamr_2d_compile" > readamr_2d.diff
	./Config.pl -default -g=4,4,1 -r=2,2,1 -ng=2
	@(make READAMR)
	echo "test_readamr_2d_run" >> readamr_2d.diff
	(cd run; ${MPIRUN} ./READAMR.exe > readamr_2d.ref)
	-@(${DIFFNUM} -t \
		run/readamr_2d.ref output/readamr_2d.ref > readamr_2d.diff)
	ls -l readamr_2d.diff

test_readamr_3d:
	echo "test_readamr_3d_compile" > readamr_3d.diff
	./Config.pl -default -single -g=4,4,4 -r=2,2,2 -ng=0
	@(make READAMR)
	echo "test_readamr_3d_run" >> readamr_3d.diff
	(cd run; ${MPIRUN} ./READAMR.exe > readamr_3d.ref)
	-@(${DIFFNUM} -t -r=1e-5 \
		run/readamr_3d.ref output/readamr_3d.ref > readamr_3d.diff)
	ls -l readamr_3d.diff

test_readamr_sph:
	echo "test_readamr_sph_compile" > readamr_sph.diff
	./Config.pl -default -g=6,4,4 -r=2,2,2 -ng=0
	@(make READAMR)
	echo "test_readamr_sph_run" >> readamr_sph.diff
	(cd run; ${MPIRUN} ./READAMR.exe > readamr_sph.ref)
	-@(${DIFFNUM} -t \
		run/readamr_sph.ref output/readamr_sph.ref > readamr_sph.diff)
	ls -l readamr_sph.diff

test_readamr_sphgenr:
	echo "test_readamr_sphgenr_compile" > readamr_sphgenr.diff
	./Config.pl -default -g=6,4,4 -r=2,2,2 -ng=2
	@(make READAMR)
	echo "test_readamr_sphgrenr_run" >> readamr_sphgenr.diff
	(cd run; ${MPIRUN} ./READAMR.exe > readamr_sphgenr.ref)
	-@(${DIFFNUM} -t \
		run/readamr_sphgenr.ref output/readamr_sphgenr.ref \
		> readamr_sphgenr.diff)
	ls -l readamr_sphgenr.diff

test_readamr_c:
	echo "test_readamr_c_compile" > readamr_c.diff
	./Config.pl -default -nompi -single -g=4,4,4 -r=2,2,2 -ng=0
	@(make NOMPI)
	@(make READAMR_C)
	echo "test_readamr_c_run" >> readamr_c.diff
	(cd run; ./READAMR_C.exe > readamr_c.ref)
	-@(${DIFFNUM} -t -r=2e-6 -a=2e-5 \
		run/readamr_c.ref output/readamr_c.ref > readamr_c.diff)
	ls -l readamr_c.diff

test_readamr_py:
	echo "test_readamr_py_compile" > readamr_py.diff
	./Config.pl -default -nompi -single -g=4,4,4 -r=2,2,2 -ng=0
	@(make WRAPAMRLIB)
	echo "test_readamr_py_run" >> readamr_py.diff
	(cd srcReadAmr; python3 read_amr_py.py > readamr_py.ref;)
	-@(${DIFFNUM} -t -a=1e-11 -r=2e-6 \
		srcReadAmr/readamr_py.ref output/readamr_py.ref \
		> readamr_py.diff)
	ls -l readamr_py.diff

clean:
	cd share; make clean
	cd util/TIMING; make clean
	cd util/NOMPI/src; make clean 
	cd src; make clean
	cd srcReadAmr; make clean
	rm -f lib/*

allclean:
	cd src; make distclean
	cd srcReadAmr; make distclean
	rm -rf run lib bin *.diff
