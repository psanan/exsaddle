#####################
# exSaddle Makefile #
#####################

# This repository includes two custom preconditioners, based on ILUPACK.
# The versions here are included for historical reasons, but you should
# instead prefer to use the plugin versions, located at
#   https://bitbucket.org/psanan/pcilupack

# If you wish to build with PCILDL or PCILUPACK, define EXSADDLE_WITH_PCILDL or EXSADDLE_WITH_PCILUPACK, respectively, e.g.
# To be extra careful, we force you to set ILUPACK_PLATFORM
#  make EXSADDLE_WITH_PCILDL=1 ILUPACK_PLATFORM=GNU64

ifdef EXSADDLE_WITH_PCILUPACK
EXSADDLE_WITH_CUSTOMPC=1
endif
ifdef EXSADDLE_WITH_PCILDL
EXSADDLE_WITH_CUSTOMPC=1
ifdef EXSADDLE_WITH_PCILUPACK
# You cannot build with both, due to linking problems with ILUPACK
$(error you can define at most one of EXSADDLE_WITH_PCILDL and EXSADDLE_WITH_PCILUPACK)
endif
endif

ifdef EXSADDLE_WITH_CUSTOMPC
ifndef ILUPACK_PLATFORM
  $(error You must provide ILUPACK_PLATFORM, for example make ILUPACK_PLATFORM=GNU64, or set ILUPACK_PLATFORM in your environment. )
endif

ILUPACK_DIR = ilupack
ILUPACK_INCLUDE=-I${ILUPACK_DIR}/include
ILUPACK_LIBS=-L${ILUPACK_DIR}/lib/${ILUPACK_PLATFORM} -lilupack_mc64  -lmetisomp -lmetis -lmetisomp -lcamd -lamd -lsuitesparseconfig -lsparspak -llapack -lblaslike -lblas

# You need to build these by hand (see ilupack/notdistributed/README).
# This will likely involve commands like:
#     gfortran -O3 -fPIC mc64d.f -o MC64D.f
ILUPACK_LIBS+=${ILUPACK_DIR}/notdistributed/MC64D.o ${ILUPACK_DIR}/notdistributed/MC21D.o ${ILUPACK_DIR}/notdistributed/MC64S.o ${ILUPACK_DIR}/notdistributed/MC21S.o

endif

# --------------------------------------------------------------------------- #

EXECUTABLES=\
            ex42mod \
            ex23mod \
            exSaddle2d exSaddle2d_lame \
            exSaddle3d exSaddle3d_lame \

all : ${EXECUTABLES}

# --------------------------------------------------------------------------- #

EX42MOD_OBJ=ex42mod.o
EX23MOD_OBJ=ex23mod.o
EXSADDLE2D_OBJ=exSaddle2d.o femixedspace2d.o exSaddle_io2d.o models2d.o
EXSADDLE2D_LAME_OBJ=exSaddle2d_lame.o femixedspace2d_lame.o exSaddle_io2d_lame.o models2d_lame.o
EXSADDLE3D_OBJ=exSaddle3d.o femixedspace3d.o exSaddle_io3d.o models3d.o
EXSADDLE3D_LAME_OBJ=exSaddle3d_lame.o femixedspace3d_lame.o exSaddle_io3d_lame.o models3d_lame.o
ifdef EXSADDLE_WITH_PCILUPACK
EX42MOD_OBJ+=pcilupack.o
EX23MOD_OBJ+=pcilupack.o
EXSADDLE2D_OBJ+=pcilupack.o
EXSADDLE2D_LAME_OBJ+=pcilupack.o
EXSADDLE3D_OBJ+=pcilupack.o
EXSADDLE3D_LAME_OBJ+=pcilupack.o
CFLAGS+=-DEXSADDLE_WITH_PCILUPACK	 -fopenmp
endif
ifdef EXSADDLE_WITH_PCILDL
EX42MOD_OBJ+=pcildl.o
EX23MOD_OBJ+=pcildl.o
EXSADDLE2D_OBJ+=pcildl.o
EXSADDLE2D_LAME_OBJ+=pcildl.o
EXSADDLE3D_OBJ+=pcildl.o
EXSADDLE3D_LAME_OBJ+=pcildl.o
CFLAGS+=-DEXSADDLE_WITH_PCILDL -fopenmp
endif

# --------------------------------------------------------------------------- #

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# --------------------------------------------------------------------------- #

ex42mod: ${EX42MOD_OBJ}
	-${CLINKER} ${PCC_LINKER_FLAGS} ${CFLAGS} -o $@ $^ ${ILUPACK_LIBS} ${PETSC_KSP_LIB}

ex42mod.o : ex42mod.c
	-${CC} ${PCC_FLAGS} ${CFLAGS} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

ex23mod: ${EX23MOD_OBJ}
	-${CLINKER} ${PCC_LINKER_FLAGS} ${CFLAGS} -o $@ $^  ${ILUPACK_LIBS} ${PETSC_KSP_LIB}

ex23mod.o : ex23mod.c
	-${CC} ${PCC_FLAGS} ${CFLAGS} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

CFLAGS_2D=-DNSD=2
CFLAGS_3D=-DNSD=3
CFLAGS_LAME=-DLAME

exSaddle2d : ${EXSADDLE2D_OBJ}
	-${CLINKER} -o $@ $^ ${ILUPACK_LIBS} ${PETSC_KSP_LIB}

exSaddle2d.o : exSaddle.c exSaddle.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_2D}  -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

exSaddle3d : ${EXSADDLE3D_OBJ}
	-${CLINKER} -o $@ $^ ${ILUPACK_LIBS} ${PETSC_KSP_LIB}

exSaddle3d.o : exSaddle.c exSaddle.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_3D}  -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

exSaddle2d_lame: ${EXSADDLE2D_LAME_OBJ}
	-${CLINKER} -o $@ $^ ${ILUPACK_LIBS} ${PETSC_KSP_LIB}

exSaddle2d_lame.o : exSaddle.c exSaddle.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_2D} ${CFLAGS_LAME} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

exSaddle3d_lame: ${EXSADDLE3D_LAME_OBJ}
	-${CLINKER} -o $@ $^ ${ILUPACK_LIBS} ${PETSC_KSP_LIB}

exSaddle3d_lame.o : exSaddle.c exSaddle.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_3D} ${CFLAGS_LAME} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

# Note that we need different versions of the objects, since they will also depend on LAME and NSD
femixedspace2d.o : femixedspace.c femixedspace.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_2D} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

femixedspace2d_lame.o : femixedspace.c femixedspace.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_2D} ${CFLAGS_LAME} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

femixedspace3d.o : femixedspace.c femixedspace.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_3D} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

femixedspace3d_lame.o : femixedspace.c femixedspace.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_3D} ${CFLAGS_LAME} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

exSaddle_io2d.o : exSaddle_io.c exSaddle_io.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_2D} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

exSaddle_io2d_lame.o : exSaddle_io.c exSaddle_io.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_2D} ${CFLAGS_LAME} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

exSaddle_io3d.o : exSaddle_io.c exSaddle_io.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_3D} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

exSaddle_io3d_lame.o : exSaddle_io.c exSaddle_io.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_3D} ${CFLAGS_LAME} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

models2d.o : models.c models.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_2D} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

models2d_lame.o : models.c models.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_2D} ${CFLAGS_LAME} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

models3d.o : models.c models.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_3D} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

models3d_lame.o : models.c models.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${CFLAGS_3D} ${CFLAGS_LAME} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -c -o $@ $<

# --------------------------------------------------------------------------- #

pcildl.o : pcildl.c pcildl.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} -c ${ILUPACK_INCLUDE} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -o $@ $<

pcilupack.o : pcilupack.c pcilupack.h
	-${CC} ${PCC_FLAGS} ${CFLAGS} ${ILUPACK_INCLUDE} -c ${ILUPACK_INCLUDE} -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -o $@ $<

# --------------------------------------------------------------------------- #
# Tests
#
# * These tests are mainly intended as regression tests. The solvers here
#   can be quite sensitive, and the examples below are NOT intended to
#   be examples of efficient, or even numerically stable, solver choices.
# * Some of these tests require UMFPACK. (Configure PETSc --download-suiteparse)
#   PETSc's built-in LU does NOT work stably as a coarse grid solver for saddle point systems
# * Do NOT configure PETSc with METIS, or things may fail
# * Subdomain-based solves (e.g. ASM) rely on a custom (element-based) definition
#   of subdomains which has not been thoroughly tested and which should be considered
#   very experimental
# * The diff-based approach here is obviously suboptimal and non-portable. It
#   would be better to define tests with a tool like SciATH (github.com/sciath/sciath)

# define this  to 1 to copy all test output to the corresponding reference file (Careful!)
COPY_TEST_OUTPUT=0
# e.g. to update all test output (Careful!) you could do
#   make COPY_TEST_OUTPUT=1 test
# This will produce the diffs, and then copy.

test : test_stokes test_lame
test_stokes : test_stokes2d test_stokes3d
test_lame : test_lame2d test_lame3d
test_stokes2d : \
  test_exSaddle2d_1         test_exSaddle2d_2 \
  test_exSaddle2d_mg_1      test_exSaddle2d_mg_2 \
  test_exSaddle2d_fs_1      test_exSaddle2d_fs_2  \
	test_exSaddle2d_asm_1 \
	test_exSaddle2d_mms_1 \
	test_exSaddle2d_ar_1
test_lame2d : \
  test_exSaddle2d_lame_1    test_exSaddle2d_lame_2 \
  test_exSaddle2d_lame_mg_1 test_exSaddle2d_lame_mg_2 \
  test_exSaddle2d_lame_fs_1 test_exSaddle2d_lame_fs_2
TESTS_STOKES3D = \
  test_exSaddle3d_1         test_exSaddle3d_2 \
  test_exSaddle3d_mg_1      test_exSaddle3d_mg_2 \
  test_exSaddle3d_fs_1      test_exSaddle3d_fs_2 \
  test_exSaddle3d_mg_fs_coarse_1   \
	test_exSaddle3d_asm_1  \
	test_exSaddle3d_mg_asm_1 \
	test_exSaddle3d_ar_1 \
	test_exSaddle3d_pseudoice_1
ifdef EXSADDLE_WITH_PCILDL
TESTS_STOKES3D += test_exSaddle3d_ildl_1
endif
ifdef EXSADDLE_WITH_PCILUPACK
TESTS_STOKES3D += test_exSaddle3d_ilupack_1
endif
test_stokes3d : ${TESTS_STOKES3D}
test_lame3d : \
  test_exSaddle3d_lame_1    test_exSaddle3d_lame_2 \
	test_exSaddle3d_lame_3    test_exSaddle3d_lame_4 \
	test_exSaddle3d_lame_5 \
  test_exSaddle3d_lame_mg_1 test_exSaddle3d_lame_mg_2 \
  test_exSaddle3d_lame_fs_1 test_exSaddle3d_lame_fs_2 \

.PHONY : \
  test test_stokes test_lame test_stokes2d test_stokes3d test_lame2d test_lame3d \
  test_exSaddle2d_1         test_exSaddle2d_2 \
  test_exSaddle2d_mg_1      test_exSaddle2d_mg_2 \
  test_exSaddle2d_fs_1      test_exSaddle2d_fs_2 \
	test_exSaddle2d_asm_1 \
	test_exSaddle2d_mms_1 \
	test_exSaddle2d_ar_1 \
  test_exSaddle2d_lame_1    test_exSaddle2d_lame_2 \
  test_exSaddle2d_lame_mg_1 test_exSaddle2d_lame_mg_2 \
  test_exSaddle2d_lame_fs_1 test_exSaddle2d_lame_fs_2 \
  test_exSaddle3d_1         test_exSaddle3d_2 \
  test_exSaddle3d_mg_1      test_exSaddle3d_mg_2 \
  test_exSaddle3d_fs_1      test_exSaddle3d_fs_2 \
  test_exSaddle3d_mg_fs_coarse_1 \
	test_exSaddle3d_asm_1 \
	test_exSaddle3d_mg_asm_1 \
	test_exSaddle3d_ar_1 \
  test_exSaddle3d_ildl_1 \
  test_exSaddle3d_lame_1    test_exSaddle3d_lame_2 \
	test_exSaddle3d_lame_3    test_exSaddle3d_lame_4 \
	test_exSaddle3d_lame_5 \
  test_exSaddle3d_lame_mg_1 test_exSaddle3d_lame_mg_2 \
  test_exSaddle3d_lame_fs_1 test_exSaddle3d_lame_fs_2 \
	test_exSaddle3d_pseudoice_1 \

test_exSaddle2d_1 :
	-@${MPIEXEC} -n 1 ./exSaddle2d  -model 0 -mx 4 -diagnostics -saddle_ksp_max_it 100 -saddle_ksp_converged_reason -saddle_pc_type jacobi > exSaddle2d_1.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle2d_1.ref exSaddle2d_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_1.tmp testref/exSaddle2d_1.ref; fi; \
	   ${RM} -f exSaddle2d_1.tmp

test_exSaddle2d_2 :
	-@${MPIEXEC} -n 2 ./exSaddle2d -model 0 -mx 4 -diagnostics -saddle_ksp_max_it 100 -saddle_ksp_converged_reason -saddle_pc_type jacobi > exSaddle2d_2.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle2d_2.ref exSaddle2d_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_2, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_2.tmp testref/exSaddle2d_2.ref; fi; \
	   ${RM} -f exSaddle2d_2.tmp

test_exSaddle2d_mg_1 :
	-@${MPIEXEC} -n 1  ./exSaddle2d -model 0 -mx 16 -mg -nlevels 3 -diagnostics -saddle_ksp_type fgmres -saddle_mg_levels_ksp_type gmres -saddle_mg_levels_pc_type jacobi -saddle_mg_levels_ksp_max_it 10 -saddle_ksp_monitor_short -saddle_mg_coarse_pc_factor_mat_solver_type umfpack > exSaddle2d_mg_1.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle2d_mg_1.ref exSaddle2d_mg_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_mg_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_mg_1.tmp testref/exSaddle2d_mg_1.ref; fi; \
	   ${RM} -f exSaddle2d_mg_1.tmp

test_exSaddle2d_mg_2 :
	-@${MPIEXEC} -n 2 ./exSaddle2d  -model 0 -mx 16 -mg -nlevels 3 -diagnostics -saddle_ksp_type fgmres -saddle_mg_levels_ksp_type gmres -saddle_mg_levels_pc_type jacobi -saddle_mg_levels_ksp_max_it 10 -saddle_ksp_monitor_short -saddle_mg_coarse_redundant_pc_factor_mat_solver_type umfpack > exSaddle2d_mg_2.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle2d_mg_2.ref exSaddle2d_mg_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_mg_2, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_mg_2.tmp testref/exSaddle2d_mg_2.ref; fi; \
	   ${RM} -f exSaddle2d_mg_2.tmp

test_exSaddle2d_fs_1 :
	-@${MPIEXEC} -n 1 ./exSaddle2d  -model 0 -fs -mx 6 -diagnostics -saddle_ksp_monitor_short > exSaddle2d_fs_1.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle2d_fs_1.ref exSaddle2d_fs_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_fs_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_fs_1.tmp testref/exSaddle2d_fs_1.ref; fi; \
	   ${RM} -f exSaddle2d_fs_1.tmp

test_exSaddle2d_fs_2 :
	-@${MPIEXEC} -n 2 ./exSaddle2d -model 0 -fs -mx 6 -diagnostics -saddle_ksp_monitor_short > exSaddle2d_fs_2.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle2d_fs_2.ref exSaddle2d_fs_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_fs_2, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_fs_2.tmp testref/exSaddle2d_fs_2.ref; fi; \
	   ${RM} -f exSaddle2d_fs_2.tmp

# This is a very unstable solve, only retained here for regression testing purposes:
test_exSaddle2d_asm_1 :
	-@${MPIEXEC} -n 9 ./exSaddle2d -mx 12 -saddle_pc_type asm -saddle_pc_asm_dm_subdomains -set_ksp_dm -options_left -saddle_ksp_monitor_short -saddle_sub_ksp_type preonly -saddle_sub_pc_type lu -saddle_sub_pc_factor_mat_solver_type umfpack -dmdafe_overlap 1 -saddle_ksp_rtol 1e-4 > exSaddle2d_asm_1.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle2d_asm_1.ref exSaddle2d_asm_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_asm_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_asm_1.tmp testref/exSaddle2d_asm_1.ref; fi; \
	   ${RM} -f exSaddle2d_asm_1.tmp

test_exSaddle2d_mms_1 :
	-@${MPIEXEC} ./exSaddle2d -saddle_pc_type lu -saddle_pc_factor_mat_solver_type umfpack -model 101 -check_solution -saddle_ksp_monitor_short  -mx 16 -constant_pressure_nullspace > exSaddle2d_mms_1.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle2d_mms_1.ref exSaddle2d_mms_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_mms_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_mms_1.tmp testref/exSaddle2d_mms_1.ref; fi; \
	   ${RM} -f exSaddle2d_mms_1.tmp

test_exSaddle2d_ar_1 :
	-@${MPIEXEC} ./exSaddle2d -options_file abf.opts  -saddle_ksp_monitor_short -model 0 -mx 32 -my 32 -options_left -size_y 0.1 > exSaddle2d_ar_1.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle2d_ar_1.ref exSaddle2d_ar_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_ar_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_ar_1.tmp testref/exSaddle2d_ar_1.ref; fi; \
	   ${RM} -f exSaddle2d_ar_1.tmp


test_exSaddle2d_lame_1 :
	-@${MPIEXEC} -n 1 ./exSaddle2d_lame -model 6 -saddle_pc_type jacobi -saddle_ksp_converged_reason -mx 8 -diagnostics > exSaddle2d_lame_1.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle2d_lame_1.ref exSaddle2d_lame_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_lame_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_lame_1.tmp testref/exSaddle2d_lame_1.ref; fi; \
	   ${RM} -f exSaddle2d_lame_1.tmp

test_exSaddle2d_lame_2 :
	-@${MPIEXEC} -n 2 ./exSaddle2d_lame -model 6 -saddle_pc_type jacobi -saddle_ksp_converged_reason -mx 8 -diagnostics > exSaddle2d_lame_2.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle2d_lame_2.ref exSaddle2d_lame_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_lame_2, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_lame_2.tmp testref/exSaddle2d_lame_2.ref; fi; \
	   ${RM} -f exSaddle2d_lame_2.tmp

test_exSaddle2d_lame_mg_1 :
	-@${MPIEXEC} -n 1 ./exSaddle2d_lame  -mx 16 -mg -nlevels 3 -diagnostics -saddle_ksp_type fgmres -saddle_mg_levels_ksp_type gmres -saddle_mg_levels_pc_type jacobi -saddle_mg_levels_ksp_max_it 10 -saddle_ksp_monitor_short -saddle_mg_coarse_pc_factor_mat_solver_type umfpack > exSaddle2d_lame_mg_1.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle2d_lame_mg_1.ref exSaddle2d_lame_mg_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_lame_mg_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_lame_mg_1.tmp testref/exSaddle2d_lame_mg_1.ref; fi; \
	   ${RM} -f exSaddle2d_lame_mg_1.tmp

test_exSaddle2d_lame_mg_2 :
	-@${MPIEXEC} -n 2 ./exSaddle2d_lame  -mx 16 -mg -nlevels 3 -diagnostics -saddle_ksp_type fgmres -saddle_mg_levels_ksp_type gmres -saddle_mg_levels_pc_type jacobi -saddle_mg_levels_ksp_max_it 10 -saddle_ksp_monitor_short -saddle_mg_coarse_redundant_pc_factor_mat_solver_type umfpack > exSaddle2d_lame_mg_2.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle2d_lame_mg_2.ref exSaddle2d_lame_mg_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_lame_mg_2, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_lame_mg_2.tmp testref/exSaddle2d_lame_mg_2.ref; fi; \
	   ${RM} -f exSaddle2d_lame_mg_2.tmp

test_exSaddle2d_lame_fs_1 :
	-@${MPIEXEC} -n 1 ./exSaddle2d_lame  -model 6 -fs -mx 6 -diagnostics -saddle_ksp_monitor_short -saddle_ksp_converged_reason > exSaddle2d_lame_fs_1.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle2d_lame_fs_1.ref exSaddle2d_lame_fs_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_lame_fs_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_lame_fs_1.tmp testref/exSaddle2d_lame_fs_1.ref; fi; \
	   ${RM} -f exSaddle2d_lame_fs_1.tmp

test_exSaddle2d_lame_fs_2 :
	-@${MPIEXEC} -n 2 ./exSaddle2d_lame -model 6 -fs -mx 6 -diagnostics -saddle_ksp_monitor_short -saddle_ksp_converged_reason > exSaddle2d_lame_fs_2.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle2d_lame_fs_2.ref exSaddle2d_lame_fs_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle2d_lame_fs_2, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle2d_lame_fs_2.tmp testref/exSaddle2d_lame_fs_2.ref; fi; \
	   ${RM} -f exSaddle2d_lame_fs_2.tmp

test_exSaddle3d_1 :
	-@${MPIEXEC} -n 1 ./exSaddle3d -saddle_pc_type jacobi -diagnostics -model 1 -saddle_ksp_converged_reason -mx 4 -my 7 -mz 5 -saddle_ksp_max_it 10 > exSaddle3d_1.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle3d_1.ref exSaddle3d_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_1.tmp testref/exSaddle3d_1.ref; fi; \
	   ${RM} -f exSaddle3d_1.tmp

test_exSaddle3d_2 :
	-@${MPIEXEC} -n 2 ./exSaddle3d -saddle_pc_type jacobi -diagnostics -model 1 -saddle_ksp_converged_reason -mx 4 -my 7 -mz 5 -saddle_ksp_max_it 10 > exSaddle3d_2.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle3d_2.ref exSaddle3d_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_2, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_2.tmp testref/exSaddle3d_2.ref; fi; \
	   ${RM} -f exSaddle3d_2.tmp

test_exSaddle3d_mg_1 :
	-@${MPIEXEC} -n 1 ./exSaddle3d -model 2 -sinker_n 1 -mx 8 -mg -nlevels 2 -diagnostics -saddle_ksp_type fgmres -saddle_mg_levels_ksp_type gmres -saddle_mg_levels_pc_type jacobi -saddle_mg_levels_ksp_max_it 10 -saddle_ksp_monitor_short -saddle_mg_coarse_pc_factor_mat_solver_type umfpack > exSaddle3d_mg_1.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle3d_mg_1.ref exSaddle3d_mg_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_mg_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_mg_1.tmp testref/exSaddle3d_mg_1.ref; fi; \
	   ${RM} -f exSaddle3d_mg_1.tmp

test_exSaddle3d_mg_2 :
	-@${MPIEXEC} -n 2 ./exSaddle3d -model 2 -sinker_n 1 -mx 8 -mg -nlevels 2 -diagnostics -saddle_ksp_type fgmres -saddle_mg_levels_ksp_type gmres -saddle_mg_levels_pc_type jacobi -saddle_mg_levels_ksp_max_it 10 -saddle_ksp_monitor_short -saddle_mg_coarse_redundant_pc_factor_mat_solver_type umfpack > exSaddle3d_mg_2.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle3d_mg_2.ref exSaddle3d_mg_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_mg_2, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_mg_2.tmp testref/exSaddle3d_mg_2.ref; fi; \
	   ${RM} -f exSaddle3d_mg_2.tmp

test_exSaddle3d_mg_fs_coarse_1 :
	-@${MPIEXEC} -n 1 ./exSaddle3d -saddle_ksp_view -mg -nlevels 2 -fs_coarse -saddle_mg_levels_ksp_type gmres -saddle_ksp_monitor_short -saddle_ksp_converged_reason -saddle_mg_coarse_fieldsplit_u_pc_type jacobi -saddle_mg_coarse_fieldsplit_p_pc_type jacobi -saddle_mg_coarse_ksp_type fgmres -saddle_mg_coarse_ksp_convergence_test default  -saddle_mg_levels_pc_type jacobi > exSaddle3d_mg_fs_coarse_1.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle3d_mg_fs_coarse_1.ref exSaddle3d_mg_fs_coarse_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_mg_fs_coarse_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_mg_fs_coarse_1.tmp testref/exSaddle3d_mg_fs_coarse_1.ref; fi; \
	   ${RM} -f exSaddle3d_mg_fs_coarse_1.tmp

test_exSaddle3d_fs_1 :
	-@${MPIEXEC} -n 1 ./exSaddle3d -model 2 -sinker_n 1 -fs -mx 4 -diagnostics -saddle_ksp_monitor_short > exSaddle3d_fs_1.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle3d_fs_1.ref exSaddle3d_fs_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_fs_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_fs_1.tmp testref/exSaddle3d_fs_1.ref; fi; \
	   ${RM} -f exSaddle3d_fs_1.tmp

test_exSaddle3d_fs_2 :
	-@${MPIEXEC} -n 2 ./exSaddle3d -model 2 -sinker_n 1 -fs -mx 4 -diagnostics -saddle_ksp_monitor_short > exSaddle3d_fs_2.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle3d_fs_2.ref exSaddle3d_fs_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_fs_2, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_fs_2.tmp testref/exSaddle3d_fs_2.ref; fi; \
	   ${RM} -f exSaddle3d_fs_2.tmp

test_exSaddle3d_asm_1 :
	-@${MPIEXEC} -n 8 ./exSaddle3d -mx 6 -saddle_pc_type asm -saddle_pc_asm_dm_subdomains -set_ksp_dm -options_left -saddle_ksp_monitor_short -saddle_sub_ksp_type preonly -saddle_sub_pc_type lu -saddle_sub_pc_factor_mat_solver_type umfpack> exSaddle3d_asm_1.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle3d_asm_1.ref exSaddle3d_asm_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_asm_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_asm_1.tmp testref/exSaddle3d_asm_1.ref; fi; \
	   ${RM} -f exSaddle3d_asm_1.tmp

test_exSaddle3d_mg_asm_1 :
	-@${MPIEXEC} -n 4 ./exSaddle3d -options_left -mg -nlevels 2 -saddle_mg_levels_ksp_type gmres -saddle_mg_levels_pc_type asm -saddle_mg_levels_pc_asm_dm_subdomains -dmdafe_overlap 1 -saddle_ksp_monitor_short -saddle_ksp_pc_side right -saddle_mg_coarse_redundant_pc_factor_mat_solver_type umfpack -saddle_mg_levels_sub_pc_type lu -saddle_mg_levels_sub_pc_factor_mat_solver_type umfpack -mx 6 -my 4 -mz 4 > exSaddle3d_mg_asm_1.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle3d_mg_asm_1.ref exSaddle3d_mg_asm_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_mg_asm_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_mg_asm_1.tmp testref/exSaddle3d_mg_asm_1.ref; fi; \
	   ${RM} -f exSaddle3d_mg_asm_1.tmp

test_exSaddle3d_ar_1 :
	-@${MPIEXEC} ./exSaddle3d -options_file abf.opts  -saddle_ksp_monitor_short -model 0 -mx 6 -my 6 -mz 6 -options_left  -saddle_fieldsplit_u_ksp_converged_reason -size_z 0.1 > exSaddle3d_ar_1.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle3d_ar_1.ref exSaddle3d_ar_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_ar_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_ar_1.tmp testref/exSaddle3d_ar_1.ref; fi; \
	   ${RM} -f exSaddle3d_ar_1.tmp

test_exSaddle3d_lame_1 :
	-@${MPIEXEC} -n 1 ./exSaddle3d_lame -model 6 -saddle_pc_type jacobi -saddle_ksp_converged_reason -mx 4 -diagnostics > exSaddle3d_lame_1.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle3d_lame_1.ref exSaddle3d_lame_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_lame_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_lame_1.tmp testref/exSaddle3d_lame_1.ref; fi; \
	   ${RM} -f exSaddle3d_lame_1.tmp

test_exSaddle3d_lame_2 :
	-@${MPIEXEC} -n 2 ./exSaddle3d_lame -model 6 -saddle_pc_type jacobi -saddle_ksp_converged_reason -mx 4 -diagnostics > exSaddle3d_lame_2.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle3d_lame_2.ref exSaddle3d_lame_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_lame_2, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_lame_2.tmp testref/exSaddle3d_lame_2.ref; fi; \
	   ${RM} -f exSaddle3d_lame_2.tmp

test_exSaddle3d_lame_3 :
	-@${MPIEXEC} -n 1 ./exSaddle3d_lame -model 2 -lambda1 10 -mu1 100 -saddle_ksp_monitor_short -saddle_pc_type jacobi -saddle_ksp_max_it 10 -saddle_ksp_type gmres -saddle_ksp_pc_side right -mx 4 -diagnostics > exSaddle3d_lame_3.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle3d_lame_3.ref exSaddle3d_lame_3.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_lame_3, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_lame_3.tmp testref/exSaddle3d_lame_3.ref; fi; \
	   ${RM} -f exSaddle3d_lame_3.tmp

test_exSaddle3d_lame_4 :
	-@${MPIEXEC} -n 1 ./exSaddle3d_lame -model 8 -lambda1 10 -lambda0 10 -saddle_ksp_monitor_short -saddle_pc_type jacobi -saddle_ksp_max_it 10 -saddle_ksp_type gmres -saddle_ksp_pc_side right -mx 4 -diagnostics > exSaddle3d_lame_4.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle3d_lame_4.ref exSaddle3d_lame_4.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_lame_4, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_lame_4.tmp testref/exSaddle3d_lame_4.ref; fi; \
	   ${RM} -f exSaddle3d_lame_4.tmp

test_exSaddle3d_lame_5 :
	-@${MPIEXEC} -n 1 ./exSaddle3d_lame -model 9 -saddle_ksp_monitor_short -saddle_pc_type jacobi -saddle_ksp_max_it 10 -saddle_ksp_type gmres -saddle_ksp_pc_side right -mx 4 -diagnostics > exSaddle3d_lame_5.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle3d_lame_5.ref exSaddle3d_lame_5.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_lame_5, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_lame_5.tmp testref/exSaddle3d_lame_5.ref; fi; \
	   ${RM} -f exSaddle3d_lame_5.tmp

test_exSaddle3d_lame_mg_1 :
	-@${MPIEXEC} -n 1 ./exSaddle3d_lame -model 6 -mx 6 -mg -nlevels 2 -diagnostics -saddle_ksp_type fgmres -saddle_mg_levels_ksp_type gmres -saddle_mg_levels_pc_type jacobi -saddle_mg_levels_ksp_max_it 10 -saddle_ksp_monitor_short -saddle_mg_coarse_pc_factor_mat_solver_type umfpack > exSaddle3d_lame_mg_1.tmp 2>&1;\
	   if (${DIFF} testref/exSaddle3d_lame_mg_1.ref exSaddle3d_lame_mg_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_lame_mg_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_lame_mg_1.tmp testref/exSaddle3d_lame_mg_1.ref; fi; \
	   ${RM} -f exSaddle3d_lame_mg_1.tmp

test_exSaddle3d_lame_mg_2 :
	-@${MPIEXEC} -n 2 ./exSaddle3d_lame -model 6 -mx 6 -mg -nlevels 2 -diagnostics -saddle_ksp_type fgmres -saddle_mg_levels_ksp_type gmres -saddle_mg_levels_pc_type jacobi -saddle_mg_levels_ksp_max_it 10 -saddle_ksp_monitor_short -saddle_mg_coarse_redundant_pc_factor_mat_solver_type umfpack > exSaddle3d_lame_mg_2.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle3d_lame_mg_2.ref exSaddle3d_lame_mg_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_lame_mg_2, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_lame_mg_2.tmp testref/exSaddle3d_lame_mg_2.ref; fi; \
	   ${RM} -f exSaddle3d_lame_mg_2.tmp

test_exSaddle3d_lame_fs_1 :
	-@${MPIEXEC} -n 1 ./exSaddle3d_lame -model 6 -fs -mx 4 -diagnostics -saddle_ksp_monitor_short -saddle_fieldsplit_u_ksp_max_it 10 -saddle_fieldsplit_p_ksp_type preonly -saddle_ksp_max_it 10 > exSaddle3d_lame_fs_1.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle3d_lame_fs_1.ref exSaddle3d_lame_fs_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_lame_fs_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_lame_fs_1.tmp testref/exSaddle3d_lame_fs_1.ref; fi; \
	   ${RM} -f exSaddle3d_lame_fs_1.tmp

test_exSaddle3d_lame_fs_2 :
	-@${MPIEXEC} -n 2 ./exSaddle3d_lame -model 6 -fs -mx 4 -diagnostics -saddle_ksp_monitor_short -saddle_fieldsplit_u_ksp_max_it 10 -saddle_fieldsplit_p_ksp_type preonly -saddle_ksp_max_it 10 > exSaddle3d_lame_fs_2.tmp 2>&1;	  \
	   if (${DIFF} testref/exSaddle3d_lame_fs_2.ref exSaddle3d_lame_fs_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_lame_fs_2, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_lame_fs_2.tmp testref/exSaddle3d_lame_fs_2.ref; fi; \
	   ${RM} -f exSaddle3d_lame_fs_2.tmp

test_exSaddle3d_ildl_1 :
	-@${MPIEXEC} -n 1 ./exSaddle3d -diagnostics -mx 8 -model 6 -eta1 100 -eta0 1 -saddle_ksp_monitor_short -saddle_pc_type ildl -options_left -saddle_pc_ildl_droptol 1e-3 -saddle_ksp_pc_side right > exSaddle3d_ildl_1.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle3d_ildl_1.ref exSaddle3d_ildl_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_ildl_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_ildl_1.tmp testref/exSaddle3d_ildl_1.ref; fi; \
	   ${RM} -f exSaddle3d_ildl_1.tmp

test_exSaddle3d_ilupack_1 :
	-@${MPIEXEC} -n 1 ./exSaddle3d -saddle_ksp_view -saddle_pc_type ilupack -saddle_pc_ilupack_droptol 1e-3 -saddle_pc_ilupack_condest 100 -saddle_pc_ilupack_droptolS 1e-4 -mx 4 -saddle_ksp_monitor_short > exSaddle3d_ilupack_1.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle3d_ilupack_1.ref exSaddle3d_ilupack_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_ilupack_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_ilupack_1.tmp testref/exSaddle3d_ilupack_1.ref; fi; \
	   ${RM} -f exSaddle3d_ilupack_1.tmp

test_exSaddle3d_pseudoice_1 :
	-@${MPIEXEC} -n 1 ./exSaddle3d -saddle_ksp_view -options_file abf.opts -model 11 -size_x 0.1 -mx 6 -saddle_ksp_monitor_short > exSaddle3d_pseudoice_1.tmp 2>&1; \
	   if (${DIFF} testref/exSaddle3d_pseudoice_1.ref exSaddle3d_pseudoice_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with exSaddle3d_pseudoice_1, diffs above\n=========================================\n"; fi; \
		 if [ ${COPY_TEST_OUTPUT} -eq 1 ] ; then cp exSaddle3d_pseudoice_1.tmp testref/exSaddle3d_pseudoice_1.ref; fi; \
	   ${RM} -f exSaddle3d_pseudoice_1.tmp


# --------------------------------------------------------------------------- #

clean ::
	rm -f *.o ${EXECUTABLES}

clean_output :
	rm -f *.vts *.pvts *.petscbin *.petscbin.info *.gp

.PHONY : clean clean_output
