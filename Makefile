#
# Makefile for selecting parallel/ensemble average/single CPU version
# of NBODY6++  June 2003 R.Sp.
#
# Usage notes: these examples show how to adapt the Makefile targets to
# specific computer architectures or compilers. The targets defined here
# should cover the most important cases and can be used as templates for
# other machines and compilers.
#
# make, make standard, make nbody6         (single CPU version)
# make pgf                         (single CPU pgf77-Comp. generic)
# make pgfp3                       (single CPU pgf77-Comp. Pentium 3 optimised)
# make pgfp4                       (single CPU pgf77-Comp. Pentium 4 optimised)
# make mpich                       (Beowulf PC clusters, standard mpif77)
# make mpichpgf                    (Beowulf PC clusters, using pgf77)
# make sun                         (single CPU sun-Comp. generic)
# make sun_e10k			   (parallel Sun Enterprise)
# make sun_e10k_fast               (parallel Sun Enterprise optimised)
# make crayt3e NPE=n              (fully parallel version using SHMEM CRAY T3E)
# make crayt3e_noshmem NPE=n      (fully parallel version no SHMEM CRAY T3E)
# make jump                       (fully parallel version on IBM Jump with AIX)
#
# make mpich_ensemble, mpichpgf_ensemble, sun_e10k_ensemble, 
#      sun_e10k_fast_ensemble, crayt3e_ensemble
#                                 (parallel run, but for ensemble job farming)
# make visit, pgfvisit, mpichvisit: (including the visit online visualisation
#                                     library, experimental, need AVS license)
#
# make clean                       (remove all .o and .f which should not exist)
# make cleanpar                    (remove all .o originating from .F)
#
#
LOAD = $(FC)
CC = cc
RESULT = nbody6
LFLAGS = -o $(RESULT)
LFLAGS2 =
CUDA_PATH=/usr/local/cuda
SDK_PATH=/usr/local/cuda_sdk/C
VISITFLAGS = -DVISIT
VISITLD = -L ./lvisit -llvisit_nbody -L/work/Tuc/spurzem/visit/visit/lvisit/lib -llvisit -L/work/Tuc/spurzem/visit/visit/visit20/lib -lvisit

#----------------------------------------------------------------------------
# Workstation single CPU g77 Linux -O6
#----------------------------------------------------------------------------
standard: 
	cp mpif.null.h mpif.rainer.h
	$(MAKE) $(RESULT) "FC = gfortran" "FFLAGS = -O3" "SOURCE = energy.f $(SOURCE)"
gpu:
	cp mpif.null.h mpif.rainer.h
	$(MAKE) nbody6.gpu "FC = gfortran" "LFLAGS = -o nbody6.gpu -L$(CUDA_PATH)/lib64 -lcudart" \
		"FFLAGS = -O3 -g -fbounds-check -D GPU" \
		"SOURCE = energy.f util_gpu.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation single CPU pgf77 Linux generic
#----------------------------------------------------------------------------
pgf:
	cp mpif.null.h mpif.rainer.h
	$(MAKE) $(RESULT) "FC = pgf77" "FFLAGS = -O4" \
		"SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation single CPU pgf77 Linux (PIII, PII, PPro)
#----------------------------------------------------------------------------
pgfp3:
	cp mpif.null.h mpif.rainer.h
	$(MAKE) $(RESULT) "FC = pgf77" "FFLAGS = -O4 -tp p6 -pc 64 -fast" \
		"SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation single CPU pgf77 Linux (P4)
#----------------------------------------------------------------------------
pgfp4:
	cp mpif.null.h mpif.rainer.h
	$(MAKE) $(RESULT) "FC = pgf77" "FFLAGS = -O4 -tp p7 -pc 64 -fast -fastsse" \
		"SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation parallel CPU f77 Linux (Beowulf generic mpif77)
#----------------------------------------------------------------------------
mpich:
	cp mpif.mpich.h mpif.rainer.h
	$(MAKE) $(RESULT) "FC = mpif77" "FFLAGS = -O4 -D PUREMPI -D PARALLEL" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
mpich_ensemble:
	cp mpif.mpich.h mpif.rainer.h
	$(MAKE) $(RESULT) "FC = mpif77" "FFLAGS = -O4 -D ENSEMBLE" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
mpichgpu:
	cp mpif.mpich.h mpif.rainer.h
	$(MAKE) nbody6.gpu "FC = mpif77" "LFLAGS = -o nbody6.gpu -L$(CUDA_PATH)/lib64 -lcudart" \
		"FFLAGS = -O4 -D PUREMPI -D PARALLEL -D GPU" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f util_gpu.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation parallel CPU pgf77 Linux (Beowulf pgf77 P4 optimised)
#----------------------------------------------------------------------------
mpichpgf:
	cp mpif.mpich.h mpif.rainer.h
	$(MAKE) $(RESULT) "LOAD = mpif77-pgi" "FC = pgf77" "FFLAGS = -O4 -tp p7 -pc 64 -fast -fastsse -D PUREMPI -D PARALLEL" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
mpichpgf_ensemble:
	cp mpif.mpich.h mpif.rainer.h
	$(MAKE) $(RESULT) "LOAD = mpif77-pgi" "FC = pgf77" "FFLAGS = -O4 -tp p7 -pc 64 -fast -fastsse -D ENSEMBLE" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation single CPU f77 Solaris 
#----------------------------------------------------------------------------
sun:
	cp mpif.null.h mpif.rainer.h
	$(MAKE) $(RESULT) "FC = f77" "FFLAGS = -O4 " "SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Parallel Sun
#----------------------------------------------------------------------------
sun_e10k:
	cp mpif.sun.h mpif.rainer.h
	$(MAKE) $(RESULT) "FFLAGS = -D PUREMPI -D PARALLEL" \
		"FC = tmf77" "LFLAGS = -lmpi $(LFLAGS)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
sun_e10k_ensemble:
	cp mpif.sun.h mpif.rainer.h
	$(MAKE) $(RESULT) "FFLAGS = -D ENSEMBLE" \
		"FC = tmf77" "LFLAGS = -lmpi $(LFLAGS)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
sun_e10k_fast:
	cp mpif.sun.h mpif.rainer.h
	$(MAKE) $(RESULT) "FFLAGS = -fast -xarch=v9a -D PUREMPI -D PARALLEL" \
		"FC = tmf77" "LFLAGS = -lmpi $(LFLAGS)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
sun_e10k_fast_ensemble:
	cp mpif.sun.h mpif.rainer.h
	$(MAKE) $(RESULT) "FFLAGS = -fast -xarch=v9a -D ENSEMBLE" \
		"FC = tmf77" "LFLAGS = -lmpi $(LFLAGS)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
#----------------------------------------------------------------------------
# IBM JUMP Fully parallel Mode   (-bnoquiet)
#----------------------------------------------------------------------------
jump:
	sed 's/APPEND/SEQUENTIAL/g' file_init.F > x.F
	mv x.F file_init.F
	cp flush.t3e.f flush.f
	cp mpif.jump.h mpif.rainer.h
	$(MAKE) $(RESULT) "FFLAGS = -I . -q64 -qfixed -O3 -qstrict -qsave" \
		"FC = mpxlf_r" "CPPFLAGS = -DPUREMPI -DPARALLEL" \
		"SOURCE = energy.f energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f flush.f $(SOURCE)" \
		"LFLAGS = -bmaxdata:500000000 $(LFLAGS)"
#----------------------------------------------------------------------------
# CRAY T3E Fully parallel Mode including SHMEM
#----------------------------------------------------------------------------
crayt3e:
	sed 's/APPEND/SEQUENTIAL/g' file_init.F > x.F
	mv x.F file_init.F
	cp mpif.t3e.h mpif.rainer.h
	cp flush.t3e.f flush.f
	$(MAKE) $(RESULT) "FFLAGS = -d p -D SHMEM -D PARALLEL" \
		"FC = f90" "LFLAGS = -X$(NPE) $(LFLAGS)-$(NPE)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f flush.f $(SOURCE)"
crayt3e_ensemble:
	sed 's/APPEND/SEQUENTIAL/g' file_init.F > x.F
	mv x.F file_init.F
	cp mpif.t3e.h mpif.rainer.h
	cp flush.t3e.f flush.f
	$(MAKE) $(RESULT) "FFLAGS = -d p -D ENSEMBLE" \
		"FC = f90" "LFLAGS = -X$(NPE) $(LFLAGS)-$(NPE)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f flush.f $(SOURCE)"
crayt3e_noshmem:
	sed 's/APPEND/SEQUENTIAL/g' file_init.F > x.F
	mv x.F file_init.F
	cp mpif.t3e.h mpif.rainer.h
	cp flush.t3e.f flush.f
	$(MAKE) $(RESULT) "FFLAGS = -d p -D PUREMPI -D PARALLEL" \
		"FC = f90" "LFLAGS = -X$(NPE) $(LFLAGS)-$(NPE)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f flush.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation single CPU f77 Linux
# with visit support
#----------------------------------------------------------------------------
visit:
	cp mpif.null.h mpif.rainer.h
	$(MAKE) $(RESULT) "FC = g77" "FFLAGS = -O4 $(VISITFLAGS)" "LFLAGS2 = $(VISITLD)"\
		"SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation single CPU pgf77 Linux -O4 -tp p6 -pc 64 (PPro, PII, PIII...)
# with visit support
#----------------------------------------------------------------------------
pgfvisit:
	cp mpif.null.h mpif.rainer.h
	$(MAKE) $(RESULT) "FC = pgf77" "FFLAGS = -O4 -tp p6 -pc 64 $(VISITFLAGS)"  "LFLAGS2 = $(VISITLD)"\
		"SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation parallel CPU f77 Linux -O4 (Beowulf et al, experimental)
# with visit support
#----------------------------------------------------------------------------
mpichvisit:
	cp mpif.mpich.h mpif.rainer.h
	$(MAKE) $(RESULT) "FC = mpif77" "FFLAGS = -O4 -D PUREMPI -D PARALLEL $(VISITFLAGS)" "LFLAGS2 = $(VISITLD)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
#----------------------------------------------------------------------------


INC = params.h common6.h

SOURCE = nbody6.f file_init.f ellan.f eigenvalue.f indexx.f \
adjust.f assess.f bindat.f binev.f binout.f binpop.f block.f bodies.f \
brake.f brake2.f brake3.f bsetid.f chaos0.f chaos.f \
check.f checkl.f chrect.f clint.f cloud.f cloud0.f \
cmbody.f cmcorr.f cmfirr.f cmfreg.f coal.f comenv.f core.f corerd.f \
cputim.f data.f decide.f define.f deform.f degen.f delay.f \
dgcore.f dtchck.f eccmod.f ecirc.f edot.f efac2.f efac3.f \
expel.f escape.f events.f evolve.f expand.f fclose.f \
fcloud.f fcorr.f fdisk.f fhalo.f ficorr.f findj.f findm.f \
flyby.f fnuc.f fpcorr.f fpert.f fpoly1.f fpoly2.f freeze.f \
gcinit.f gcint.f giant.f giant3.f gntage.f grrad.f hcorr.f \
hiarch.f hicirc.f hidat.f higrow.f himax.f himax2.f himod.f \
hipop.f hirect.f histab.f hivel.f hmdot.f hmdot2.f hotsys.f \
hrdiag.f hrplot.f hut.f hut2.f iblock.f imf.f imfbd.f imf2.f \
impact.f induce.f input.f insert.f instar.f intgrt.f \
jacobi.f kick.f kick2.f ksapo.f kscorr.f \
ksin2.f ksinit.f ksint.f kslist.f ksmod.f ksperi.f kspert.f \
kspoly.f kspred.f ksrect.f ksreg.f ksres.f ksres2.f ksterm.f \
kstide.f lagr.f levels.f magbrk.f matrix.f mdot.f merge.f \
merge2.f mix.f mloss.f mlwind.f modify.f mrenv.f mtrace.f mydump.f \
nbint.f nblist.f nbpot.f nbrem.f nbrest.f nbsort.f nbtide.f \
newtev.f nstab.f offset.f orbit.f output.f peri.f permit.f \
pfac.f poti.f proto_star.f qtides.f ran2.f reflct.f regint.f \
remove.f rename.f reset.f reset2.f resolv.f rkint.f rl.f roche.f \
rpmax.f rpmax2.f rpmin.f scale.f search.f setup.f setup2.f short.f shrink.f \
sort1.f spiral.f stability.f star.f start.f stepk.f steps.f stumpf.f \
subint.f swap.f sweep.f synch.f tcirc.f tides.f tides2.f \
tides3.f touch.f tpert.f trdot.f trdot2.f trflow.f tstab.f tstep.f \
units.f unpert.f update.f verify.f xtrnl0.f xtrnld.f xtrnlf.f xtrnlp.f \
xtrnlv.f xvpred.f zare.f zcnsts.f zero.f zfuncs.f \
triple.f derqp3.f difsy3.f erel3.f extend.f qpmod3.f stabl3.f \
stablz.f start3.f subsys.f tperi.f trans3.f \
quad.f derqp4.f difsy4.f endreg.f erel4.f ichain.f newreg.f newsys.f \
qpmod4.f rchain.f rsort.f stabl4.f start4.f status.f trans4.f \
cfuncs.f chain.f chstab.f const.f cstab2.f cstab3.f cstab4.f cstab5.f \
derqp.f difsy1.f erel.f hpsort.f inclin.f invert.f ksphys.f physks.f \
qforce.f qpmod.f r2sort.f recoil.f redraw.f select.f slow.f stablc.f \
swcond.f switch.f transk.f transq.f transx.f vector.f xtf.f xtrnlu.f \
ycopy.f ysave.f \
absorb.f chaos2.f chdata.f chfind.f chfirr.f chinit.f chlist.f chmod.f \
chpert.f chpot.f chterm.f expel2.f fchain.f ghost.f giant2.f kcpert.f \
reduce.f reinit.f renew.f setsys.f tchain.f xcpred.f xtpert.f

CUDASOURCE = gpunb.gpu.cu
CUDAOBJECTS = gpunb.gpu.o
#CUDASOURCE = cu_nbody.cu
#CUDAOBJECTS = cu_nbody.o


.SUFFIXES : .o .F

OBJECTS = $(SOURCE:.f=.o)

nbody6:	$(OBJECTS)
	$(LOAD) $(FFLAGS) $(LFLAGS) $(OBJECTS) $(LFLAGS2)

nbody6.gpu: $(OBJECTS) $(CUDAOBJECTS)
	$(LOAD) $(FFLAGS) $(LFLAGS) $(OBJECTS) $(CUDAOBJECTS) $(LFLAGS2)

$(OBJECTS): $(INC)

start.o: energy.f energy_mpi.f fpoly1.f fpoly2.f fpoly1_mpi.f fpoly2_mpi.f

cu_nbody.o: cu_nbody.cu
	nvcc -c -g -I$(SDK_PATH)/common/inc cu_nbody.cu

gpunb.gpu.o: gpunb.gpu.cu
	nvcc -c -g -I$(SDK_PATH)/common/inc gpunb.gpu.cu

clean: 
	\rm *.o adjust.f  binpop.f  cputim.f  data.f  file_init.f  hipop.f  input.f  intgrt.f  modify.f  mydump.f  nbody6.f  output.f  scale.f  start.f  xtrnl0.f

