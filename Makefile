#
# Makefile for selecting parallel/ensemble average/single CPU version
# of NBODY6++  May 98 R.Sp.
#
# USAGE:
# make nbody6                              (single CPU version)
# make parallel                            (full parallel version using SHMEM)
# make parnoshmem                          (full parallel version no SHMEM)
# make ensemble                            (parallel version for ensemble)
# make clean                               (remove all .o)
# make cleanpar                            (remove all .o originating from .F)
#
# and select in the following the proper FFLAGS/LFLAG 
#
#FFLAGS = -O3 -d p # -g -p -i4 -r8
#FFLAGS = -O3 -d p -D SHMEM -D PARALLEL 
#FFLAGS = -O3 -d p -D PARALLEL
#FFLAGS = -O3 -d p -D ENSEMBLE
# FFLAG Linux PGI
#FFLAGS = -O4 -tp p6 -pc64
# FFLAG Linux:
FFLAGS = -O6
#FFLAGS SUN
#FFLAGS = -O2 -xtarget=ss2
#
# Select where executable goes and how many PE's:
LFLAGS = -o Run/nbody6
#LFLAGS = -X8 -o Run/nbody6
#LFLAGS = -X16 -o Run/nbody6-16
#LFLAGS = -X32 -o Run/nbody6-32
#LFLAGS = -X64 -o Run/nbody6-64
#LFLAGS = -X128 -o Run/nbody6-128
#LFLAGS = -X256 -o Run/nbody6-256
#LFLAGS = -X512 -o Run/nbody6-512

#FC = f90
#FC = pgf77
FC = g77

INC = params.h common6.h

SOURCE1 = nbody6_pk.f ellan.f eigenvalue.f indexx.f \
adjust_pk.f bindat_pk.f binout.f binpop_pk.f block.f bodies.f \
checkl.f clint.f cloud.f cloud0.f \
cmbody.f cmcorr.f core.f cputim.f data_pk.f define_pk.f \
efac2.f efac3.f escape.f events.f evolve.f fclose.f \
fcloud.f fcorr.f ficorr.f fpert.f fpoly1.f fpoly2.f freeze.f hotsys.f \
hrdiag.f iblock.f imf.f imf2_pk.f impact.f input_pk.f instar.f intide.f \
kepler.f kscorr.f ksinit.f ksint_pk.f kslist.f ksmod.f ksperi.f kspert.f \
kspoly.f kspred.f ksreg.f ksres.f ksres2.f ksterm.f kstide.f lagr.f \
levels.f merge.f mdot.f mloss.f modify.f mydump_pk.f nbpot.f \
nbrem.f nbrest.f nbtide.f nlmod.f output_pk.f peri.f ran2_pk.f reflct.f \
remove.f rename.f reset.f resolv.f scale.f search.f signal.f sort1.f \
sort2.f star.f start_pk.f tides.f tpert.f tstep.f stepi.f unpert.f update.f \
verify.f xtrnl0_pk.f xtrnld.f xtrnlf.f xtrnlp.f xtrnlv.f xvpred.f zero.f \
triple.f derqp3.f difsy3.f erel3.f extend.f qpmod3.f stabl3.f stablz.f \
start3.f subsys.f tperi.f trans3.f \
quad.f derqp4.f difsy4.f endreg.f erel4.f ichain.f \
newreg.f newsys.f qpmod4.f rchain.f rsort.f stabl4.f start4.f \
status.f trans4.f \
chfind.f chfirr.f chlist.f chpot.f fchain.f tchain.f kcpert.f xcpred.f \
chinit.f chterm.f absorb.f chmod.f ghost.f reduce.f reinit.f setsys.f \
chain.f const.f derqp.f difsy1.f hpsort.f ksphys.f matrix.f physks.f \
qforce.f recoil.f redraw.f r2sort.f select.f swcond.f switch.f transq.f \
transx.f xtf.f xtpert.f xtrnlu.f ycopy.f ysave.f \
stepk.f inext.f resort.f sort3.f nbsort.f check.f cmfirr.f cmfreg.f \
intgrt.f nbint.f regint.f steps.f fpcorr.f nblist.f

SOURCE2 = energy.f

SOURCE2_PARA = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f

SOURCE = $(SOURCE1) $(SOURCE2)

SOURCE_PARA = $(SOURCE1) $(SOURCE2_PARA)

OBJECTS = $(SOURCE:.f=.o)

OBJECTS_PARA = $(SOURCE_PARA:.f=.o)

nbody6:	$(OBJECTS)
	echo "Compilation for single CPU (workstation) done"
	$(FC) $(FFLAGS) $(LFLAGS) $(OBJECTS)

parallel: $(OBJECTS_PARA)
	echo "Compilation for big parallel run done"
	$(FC) $(FFLAGS) $(LFLAGS) $(OBJECTS_PARA)

ensemble: $(OBJECTS)
	echo "Compilation for ensemble parallel runs done"
	$(FC) $(FFLAGS) $(LFLAGS) $(OBJECTS)

$(OBJECTS): $(INC)

$(OBJECTS_PARA): $(INC)

clean: 
	\rm *.o

cleanpar: 
	\rm adjust.o cputim.o data.o input.o lagr.o nbody6.o output.o \
intgrt.o start.o scale.o fpoly1_mpi.o fpoly2_mpi.o fpoly1.o fpoly2.o \
energy_mpi.o energy.o 
