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
# make crayt3e NPE=n              (full parallel version using SHMEM CRAY T3E)
# make crayt3e_noshmem NPE=n      (full parallel version no SHMEM CRAY T3E)
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
VISITFLAGS = -DVISIT
VISITLD = -L ./lvisit -llvisit_nbody -L/work/Tuc/spurzem/visit/visit/lvisit/lib -llvisit -L/work/Tuc/spurzem/visit/visit/visit20/lib -lvisit

#----------------------------------------------------------------------------
# Workstation single CPU g77 Linux -O6
#----------------------------------------------------------------------------
standard: 
	cp mpif.null.h mpif.h
	$(MAKE) $(RESULT) "FC = g77" "FFLAGS = -O6" "SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation single CPU pgf77 Linux generic
#----------------------------------------------------------------------------
pgf:
	cp mpif.null.h mpif.h
	$(MAKE) $(RESULT) "FC = pgf77" "FFLAGS = -O4" \
		"SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation single CPU pgf77 Linux (PIII, PII, PPro)
#----------------------------------------------------------------------------
pgfp3:
	cp mpif.null.h mpif.h
	$(MAKE) $(RESULT) "FC = pgf77" "FFLAGS = -O4 -tp p6 -pc 64 -fast" \
		"SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation single CPU pgf77 Linux (P4)
#----------------------------------------------------------------------------
pgfp4:
	cp mpif.null.h mpif.h
	$(MAKE) $(RESULT) "FC = pgf77" "FFLAGS = -O4 -tp p7 -pc 64 -fast -fastsse" \
		"SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation parallel CPU f77 Linux (Beowulf generic mpif77)
#----------------------------------------------------------------------------
mpich:
	cp mpif.mpich.h mpif.h
	$(MAKE) $(RESULT) "FC = mpif77" "FFLAGS = -O4 -D PUREMPI -D PARALLEL" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
mpich_ensemble:
	cp mpif.mpich.h mpif.h
	$(MAKE) $(RESULT) "FC = mpif77" "FFLAGS = -O4 -D ENSEMBLE" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation parallel CPU pgf77 Linux (Beowulf pgf77 P4 optimised)
#----------------------------------------------------------------------------
mpichpgf:
	cp mpif.mpich.h mpif.h
	$(MAKE) $(RESULT) "LOAD = mpif77-pgi" "FC = pgf77" "FFLAGS = -O4 -tp p7 -pc 64 -fast -fastsse -D PUREMPI -D PARALLEL" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
mpichpgf_ensemble:
	cp mpif.mpich.h mpif.h
	$(MAKE) $(RESULT) "LOAD = mpif77-pgi" "FC = pgf77" "FFLAGS = -O4 -tp p7 -pc 64 -fast -fastsse -D ENSEMBLE" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation single CPU f77 Solaris 
#----------------------------------------------------------------------------
sun:
	cp mpif.null.h mpif.h
	$(MAKE) $(RESULT) "FC = f77" "FFLAGS = -O4 " "SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Parallel Sun
#----------------------------------------------------------------------------
sun_e10k:
	cp mpif.sun.h mpif.h
	$(MAKE) $(RESULT) "FFLAGS = -D PUREMPI -D PARALLEL" \
		"FC = tmf77" "LFLAGS = -lmpi $(LFLAGS)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
sun_e10k_ensemble:
	cp mpif.sun.h mpif.h
	$(MAKE) $(RESULT) "FFLAGS = -D ENSEMBLE" \
		"FC = tmf77" "LFLAGS = -lmpi $(LFLAGS)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
sun_e10k_fast:
	cp mpif.sun.h mpif.h
	$(MAKE) $(RESULT) "FFLAGS = -fast -xarch=v9a -D PUREMPI -D PARALLEL" \
		"FC = tmf77" "LFLAGS = -lmpi $(LFLAGS)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
sun_e10k_fast_ensemble:
	cp mpif.sun.h mpif.h
	$(MAKE) $(RESULT) "FFLAGS = -fast -xarch=v9a -D ENSEMBLE" \
		"FC = tmf77" "LFLAGS = -lmpi $(LFLAGS)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
#----------------------------------------------------------------------------
# CRAY T3E Fully parallel Mode including SHMEM
#----------------------------------------------------------------------------
crayt3e:
	cp mpif.t3e.h mpif.h
	cp flush.t3e.f flush.f
	$(MAKE) $(RESULT) "FFLAGS = -d p -D SHMEM -D PARALLEL" \
		"FC = f90" "LFLAGS = -X$(NPE) $(LFLAGS)-$(NPE)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f flush.f $(SOURCE)"
crayt3e_ensemble:
	cp mpif.t3e.h mpif.h
	cp flush.t3e.f flush.f
	$(MAKE) $(RESULT) "FFLAGS = -d p -D ENSEMBLE" \
		"FC = f90" "LFLAGS = -X$(NPE) $(LFLAGS)-$(NPE)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f flush.f $(SOURCE)"
crayt3e_noshmem:
	cp mpif.t3e.h mpif.h
	cp flush.t3e.f flush.f
	$(MAKE) $(RESULT) "FFLAGS = -d p -D PUREMPI -D PARALLEL" \
		"FC = f90" "LFLAGS = -X$(NPE) $(LFLAGS)-$(NPE)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f flush.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation single CPU f77 Linux
# with visit support
#----------------------------------------------------------------------------
visit:
	cp mpif.null.h mpif.h
	$(MAKE) $(RESULT) "FC = g77" "FFLAGS = -O4 $(VISITFLAGS)" "LFLAGS2 = $(VISITLD)"\
		"SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation single CPU pgf77 Linux -O4 -tp p6 -pc 64 (PPro, PII, PIII...)
# with visit support
#----------------------------------------------------------------------------
pgfvisit:
	cp mpif.null.h mpif.h
	$(MAKE) $(RESULT) "FC = pgf77" "FFLAGS = -O4 -tp p6 -pc 64 $(VISITFLAGS)"  "LFLAGS2 = $(VISITLD)"\
		"SOURCE = energy.f $(SOURCE)"
#----------------------------------------------------------------------------
# Workstation parallel CPU f77 Linux -O4 (Beowulf et al, experimental)
# with visit support
#----------------------------------------------------------------------------
mpichvisit:
	cp mpif.mpich.h mpif.h
	$(MAKE) $(RESULT) "FC = mpif77" "FFLAGS = -O4 -D PUREMPI -D PARALLEL $(VISITFLAGS)" "LFLAGS2 = $(VISITLD)" \
		"SOURCE = energy_mpi.f fpoly1_mpi.f fpoly2_mpi.f $(SOURCE)"
#----------------------------------------------------------------------------


INC = params.h common6.h

SOURCE = nbody6.f file_init.f ellan.f eigenvalue.f indexx.f \
adjust.f bindat.f binout.f binpop.f block.f bodies.f \
brake.f check.f checkl.f clint.f cloud.f cloud0.f \
cmbody.f cmcorr.f cmfirr.f cmfreg.f core.f cputim.f data.f define.f \
delay.f efac2.f efac3.f escape.f events.f evolve.f expand.f \
fclose.f fcloud.f fcorr.f ficorr.f findj.f flyby.f fpcorr.f fpert.f \
fpoly1.f fpoly2.f freeze.f hcorr.f hiarch.f hidat.f himax.f \
hipop.f hivel.f histab.f hmdot.f hmdot2.f hotsys.f \
hrdiag.f hrplot.f iblock.f imf.f imf2.f imfbd.f impact.f input.f \
insert.f instar.f intide.f \
intgrt.f kepler.f kick.f ksapo.f \
kscorr.f ksin2.f ksinit.f ksint.f kslist.f ksmod.f ksperi.f kspert.f \
kspoly.f kspred.f ksrect.f ksreg.f ksres.f ksres2.f ksterm.f kstide.f lagr.f \
levels.f matrix.f mdot.f merge.f merge2.f mloss.f mlwind.f \
modify.f mydump.f nbint.f nblist.f nbpot.f \
nbrem.f nbrest.f nbsort.f nbtide.f nlmod.f offset.f orbit.f output.f peri.f \
permit.f pfac.f quad.f ran2.f reflct.f regint.f remove.f rename.f \
reset.f reset2.f resolv.f scale.f search.f setup.f short.f shrink.f \
subint.f sort1.f sort3.f star.f start.f stepi.f stepk.f steps.f stumpf.f \
tcirc.f tides.f tpert.f triple.f tstep.f units.f unpert.f update.f \
verify.f xtrnl0.f xtrnld.f xtrnlf.f xtrnlp.f xtrnlv.f xvpred.f \
zcnsts.f zdata.f zfuncs.f zero.f \
derqp3.f difsy3.f erel3.f extend.f qpmod3.f stabl3.f stablz.f start3.f \
subsys.f tperi.f trans3.f \
derqp4.f difsy4.f endreg.f erel4.f ichain.f newreg.f newsys.f qpmod4.f \
rchain.f rsort.f stabl4.f start4.f status.f trans4.f \
absorb.f cfuncs.f chfirr.f chfind.f chinit.f chlist.f chmod.f \
chpot.f chterm.f fchain.f ghost.f kcpert.f reduce.f \
reinit.f renew.f setsys.f tchain.f xcpred.f xtpert.f \
chain.f chdata.f chstab.f const.f cstab2.f cstab3.f cstab4.f cstab5.f derqp.f \
difsy1.f erel.f hpsort.f inclin.f invert.f ksphys.f physks.f qforce.f qpmod.f \
recoil.f redraw.f r2sort.f select.f slow.f stablc.f swcond.f \
switch.f transk.f transq.f transx.f vector.f xtf.f xtrnlu.f \
ycopy.f ysave.f zare.f


OBJECTS = $(SOURCE:.f=.o)

nbody6:	$(OBJECTS)
	$(LOAD) $(FFLAGS) $(LFLAGS) $(OBJECTS) $(LFLAGS2)

$(OBJECTS): $(INC)

start.o: energy.f energy_mpi.f fpoly1.f fpoly2.f fpoly1_mpi.f fpoly2_mpi.f

clean: 
	\rm *.o adjust.f  cloud0.f  escape.f     hotsys.f  intide.f  modify.f  output.f bindat.f  cputim.f  file_init.f  input.f   ksint.f   mydump.f  scale.f binpop.f  data.f    hipop.f      intgrt.f  lagr.f    nbody6.f  start.f

cleanpar: 
	\rm adjust.o binpop.o cloud0.f cputim.o data.o hipop.o hotsys.o input.o intide.o lagr.o mydump.o nbody6.o output.o \
intgrt.o start.o scale.o fpoly1_mpi.o fpoly2_mpi.o fpoly1.o fpoly2.o bindat.o \
energy_mpi.o energy.o ksint.o escape.o
