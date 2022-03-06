#
# Makefile for selecting parallel/ensemble average/single CPU version
# of NBODY6++  June 2003 R.Sp.
#
# modifications: M.A.Hermanns
#
# Defaults for some applications
TAR             = tar
GZIP            = gzip

include Makefile.defs

RESULT          = nbody6
PACKAGE_NAME    = nbody6_package

#VISITFLAGS = -DVISIT
#VISITLD = -L ./lvisit -llvisit_nbody -L/work/Tuc/spurzem/visit/visit/lvisit/lib -llvisit -L/work/Tuc/spurzem/visit/visit/visit20/lib -lvisit

HEADER_FILES    = params.h common6.h

SOURCE_MPI      = energy_mpi.F fpoly1_mpi.F fpoly2_mpi.F flush.F

SOURCE_BASE     = nbody6.F file_init.F ellan.F eigenvalue.F indexx.F \
    adjust.F bindat.F binout.F binpop.F block.F bodies.F  brake.F \
    check.F checkl.F clint.F cloud.F cloud0.F  cmbody.F cmcorr.F \
    cmfirr.F cmfreg.F core.F cputim.F  data.F define.F delay.F efac2.F \
    efac3.F escape.F  events.F evolve.F expand.F fclose.F fcloud.F \
    fcorr.F  ficorr.F findj.F flyby.F fpcorr.F fpert.F fpoly1.F \
    fpoly2.F freeze.F gntage.F hcorr.F hiarch.F hidat.F  himax.F \
    himax2.F hipop.F hivel.F histab.F hmdot.F  hmdot2.F hotsys.F \
    hrdiag.F hrplot.F iblock.F imf.F  imf2.F imfbd.F impact.F input.F \
    insert.F instar.F  intide.F intgrt.F kepler.F kick.F ksapo.F \
    kscorr.F  ksin2.F ksinit.F ksint.F kslist.F ksmod.F ksperi.F \
    kspert.F kspoly.F kspred.F ksrect.F ksreg.F ksres.F  ksres2.F \
    ksterm.F kstide.F lagr.F levels.F matrix.F  mdot.F merge.F merge2.F \
    mloss.F mlwind.F modify.F  mydump.F nbint.F nblist.F nbpot.F \
    nbrem.F nbrest.F  nbtide.F nlmod.F offset.F orbit.F output.F peri.F \
     permit.F pfac.F quad.F ran2.F reflct.F regint.F  remove.F \
    rename.F reset.F reset2.F resolv.F scale.F  search.F setup.F \
    short.F shrink.F subint.F sort1.F  sort3.F star.F start.F stepi.F \
    stepk.F steps.F stumpf.F  tcirc.F tides.F tpert.F triple.F tstep.F \
    units.F  unpert.F update.F verify.F xtrnl0.F xtrnld.F xtrnlf.F \
    xtrnlp.F xtrnlv.F xvpred.F  zcnsts.F zdata.F zfuncs.F  zero.F \
    derqp3.F difsy3.F erel3.F extend.F qpmod3.F  stabl3.F stablz.F \
    start3.F  subsys.F tperi.F trans3.F  derqp4.F difsy4.F endreg.F \
    erel4.F ichain.F newreg.F  newsys.F qpmod4.F rchain.F rsort.F \
    stabl4.F start4.F  status.F trans4.F absorb.F cfuncs.F chfirr.F \
    chfind.F  chinit.F chlist.F chmod.F  chpot.F chterm.F fchain.F \
    ghost.F kcpert.F reduce.F  reinit.F renew.F setsys.F  tchain.F \
    xcpred.F xtpert.F  chain.F chdata.F chstab.F  const.F cstab2.F \
    cstab3.F cstab4.F cstab5.F derqp.F  difsy1.F erel.F hpsort.F \
    inclin.F invert.F ksphys.F  physks.F qforce.F qpmod.F  recoil.F \
    redraw.F r2sort.F  select.F slow.F stablc.F swcond.F  switch.F \
    transk.F  transq.F transx.F vector.F xtf.F xtrnlu.F  ycopy.F \
    ysave.F zare.F nbsort.F nbcleanup.F 

OBJECTS = $(SOURCE:.F=.o)
GEN_SRC = $(SOURCE:.F=.f)
	
.SUFFIXES :
.SUFFIXES : .o .F

# Targets

nbody6:	$(OBJECTS)
	$(LINKER) $(FFLAGS) $(LFLAGS_PRE) $(OBJECTS) $(LFLAGS_POST)

test:
	if [ ! -z "$(POSTPROC)" ]; then eval $(POSTPROC); fi

file_init.F: file_init.F.tmpl
	sed 's/APPEND/SEQUENTIAL/g' file_init.F.tmpl > file_init.F 

.F.o:
	sed -e "s/\(^[\*C].*\)//" $< | \
	$(CPP) $(CPPFLAGS) | ./padding.pl > $*.f
	$(FC) $(FFLAGS) -c $*.f -o $@
	rm -f $*.f

.c.o:
	$(CC) $(CFLAGS) 
	
$(OBJECTS): $(HEADER_FILES)

start.o: energy.F energy_mpi.F fpoly1.F fpoly2.F fpoly1_mpi.F fpoly2_mpi.F

clean: 
	rm -f *.o $(RESULT) $(GEN_SRC)

clean-F:
	rm -f adjust.F  cloud0.F  escape.F hotsys.F  intide.F  modify.F  output.F bindat.F  cputim.F  file_init.F  input.F   ksint.F   mydump.F  scale.F binpop.F  data.F    hipop.F      intgrt.F  lagr.F    nbody6.F  start.F
	
cleanpar: 
	rm -f adjust.o binpop.o cloud0.F cputim.o data.o hipop.o hotsys.o input.o intide.o lagr.o mydump.o nbody6.o output.o \
	intgrt.o start.o scale.o fpoly1_mpi.o fpoly2_mpi.o fpoly1.o fpoly2.o bindat.o \
	energy_mpi.o energy.o ksint.o escape.o

dist:
	$(TAR) cf ../$(PACKAGE_NAME).tar ../run ../src/*.F ../src/*.h ../src/*.tmpl ../src/padding.pl ../src/Makefile.defs.* ../src/Makefile
	$(GZIP) ../$(PACKAGE_NAME).tar
