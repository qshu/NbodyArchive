.KEEP_STATE:
FFLAGS = -O3
#FFLAGS = -pg
FC = g77
#FC = f77

VPATH = ../Chain:../Nchain

SOURCE = \
nbody6.f adjust.f bindat.f binout.f binpop.f block.f bodies.f \
brake.f check.f checkl.f clint.f cloud.f cloud0.f cmbody.f \
cmcorr.f cmfirr.f cmfreg.f core.f cputim.f data.f define.f \
delay.f efac2.f efac3.f energy.f escape.f events.f evolve.f \
expand.f fclose.f fcloud.f fcorr.f ficorr.f findj.f flyby.f \
fpcorr.f fpert.f fpoly1.f fpoly2.f freeze.f hcorr.f hiarch.f \
hidat.f himax.f himax2.f hipop.f histab.f hivel.f hmdot.f hmdot2.f \
hotsys.f hrdiag.f hrplot.f iblock.f imf.f imfbd.f imf2.f impact.f \
inext.f input.f insert.f instar.f intgrt.f intide.f kepler.f kick.f \
ksapo.f kscorr.f ksin2.f ksinit.f ksint.f kslist.f ksmod.f \
ksperi.f kspert.f kspoly.f kspred.f ksrect.f ksreg.f ksres.f \
ksres2.f ksterm.f kstide.f lagr.f levels.f matrix.f merge.f \
merge2.f mdot.f mloss.f mlwind.f modify.f mydump.f nbint.f nblist.f \
nbpot.f nbrem.f nbrest.f nbtide.f nbsort.f nlmod.f offset.f output.f \
orbit.f peri.f permit.f pfac.f ran2.f reflct.f regint.f remove.f \
rename.f reset.f reset2.f resolv.f scale.f search.f setup.f \
shrink.f sort1.f sort2.f star.f start.f stepi.f stepk.f steps.f \
stumpf.f subint.f tcirc.f tides.f tstep.f tpert.f units.f unpert.f \
update.f verify.f xtrnl0.f xtrnld.f xtrnlf.f xtrnlp.f xtrnlv.f \
xvpred.f zare.f zcnsts.f zdata.f zfuncs.f zero.f \
triple.f derqp3.f difsy3.f erel3.f extend.f qpmod3.f stabl3.f \
stablz.f start3.f subsys.f tperi.f trans3.f \
quad.f derqp4.f difsy4.f endreg.f erel4.f ichain.f newreg.f newsys.f \
qpmod4.f rchain.f rsort.f stabl4.f start4.f status.f trans4.f \
cfuncs.f chain.f chstab.f const.f cstab2.f cstab3.f cstab4.f cstab5.f \
derqp.f difsy1.f erel.f hpsort.f inclin.f invert.f ksphys.f physks.f \
qforce.f qpmod.f r2sort.f recoil.f redraw.f select.f slow.f stablc.f \
swcond.f switch.f transk.f transq.f transx.f vector.f xtf.f xtrnlu.f \
ycopy.f ysave.f \
absorb.f chdata.f chfind.f chfirr.f chinit.f chlist.f chmod.f chpot.f \
chterm.f fchain.f ghost.f kcpert.f reduce.f reinit.f renew.f setsys.f \
tchain.f xcpred.f xtpert.f

OBJECTS = $(SOURCE:.f=.o)

nbody6:	$(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o nbody6

print:
	@- \rm -f NBODY6.TEXT
	@cat $(SOURCE) > NBODY6.TEXT
