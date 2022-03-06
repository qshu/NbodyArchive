FFLAGS = -O3
FORTRAN = f77
#Note: 'make printh' requires a special function
#for making nice listings of individual routines.

SOURCE = \
absorb.f chfirr.f chfind.f chinit.f chlist.f chmod.f \
chpot.f chterm.f fchain.f ghost.f kcpert.f reduce.f \
reinit.f renew.f setsys.f tchain.f xcpred.f xtpert.f

print:
	@- \rm -f NCHAIN.TEXT
	@cat $(SOURCE) > NCHAIN.TEXT
