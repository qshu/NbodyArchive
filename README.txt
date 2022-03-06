The code nbody6.tar.gz is the standard version of NBODY6.

New way to calculate forces for c.m. bodies 2 Feb 2004.

A few rare bugs fixed 1/2 April and 26 April 2004.
Major improvements for the 21 May version and a few more 30 June.

PS. 20 June: It turns out there is a compiler problem (also with g77)
in select.f. This can be bypassed by including SAVE.

The manual for NBODY6 can be found in the postscript file man6.ps.
The new stability criterion referred to is not yet available.

22 September: new version includes integration of tidal tail. It is
completely compatible with the previous version. The manual (man6.ps)
has also been updated to describe the new features in detail.

22 October: external Plummer potential added with #14 = 4 and manual
updated. Also a few minor improvements. The compiler bug with f77 in
routine transq.f showed up again and was fixed by an artificial stop.

21 December: several improvements; Plummer potential included for 3D.

19 February 2005: a few bug fixes; one in regint.f is important.

1 April 2005: a large number of bug fixes and some improvements.

27 October 2005: important bug fixes for chain regularization.

30 November 2005: improved scaling for cluster inside Plummer model.

6 February 2006: expanded to extra input models as in public NBODY4.

8 September 2006: several improvements, extra input, larger common.

13 November 2006: bug fix re virial theorem for 3D potential.
 
1 Dec. 2007: new expanded NBODY6 in nbody6.tar.gz (old old6.tar.gz).

13 May 2008: bug fix for initial time-steps (routine steps.f).

23 June 2008: slightly improved standard version also for GPU or SSE.

1 Sept. 2008: a few bug fixes for standard version.

19 Sept. 2008: implementation of NS and BH binaries (option #28).

7 Oct. 2008: routine SWEEP for regularization of wider binaries (#8).

24 Oct. 2008: improvements and bug fixes for circularized binaries.

28 Nov. 2008: several changes concerning stellar evolution procedures.
 
31 Dec. 2008: two improvements related to common envelope evolution.

13 Jan. 2009: two important bug fixes in chain regularization.

8 May 2009: improved random number generator (Press et al., 2nd ed).

10 May 2009: re-design of chain coalescence with mass loss.

2 Sep. 2009: several bug fixes included, especially in chain.f.
