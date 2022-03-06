      SUBROUTINE DEFINE
*
*
*       Definition of input parameters, options & counters.
*       ---------------------------------------------------
*
*
*       Input parameters
*       ****************
*
*       ---------------------------------------------------------------------
***
* NBODY6:
*
*       KSTART  Control index (1: new run; >1: restart; 3, 4, 5: new params).
*       TCOMP   Maximum computing time in minutes (saved in CPU).
***
* INPUT:
*
*       N       Total particle number (singles + binary c.m.; < NMAX - 2).
*       NFIX    Output frequency of data save or binaries (options 3 & 6).
*       NCRIT   Final particle number (alternative termination criterion).
*       NRAND   Random number sequence skip.
*       NNBMAX  Maximum number of neighbours (< LMAX - 2).
*       NRUN    Run identification index.
*
*       ETAI    Time-step parameter for irregular force polynomial.
*       ETAR    Time-step parameter for regular force polynomial.
*       RS0     Initial radius of neighbour sphere.
*       DTADJ   Time interval for parameter adjustment (N-body units).
*       DELTAT  Output time interval (N-body units).
*       TCRIT   Termination time (N-body units).
*       QE      Energy tolerance (restart if DE/E > 5*QE & KZ(2) > 1).
*       RBAR    Virial cluster radius in pc (set = 1 for isolated cluster).
*       ZMBAR   Mean mass in solar units (1.0 assumed if 0).
*
*       KZ(J)   Non-zero options for alternative paths (see table).
*
*       DTMIN   Time-step criterion for regularization search.
*       RMIN    Distance criterion for regularization search.
*       ETAU    Regularized time-step parameter (6.28/ETAU steps/orbit).
*       ECLOSE  Binding energy per unit mass for hard binary (positive).
*       GMIN    Relative two-body perturbation for unperturbed motion.
*       GMAX    Secondary termination parameter for soft KS binaries.
***
* INPUT: if (kz(4).gt.0)
*
*       DELTAS  Output interval for binary search (in TCR; #4, suppressed).
*       ORBITS  Minimum periods for binary output (level 1).
*       GPRINT  Perturbation thresholds for binary output (9 levels).
***
* DATA:
*
*       ALPHAS  Power-law index for initial mass function (used if #20 < 2).
*       BODY1   Maximum particle mass before scaling (KZ(20): solar mass).
*       BODYN   Minimum particle mass before scaling.
*       NBIN0   Number of primordial binaries (for IMF2 with KZ(20) > 1).
*       ZMET    Metal abundance (in range 0.03 - 0.0001).
*       EPOCH0  Evolutionary epoch (in 10**6 yrs).
***
* SCALE:
*
*       Q       Virial ratio (routine SCALE; Q = 0.5 for equilibrium).
*       VXROT   XY-velocity scaling factor (> 0 for solid-body rotation).
*       VZROT   Z-velocity scaling factor (not used if VXROT = 0).
*       RTIDE   Unscaled tidal radius (#14 >= 2; otherwise copied to RSPH2).
***
* XTRNL0: if (kz(14).eq.2)
*
*       GMG     Point-mass galaxy (solar masses, linearized circular orbit).
*       RG0     Central distance (in kpc).
*         if (kz(14).eq.3)
*       GMG     Point-mass galaxy (solar masses).
*       DISK    Mass of Miyamoto disk (solar masses).
*       A       Softening length in Miyamoto potential (in kpc).
*       B       Vertical softening length (kpc).
*       VCIRC   Galactic circular velocity (km/sec) at RCIRC (=0: no halo).
*       RCIRC   Central distance for VCIRC with logarithmic potential (kpc).
*       RG      Initial position; GMG+DISK=0, VG(3)=0: A(1+E)=RG(1), E=RG(2).
*       VG      Initial cluster velocity vector (km/sec).
***
* HOTSYS: if (kz(29).gt.0)
*
*       SIGMA0  Hot initial velocities in km/sec (routine HOTSYS; option 29).
***
* BINPOP: if (kz(8).eq.1.or.kz(8).gt.2)
*
*       NBIN    Number of primordial binaries (routine BINPOP; option 8).
*       SEMI    Semi-major axis in model units (all equal if RANGE = 0).
*       ECC     Initial eccentricity (< 0 for thermal distribution).
*       RATIO   Mass ratio M1/(M1 + M2); (= 1.0: M1 = M2 = <M>; not #20 > 1).
*       RANGE   Range in SEMI for uniform logarithmic distribution (> 0).
*       NSKIP   Binary frequency of mass spectrum (#20 < 2; body #1 first).
*       IDORM   Indicator for dormant binaries (>0: merged components).
*       ICIRC   Eigenevolution & period distribution (RANGE: minimum period).
***
* HIPOP: if (kz(8).gt.0.and.kz(11).gt.1)
*
*       NHI     Number of primordial hierarchies (routine HIPOP; #11 > 1).
*       SEMI    Semi-major axis in model units (all equal if RANGE = 0).
*       ECC     Initial eccentricity (< 0 for thermal distribution).
*       RATIO   Mass ratio (= 1.0: M1 = M2; random in [0.5-0.9]).
*       RANGE   Range in SEMI for uniform logarithmic distribution (> 0).
*       ICIRC   Circularization & collision check (not implemented yet).
***
* INTIDE:
*
*       RSTAR   Size of typical star in A.U. (routine INTIDE; option 27).
*       IMS     # idealized main-sequence stars.
*       IEV     # idealized evolved stars.
*       RMS     Scale factor for main-sequence radii (>0: fudge factor).
*       REV     Scale factor for evolved radii (initial size RSTAR).
***
* CLOUD0: if (kz(13).gt.0)
*
*       NCL     Number of interstellar clouds (routine CLOUD0; option 13).
*       RB2     Radius of cloud boundary in pc (square is saved).
*       VCL     Mean cloud velocity in km/sec.
*       SIGMA   Velocity dispersion (#13 > 1: Gaussian).
*       CLM     Individual cloud masses in solar masses (maximum MCL).
*       RCL2    Half-mass radii of clouds in pc (square is saved).
*       ---------------------------------------------------------------------
*
*
*       Options KZ(J)
*       *************
*
*       ---------------------------------------------------------------------
*       1  COMMON save on unit 1 at end of run (=2: every 100*NMAX steps).
*       2  COMMON save on unit 2 at output (=1); restart if DE/E > 5*QE (=2).
*       3  Basic data on unit 3 at output time (frequency NFIX).
*       4  Binary diagnostics on unit 4 (# threshold levels = KZ(4) < 10).
*       5  Initial conditions (#22 =0; =0: uniform & isotropic; =1: Plummer).
*       6  Output of significant & regularized binaries (=1, 2, 3 & 4).
*       7  Lagrangian radii (>0: RSCALE; =2, 3, 4: output on unit 6 & 7).
*       8  Primordial binaries (=1 & >=3; >0: BINOUT; >2: BINDAT; >3: HIDAT).
*       9  Individual bodies printed at output time (MIN(5**KZ9,NTOT)).
*      10  Diagnostic KS output (>0: begin; >1: end; >=3: each step).
*      11  Hierarchical systems (=1: diagnostics; =2: primordial; =3: both).
*      12  HR diagnostics of evolving stars (interval DTPLOT).
*      13  Interstellar clouds (=1: constant velocity; >1: Gaussian).
*      14  External force (=1: standard tidal field; =2: point-mass galaxy;
*             =3: point-mass + disk + logarithmic halo (any combination).
*      15  Triple, quad, chain (#30 > 0) or merger search (>1: full output).
*      16  Updating of regularization parameters (RMIN, DTMIN & ECLOSE).
*      17  Modification of ETAI, ETAR (>=1) and ETAU (>1) by tolerance QE.
*      18  Neighbour additions in CHECKL (>0: high-velocity; >1: all types).
*      19  Mass loss (=1: old supernova scheme; =3: Eggleton, Tout & Hurley).
*      20  Initial mass function (=1: Scalo; =2,4,6: Kroupa; =3,5: Eggleton).
*      21  Extra output line (MODEL #, TCOMP, DMIN, AMIN, RMAX & RSMIN).
*      22  Initial m,r,v on #10 (=1: output; >=2: input; >2: no scaling).
*      23  Escaper removal (>1: diagnostics in file ESC; >2: angles unit #6).
*      24  Initial conditions for subsystem (routine SCALE; KZ(24) = #).
*      25  Partial reflection of KS binary orbit (GAMMA < GMIN; suppressed).
*      26  Slow-down of two-body motion (>=1: KS binary; =2: chain binary).
*      27  Tidal circularization & collisions (R_coll = 0.75*(R_1 + R_2)).
*      28  Magnetic braking and gravitational radiation (options 19 & 27).
*      29  Boundary reflection for hot system (suppressed).
*      30  Chain regularization (>=2: main output; >2: diagnostic output).
*      31  Centre of mass correction after energy check.
*      32  Increase of output intervals (based on single particle energy).
*      33  Block-step diagnostics at main output (>=1: STEP; =2: STEPR).
*      34  Roche lobe overflow (not implemented yet).
*      35  Time offset (global time from TTOT = TIME + TOFF).
*      36  Step reduction for hierarchical systems.
*      37  Fast time-step criterion (>0: STEP; >1: STEPR).
*      38  No force polynomial corrections (I <= N; block-step version).
*      39  No unique density centre (skips velocity modification of RS(I)).
*      40  Increase of neighbour numbers if <NNB> < NNBMAX/2.
*       ---------------------------------------------------------------------
*
*
*       Output counters
*       ***************
*
*       ---------------------------------------------------------------------
*       NSTEPI  Irregular integration steps.
*       NSTEPR  Regular integration steps.
*       NSTEPU  Regularized integration steps.
*       NNPRED  Coordinate & velocity predictions of all particles.
*       NBPRED  Coordinate & velocity prediction of neighbours (NNB counted).
*       NBCORR  Force polynomial corrections.
*       NBFULL  Too many neighbours with standard criterion.
*       NBVOID  No neighbours inside 1.26 times the basic sphere radius.
*       NICONV  Irregular step reduction (force convergence test).
*       NRCONV  Regular step reduction (Hermite) or increase (block-steps).
*       NBSMIN  Retained neighbours inside 2*RS (STEP < SMIN).
*       NLSMIN  Small step neighbours selected from other neighbour lists.
*       NBDIS   Second component of recent KS pair added as neighbour (#18).
*       NBDIS2  Second component of old KS pair added as neighbour (#18 > 1).
*       NCMDER  C.m. values for force derivatives of KS component.
*       NBDER   Large F3DOT corrections not included in D3 & D3R.
*       NFAST   Fast particles included in LISTV (option 18).
*       NBFAST  Fast particles included in neighbour list (option 18).
*       NBLOCK  Number of blocks (block-step version).
*       NMDOT   Mass loss events (option 19).
*       NBSTAT  Diagnostic data on binary interactions (option 4).
*       NKSTRY  Two-body regularization attempts.
*       NKSREG  Total KS regularizations.
*       NKSHYP  Hyperbolic KS regularizations.
*       NKSPER  Unperturbed KS binary orbits.
*       NPRECT  Initialization of NKSPER after exceeding 2*10**9.
*       NKSREF  Partial reflections of KS binary (option 25; suppressed).
*       NKSMOD  Slow KS motion restarts (option 26).
*       NTTRY   Search for triple, quad & chain regularization or mergers.
*       NTRIP   Three-body regularizations (option 15).
*       NQUAD   Four-body regularizations (option 15).
*       NCHAIN  Chain regularizations (options 15 & 30).
*       NMERG   Mergers of stable triples or quadruples (option 15).
*       NEWHI   New hierarchical systems (counted by routine HIARCH).
*       NSTEPT  Triple regularization integration steps (option 15).
*       NSTEPQ  Quadruple regularization integration steps (option 15).
*       NSTEPC  Chain regularization steps (# DIFSY calls).
*       NDISS   Tidal dissipations at pericentre (option 27).
*       NTIDE   Tidal captures from hyperbolic motion (option 27).
*       NSYNC   Number of synchronous binaries (option 27).
*       NCOLL   Stellar collisions (option 27).
*       NSESC   Escaped single particles (option 23).
*       NBESC   Escaped binaries (option 23).
*       NMESC   Escaped mergers (options 15 & 23).
*       NRG     Red giants.
*       NHE     Helium stars.
*       NRS     Red supergiants.
*       NNH     Naked Helium stars.
*       NWD     White dwarfs.
*       NSN     Neutron stars.
*       NBH     Black holes.
*       NBS     Blue stragglers.
*       ---------------------------------------------------------------------
*
*
*       Stellar evolution types
*       ***********************
*
*       ---------------------------------------------------------------------
*       0       Low main sequence (M < 0.7).
*       1       Main sequence.
*       2       Hertzsprung gap (HG).
*       3       Red giant.
*       4       Core Helium burning.
*       5       First AGB.
*       6       Second AGB.
*       7       Helium main sequence.
*       8       Helium HG.
*       9       Helium GB.
*      10       Helium white dwarf.
*      11       Carbon-Oxygen white dwarf.
*      12       Oxygen-Neon white dwarf.
*      13       Neutron star.
*      14       Black hole.
*      15       Massless supernova remnant.
*      19       Circularizing binary (c.m. value).
*      20       Circularized binary.
*      21       First Roche stage (inactive).
*      22       Second Roche stage.
*       ---------------------------------------------------------------------
*
      RETURN
*
      END
