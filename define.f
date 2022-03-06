      SUBROUTINE DEFINE
*
*
*       Definition of parameters, options & counters.
*       ----------------------------------------------
*
*
*       Input parameters
*       ****************
*
*       --------------------------------------------------------------------
***
* NBODY4: 
*
*       KSTART  Control index (1: new run; >1: restart; 3, 4, 5: new params).
*       TCOMP   Maximum computing time in minutes (saved in CPU).
*       GPID    Number of boards (counting from 0).
***
* INPUT:
*
*       N       Total particle number (singles + binary c.m; < NMAX - 2).
*       NFIX    Output frequency of data save or binaries (#3 & 6).
*       NCRIT   Final particle number (alternative termination criterion).
*       NRAND   Random number sequence skip.
*       NRUN    Run identification index.
*
*       ETA     Time-step parameter for total force polynomial.
*       DTADJ   Time interval for parameter adjustment (N-body units).
*       DELTAT  Output time interval (N-body units).
*       TCRIT   Termination time (N-body units).
*       QE      Energy tolerance (restart if DE/E > 5*QE & KZ(2) > 1).
*       RBAR    Virial cluster radius in pc (set = 0 for isolated cluster).
*       ZMBAR   Mean mass in solar units for equal-mass case (1.0 used if 0).
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
* INPUT: if(kz(4).gt.0)
*
*       DELTAS  Output interval for binary search (in TCR; #4, suppressed).
*       ORBITS  Minimum periods for binary output (level 1).
*       GPRINT  Perturbation thresholds for binary output (9 levels).
***
* DATA:
*
*       ALPHA   Power-law index for initial mass function (used if KZ(20)<2)).
*       BODY1   Maximum particle mass before scaling (KZ(20): solar mass).
*       BODYN   Minimum particle mass before scaling (also after scaling).
*       NBIN0   Number of primordial binaries (can be -ve: see IMF).
*       ZMET    Metal abundance (range 0.03 - 0.0001).
*       EPOCH0  Evolutionary epoch (in Myr).
*       DTPLOT  Plotting interval for HRDIAG (N-body units).
***
* SCALEX:
*
*       Q       Virial ratio (Q = 0.5 for equilibrium).
*       VXROT   XY-velocity scaling factor (> 0 for solid-body rotation).
*       VZROT   Z-velocity scaling factor (not used if VXROT = 0).
*       RTIDE   Unscaled tidal radius (#14 > 1; otherwise copied to RSPH2).
***
* SCALEX: if(kz(24).gt.0)
*
*       X(3),V(3) Position and velocity for initial subsystems (KZ(24) lines).
***
* XTRNL0: if(kz(14).gt.1)
*
*       VG0     Cluster circular velocity (km/sec; routine XTRNL0; #14 >= 2).
*       FCW     Dimensionless parameter in Chernoff-Weinberg models.
*       ZMS     Scaled cluster mass (in solar units).
*       ECC     Eccentricity for 3D orbit (only with #14 > 2).
* but for (kz(14).eq.3) 
*       VG0     Galactic circular speed at radius FCW (km/sec). 
*       FCW     Initial galactocentric radius in kpc. 
*       ZMS     DUMMY
***
* HOTSYS: if(kz(29).gt.0) 
*
*       SIGMA0  Hot initial velocity in km/sec (routine HOTSYS; #29).
***
* BINPOP: if(kz(8).eq.1.or.kz(8).gt.2)
*
*       SEMI    Semi-major axis in model units (all equal if RANGE = 0).
*       ECC     Initial eccentricity (< 0 for thermal distribution).
*       RATIO   Mass ratio M1/(M1 + M2); (= 1.0: M1 = M2 = <M>).
*       RANGE   Range in SEMI for uniform logarithmic distribution (> 0).
*       NSKIP   Binary frequency of mass spectrum (#20 < 2; body #1 first).
*       IDORM   Indicator for dormant binaries (>0: merged components).
*       ICIRC   Indicator for eigen-evolution (RANGE: minimum period).
*        (see BINPOP for further explanation of SEMI, RATIO & RANGE uses)
*
*          if(ICIRC.gt.1)
*
*       TC      Circularization time.
*
***
* HIPOP: if(kz(8).gt.0.and.kz(18).gt.1)
*
*       NHI     Number of primordial hierarchies.
*       SEMI    Semi-major axis in model units (all equal if RANGE = 0).
*       ECC     Initial eccentricity (< 0 for thermal distribution).
*       RATIO   Mass ratio (= 1.0: M1 = M2; random in [0.5-0.9]).
*       RANGE   Range in SEMI for uniform logarithmic distribution (> 0).
*       ICIRC   Circularization & collision check (not implemented yet).
***
* INTIDE: if(kz(27).gt.0) 
*
*       RSTAR   Size of typical star in S.U..
*       IMS     # idealized main-sequence stars.
*       IEV     # idealized evolved stars.
*       RMS     Scale factor for main-sequence radii (>0: fudge factor).
*       REV     Scale factor for evolved radii (initial size RSTAR).
***
* CLOUD0: if(kz(12).lt.0)
*
*       NCL     Number of interstellar clouds.
*       RB      Boundary radius in pc.
*       VCL     Mean cloud velocity in km/sec.
*       SIGMA   Velocity dispersion (#12 = -2 or -3: Gaussian).
*       DTBIG   Interval for injecting massive cloud (#12 = -2).
*
*       CLM     Cloud mass in solar units (#12 < 0; NCL members).
*       RCL     Cloud radius (softening) in pc.
***
* MOVIE0: if(kz(40).gt.1)
*
*       XMAX    Movie frame half-size.
*       IWIN    Actual window size (in pixels).
*       DTM     Movie time interval (quantized by STEPK).
***
*       ---------------------------------------------------------------------
*
*
*       Options KZ(J)
*       *************
*
*       ---------------------------------------------------------------------
*       1  COMMON save on unit 1 at end of run (=2: every 100*NMAX steps).
*       2  COMMON save on unit 2 at output (=1); restart if DE/E > 5*QE (=2).
*       3  Basic data on unit 3 at main output time (frequency NFIX).
*       4  Binary diagnostics on unit 4 (# threshold levels = KZ(4) < 10).
*       5  Initial conditions (=0: uniform & isotropic; =1: Plummer).
*       6  Output of significant & regularized binaries (=1, 2, 3 & 4).
*       7  Lagrangian radii (>0: RSCALE; =2, 3, 4: output on unit 6 & 7).
*       8  Primordial binaries (>0: BINOUT; >2: BINDAT; > 3: DEGEN & BINEV).
*       9  Individual bodies printed at output time (MIN(5**KZ9,NTOT)).
*      10  Diagnostic KS output (>0: begin; >1: end; >=3: each step).
*    # 11  Synchronization of circular orbits (suppressed)
*      12  Disk shocks (=1: standard model) or interstellar clouds (<0).
*      13  Scaling of time (1: variable by t_cr; 2: variable by t_r;
*                 -1: constant scaling to t_r; -2: constant scaling to t_c).
*      14  External force (=1: standard tidal field; -1: cutoff. )
*      15  Triple, quad, chain (#30 > 0) or merger search (>1: full output).
*      16  Updating of regularization parameters (RMIN, DTMIN & ECLOSE).
*      17  Modification of ETA (>=1) & ETAU (>1) by tolerance QE.
*      18  Hierarchical systems (=1: diagnostics; =2: primordial; =3: both).
*      19  Mass loss (=1: old; =3: Eggleton-Tout; >4: Chernoff-Weinberg).
*      20  Initial mass function (=0,1: Salpeter; >1: various => see IMF).
*      21  Diagnostic output (>0: MODEL; >1: CENTRE; >2: MTRACE; >3: GLOBAL).
*      22  Initial conditions on unit 10 (=1: output; =2,3: input).
*      23  Escaper removal (>1: diagnostics in file ESC; >2: angles unit #6).
*      24  Initial conditions for subsystem (routine SCALE; KZ(24) = #).
*    # 25  Partial reflection of KS binary orbit (GAMMA < GMIN; suppressed).
*      25  HR diagnostics of evolving stars (>0: output; >1: BSs; =3: x&v data).
*      26  Slow-down of two-body motion (=1: KS binary; =2: chain binary).
*      27  Two-body interactions (-2: RADIUS = 0; -1: collision detection;
*          =1: sequential circ; > 0: collision).
*      28  (not used).
*    # 29  Boundary reflection for hot system (suppressed).
*      30  Chain regularization (>=2: main output; >2: diagnostic output).
*      31  Centre of mass correction after energy check.
*      32  Increase of output intervals (based on single particle energy).
*      33  Block-step diagnostics at main output (=2: active pipes).
*    # 34  Roche lobe overflow (suppressed).
*      35  Time offset (global time from TTOT = TIME + DTOFF).
*      36  Step reduction for hierarchical systems (>0: blocks; <0: Hermite).
*      37  Step reduction for encounters with high-velocity particles.
*      38  Multiple use of GRAPE-6 (sleep 1 sec after each timer check).
*      39  Full N neighbour list instead of GRAPE functions. (=1)
*      40  Block-step version (=1; >1: X11 movie).
*       ---------------------------------------------------------------------
* # currently suppressed
*       ---------------------------------------------------------------------
*
*
*       Output counters
*       ***************
*
*       ---------------------------------------------------------------------
*       NBCALL  Neighbour lists (# NBLIST calls).
*       NBESC   Escaped binaries (#23).
*       NBLOCK  Block integration steps.
*       NBKICK  Binary neutron star kicks (#19).
*       NBPREV  Indicator for enforcing new KS scheduling (SUBINT).
*       NBREF   Boundary reflections (#29; suppressed).
*       NBSTAT  Diagnostic data on binary interactions (#4; inactive).
*       NCBLK1  Active binary steps on HARP.
*       NCBLK2  Active binary steps with more than one pipe.
*       NCHAIN  Chain regularizations (#30).
*       NCOLL   Stellar collisions (#27).
*       NDIAG   Diagnostic warning counter (BRAKE, IMPACT & SPIRAL).
*       NDISS   Tidal dissipation at pericentre (#27).
*       NDUMP   Restart counter (STOP after two restarts).
*       NHIVEL  High-velocity search for all particles (#37).
*       NHLIST  Neighbour lists on HARP (not implemented).
*       NIRECT  Initialization of NSTEPI after exceeding 2*10**9.
*       NKICK   Neutron star kicks (#19).
*       NKSHYP  Hyperbolic KS regularizations.
*       NKSMOD  Slow KS motion restarts (#26).
*       NKSPER  Unperturbed KS binary orbits.
*       NKSREF  Partial reflections of KS binary (#25; suppressed).
*       NKSREG  Total KS regularizations.
*       NKSTRY  Two-body regularization attempts.
*       NMERG   Mergers of stable triples or quadruples (#15).
*       NMESC   Escaped mergers (#15 and 23).
*       NMTRY   Attempted mergers.
*       NPRECT  Initialization of NKSPER after exceeding 2*10**9.
*       NPRINT  Output counter (data bank written at NFIX, then reset).
*       NQUAD   Four-body regularizations (#15).
*       NSESC   Escaped single particles (#23).
*       NSHOCK  Tidal shocks (#12).
*       NSHORT  Shortened time-step due to high-velocity particles (#37).
*       NSPERT  Membership of perturber list (set in routine SIEVE).
*       NSTEPC  Chain regularization steps (# DIFSY calls).
*       NSTEPI  Irregular integration steps (reset to zero at 2*10**9).
*       NSTEPQ  Quadruple regularization integration steps (#15).
*       NSTEPS  Total number of time-steps (reset to zero at 2*10**9).
*       NSTEPT  Triple regularization integration steps (#15).
*       NSTEPU  Regularized integration steps.
*       NSYNC   Number of synchronous binaries (e < 0.002; #27).
*       NTIDE   Tidal captures from hyperbolic orbits (#27).
*       NTIMER  Time-step counter (reset to zero on checking CPU time).
*       NTPERT  Perturbation time-scale evaluations (includes NBLIST calls).
*       NTRECT  Time rectifications (not implemented yet).
*       NTRIP   Three-body regularizations (#15).
*       NTTRY   Search for triple, quad & chain regularization or mergers.
*       NURECT  Initialization of NSTEPU after exceeding 2*10**9.
*       NWARN   Warning messages (only first 1000 printed).
*       ---------------------------------------------------------------------
*
*       Counters in COMMON/STAR/
*       ************************
*
*       ---------------------------------------------------------------------
*       NMDOT   Calls to routine MDOT.
*       NMS     Main sequence stars (type 0/1).
*       NHG     Hertzsprung gap (type 2).
*       NRG     Red giants (type 3).
*       NHE     Helium burning (type 4).
*       NRS     red supergiants (type 5).
*       NWD     White dwarfs (type 8).
*       NSN     Neutron stars (type 9).
*       NHI     High-velocity stars (r < 3*<R>).
*       NDD     Doubly generate binaries (type >= 8).
*       NBS     Blue stragglers (type 1).
*       NTZ     Thorne-Zytkow objects (type 9 + 0/1 -> 9).
*       NAS     Accretion-induced supernovae (Roche; M2 > MCH; type 12).
*       NBH     Black holes (type 10).
*       NGB     Gamma-ray bursters.
*       NROCHE  Calls to routine ROCHE.
*       NRO     Completed Roche stages.
*       NBR     Blue Roche (TM < TPHYS & type 0/1).
*       NBRK    Calls to routine BRAKE (GR & MB).
*       NCONT   Contact binaries (Roche; enforced collisions).
*       NCOAL   Coalescence of binaries (Roche; giant collisions).
*       NHYP    Hyperbolic collisions & coalescence.
*       NSPIR   Calls to routine SPIRAL.
*       NSP     Circularization events.
*       NCIRC   Successful circularizations.
*       NSLP    Inactive circularizations (TC > 2x10**9 yr).
*       NEWHI   New hierarchical systems (counted by routine HIARCH).
*       NEINT   Runge-Kutta integration steps for eccentricity modulation.
*       NGLOB0  Initial membership of globular cluster model.
*       NGLOB   Current membership.
*       NSHOCK  Tidal shocks.
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
*       ---------------------------------------------------------------------
*
      RETURN
*
      END
