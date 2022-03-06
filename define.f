      SUBROUTINE DEFINE
*
*
*       Definition of input parameters, options & counters. (NB6++)
*       ---------------------------------------------------
*
* ********** NOTE: Contents of [] to be removed with time.*************
* ********** Not read by nb6++                            *************
*
*       Input parameters
*       ****************
*
*       ---------------------------------------------------------------------
*       KSTART  Control index (1: new run; >1: restart; 3, 4, 5: new params).
*       TCOMP   Maximum computing time in minutes (saved in CPU).
*       TCRITp  Termination time in Myrs.
*       isernb  Max size of sequential irr blocks on parallel machine 
*               for single CPU dummy
*       iserreg as isernb for reg blocks
*               for single CPU dummy
*--------
*       N       Total number of centre of masses (<NMAX - 2).
*               e.g: N=100 binaries = NBIN0=NBIN below for f=1
*                    N=100, and NBIN0=NBIN=50 for f=0.5
*       NFIX    Output frequency of data save or binaries (options 3 & 6).
*       NCRIT   Final particle number (alternative termination criterion).
*       NRAND   Random number seed; any positive integer
*      [ NNBMAX  Maximum number of neighbours (= LMAX - 3) set in params.h ].
*       NNBOPT  Desired optimal neighbour number (R.Sp.)
*       NRUN    Run identification index.
*--------
*       ETAI    Time-step parameter for irregular force polynomial.
*       ETAR    Time-step parameter for regular force polynomial.
*       RS0     Initial radius of neighbour sphere.
*       DTADJ   Time interval for parameter adjustment (in TCR if KZ(35)=0).
*       DELTAT  Output time interval (KZ(35) = 0: in TCR; >0: scaled units).
*->             NFIX=1 and DTADJ=DELTAT => OUT3 written every adjust time
*       TCRIT   Termination time in units of TCR if KZ(35)=0.
*->             The _earlier_ termination criterion becomes active
*       QE      Energy tolerance (restart if DE/E > 5*QE & KZ(2) > 1).
*       RBAR    Virial cluster radius in pc (set = 0 for isolated cluster).
*       ZMBAR   Mean mass in solar units (1.0 assumed if 0).
*--------
*       KZ(J)   Non-zero options for alternative paths (see table below).
*       BK(J)   Non-zero options for binpop_4new routine (see table below).
*--------
*       DTMIN   Time-step criterion for regularization search.
*       RMIN    Distance criterion for regularization search.
*       ETAU    Regularized time-step parameter (6.28/ETAU steps/orbit).
*       ECLOSE  Binding energy per unit mass for hard binary (positive).
*       GMIN    Relative two-body perturbation for unperturbed motion.
*       GMAX    Secondary termination parameter for soft KS binaries.
*--------
*         DELTAS  Output interval for binary search (#35=0: in TCR; option 4).
*         ORBITS  Minimum periods for binary output (level 1).
*         GPRINT  Perturbation thresholds for binary output (9 levels).
*--------
*       ALPHAS  Power-law index for initial mass function (routine DATA).
*       BODY1   Maximum particle mass before scaling (KZ(20): solar mass).
*       BODYN   Minimum particle mass before scaling (KZ(20): solar mass).
*       NBIN0   Number of primordial binaries (for IMF2 with KZ(20) > 2).
*       ZMET    Metal abundance (in range 0.03 - 0.0001).
*       EPOCH0  Evolutionary epoch (in 10**6 yrs).
*--------
*       Q       Virial ratio (routine SCALE; Q = 0.5 for equilibrium).
*       VXROT   XY-velocity scaling factor (> 0 for solid-body rotation).
*       VZROT   Z-velocity scaling factor (not used if VXROT = 0).
*       RSPH2   Radius of reflecting sphere (option 29; units of RSCALE).
*--------
*         NCL     Number of interstellar clouds (routine CLOUD0; option 13).
*         RB2     Radius of cloud boundary in pc (square is saved).
*         VCL     Mean cloud velocity in km/sec.
*         SIGMA   Gaussian velocity dispersion of clouds in km/sec.
*         CLM     Individual cloud masses in solar masses (maximum is MCL).
*         RCL2    Half-mass radii of clouds in pc (square is saved).
*--------
*         SIGMA0  Hot initial velocities in km/sec (routine HOTSYS; option 29).
*--------
*       NBIN    Number of initial binaries (routine BINPOP; option 8).
*       SEMI    Initial semi-major axis (= 0 for range of energies).
*       ECC     Initial eccentricity (for BINPOP_4NEW)
*               <=1 AND >=0 for one particular fixed ecc. for all systems
*               < 0 for thermal distribution,
*               =20 for uniform distribution,
*               =30 for f(e)=0.1765/(e*e)
*               =40 for general f(e)=a*e^b, e0<=e<=1 with a=(1+b)/(1-e0^(1+b))
*                   e0 and b must be defined in binpop routine
*       RATIO   Mass ratio M1/(M1 + M2); (= 1.0: M1 = M2 = <M>).
*       NBGR    Number of binaries in fixed energy groups.
*       REDUCE  Reduction factor in semi-major axis for each group.
*       RANGE   Energy range for uniform logarithmic distribution.
*       NSKIP   Binary frequency of mass spectrum (starting from body #1).
*       IDORM   Indicator for dormant binaries (>0: merged components).
*--------
*         RSTAR   Size of typical star in A.U. (routine INTIDE; option 27).
*         RTSTAR  Two-body separation in A.U. for tidal capture.
*         RSYNC   Size of synchronous binary orbit in A.U.
*         EPOCH   Evolutionary epoch (in 10**6 yrs) for mass-radius relation.
*--------
*         IMS     # idealized main-sequence stars (routine INTIDE; option 27).
*         IEV     # idealized evolved stars.
*         RMS     Scale factor for main-sequence radii (initial size RSTAR).
*         REV     Scale factor for evolved radii (initial size RSTAR).
*         IMF     Indicator for realistic stellar radii (scaled by RMS if > 0).
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
*       8  Primordial binaries (=1 & 3; >0: BINOUT; >2: BINDAT).
*       9  Individual bodies printed at output time (MIN(5**KZ9,NTOT)).
*      10  Diagnostic KS output (>0: begin; >1: end; >=3: each step).
*      11  Hierarchical systems (=1: diagnostics; =2: primordial; =3: both).
*      12  HR diagnostics of evolving stars (interval DTPLOT).
*      13  Interstellar clouds (<0 or >2: Gaussian velocities).
*      14  External force (=1: standard tidal field; =2: not implemented).
*      15  Triple, quad, chain (#30 > 0) or merger search (>1: full output).
*      16  Updating of regularization parameters (RMIN, DTMIN & ECLOSE).
*      17  Modification of ETAI, ETAR (>=1) and ETAU (>1) by tolerance QE.
*      18  Neighbour addition in routine CHECKL (=1: LISTV; =2: all types).
*      19  Mass loss (=1: supernova scheme; =3: Eggleton, Tout & Hurley).
*      20  Initial mass function (=1: Tout; =2,4: Kroupa; =3,5: Eggleton).
*          >3 => mass ratio distr. as defined in imf2.f
*          =2 for KTG93 IMF with random pairing (imf2.f)
*      21  Extra output line (MODEL #, TCOMP, DMIN, AMIN, RMAX & RSMIN).
*      22  Initial m,r,v on #10 (=1: output; >=2: input; >2: no scaling).
*                               (=4: starlab input format)
*      23  Removal of escapers (=1: isolated cluster; =2: diag; 
*                               =3: tidal cut       ; =4: diag).
*      24  Initial conditions for subsystem (routine SCALE; KZ(24) = #).
*      25  Partial reflection of KS binary orbit (GAMMA < GMIN; suppressed).
*      26  Slow-down of two-body motion (>=1: KS binary; =2: chain binary).
*      27  Two-body tidal interaction (n = 1.5: type 3 & 5; n = 3: others).
*      28  Magnetic braking and gravitational radiation (options 19 & 27).
*      29  Boundary reflection for hot system (suppressed).
*      30  Chain regularization (>=2: main output; >2: diagnostic output).
*      31  Centre of mass correction after energy check.
*      32  Increase of output intervals (based on single particle energy).
*      33  Block-step diagnostics at main output (>=1: STEP; =2: STEPR).
*      34  Roche lobe overflow (not implemented yet).
*      35  No scaling of output intervals (routine SCALE & MODIFY).
*          =0: DTADJ,DELTAT,TCRIT scaled by TCR (default in nb5)
*      36  Step reduction for hierarchical systems.
*      37  Fast time-step criterion (>0: STEP; >1: STEPR).
*      38  No force polynomial corrections (I <= N; block-step version).
*      39  shape analysis by routine ellan (=2) with Ch. Theis
*      40  adjust neighbour number to optimal neighbour number.
*       ---------------------------------------------------------------------
*
*
*       Options BK(J)   (for binpop_4new.f)
*       *************
*
*       ---------------------------------------------------------------------
*       1  =0: no proto-star evolution                                       
*          =1:"proto-star" evol. of ecc,period in binpop_4new.f           
*       2  = -1: use NBGR and REDUCE in binpop_pk.f
*          =0:flat distr. in semi-major axis                                 
*          =1:f=0.034388logP                                                 
*          =2:f=3.5logP/[100+(logP)**2]                                      
*             KZ(40)=1,2 are 1st and 2nd iterations                          
*          =3:f=2.3(logP-1)/[45+(logP-1)**2]                                 
*          =4:f=2.5(logP-1)/[45+(logP-1)**2] -- derived in K2
*               NOTE: in routine adjust.f KZ(40)>0 is used to adjust         
*                     neighbour number                                       
*          =5:f = Duquennoy&Mayor (1991), i.e. Gaussian in logP 
*       3  =0: only single precission output on "unit3" (i.e. "OUT3" in nb5)
*              written in fort.3000+rank.
*          =1: only double precission output on "unit 3",
*              but written in fort.4000+rank. 
*              Necessary for full binary star analysis.
*          =2: both, sngl & dbl, files above are written
*       4  =1: file peri_hyperbol.dat opened and written to (see ksint.f)
*          =0: not opened
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
*       NFAST   Fast particle included in LISTV (option 18).
*       NBFAST  Fast particle included in neighbour list (option 18).
*       NBREF   Boundary reflections (option 29; suppressed).
*       NBLOCK  Number of irregular blocks (block-step version).
*       NBLCKR  Number of regular blocks (block-step version) (R.Sp.)
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
*       NSTEPT  Triple regularization integration steps (option 15).
*       NSTEPQ  Quadruple regularization integration steps (option 15).
*       NSTEPC  Chain regularization steps (# DIFSY calls).
*       NDISS   Tidal dissipation at pericentre (option 27).
*       NTIDE   Tidal captures from hyperbolic motion (option 27).
*       NSYNC   Number of synchronous binaries (a < RSYNC; option 27).
*       NCOLL   Stellar collisions (option 27).
*       NSESC   Escaped single particles (option 23).
*       NBESC   Escaped binaries (option 23).
*       NMESC   Escaped mergers (options 15 & 23).
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
*      20       Circularized binary (c.m. value).
*       ---------------------------------------------------------------------
*
      RETURN
*
      END
