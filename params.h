*       NBODY1 parameters.
*       ------------------
*
      integer NMAX
      PARAMETER  (NMAX=300022)
*     PARAMETER  (NMAX=10022)
*
*
*       ----------------------------------------
*       NMAX    Maximum number of bodies.
*       ----------------------------------------
*
      real machine_eps, meps2
*      parameter (machine_eps=1e-5)
      parameter (machine_eps=3e-3)
      parameter (meps2=machine_eps*machine_eps)
