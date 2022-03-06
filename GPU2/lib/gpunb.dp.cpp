#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <algorithm>

#define NNBMAX 768
#define PROFILE

#if 1
#include <omp.h>
#else
static inline int omp_get_num_threads(){return 1;}
static inline int omp_get_thread_num() {return 0;}
#endif


static double *jparr __attribute__((aligned(64))); // {x, y, z, vx, vy, vz, m, pad}
static int cur_nbmax;

static int nbody, nbodymax;

static inline void *amalloc64(const size_t n){
	void *ptr;
	(void)posix_memalign(&ptr, 64, n);
	assert(ptr);
	return ptr;
}

static int num_CPUs = 1;
static inline void GPUNB_open(
		const int nbmax)
{
	nbodymax = nbmax;
	if(nbmax>cur_nbmax){
		if(jparr) free(jparr);
		jparr = (double *)amalloc64(8 * sizeof(double) * (nbmax+15));
		cur_nbmax = nbmax;
#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			if(tid==0) num_CPUs=omp_get_num_threads();
		}
	}
// #ifdef PROFILE
    fprintf(stderr, "# Open AVX regular force - threads: %d\n", num_CPUs);
// #endif
}

static inline void GPUNB_close(){
//	free(jparr); jparr = NULL;
	nbodymax = 0;

// #ifdef PROFILE
 	fprintf(stderr, "# ***********************\n");
 	fprintf(stderr, "# Closed NBODY6/AVX library\n");
 	fprintf(stderr, "# ***********************\n");
// #endif
}

static inline void GPUNB_send(
		const int nj,
		const double mj[],
		const double xj[][3],
		const double vj[][3])
{
//printf("Enter GPUNB_send, nj=%d...\n", nj); fflush(0);
	nbody = nj;
	assert(nbody <= nbodymax);
#pragma omp parallel for
	for(int j=0; j<nj; j++){
		jparr[8*j+0] = xj[j][0];
		jparr[8*j+1] = xj[j][1];
		jparr[8*j+2] = xj[j][2];
		jparr[8*j+3] = vj[j][0];
		jparr[8*j+4] = vj[j][1];
		jparr[8*j+5] = vj[j][2];
		jparr[8*j+6] = mj[j];
	}
//printf("Exit GPUNB_send, t=%g ms\n", (time_send-time_send0)*1e3); fflush(0);
}

static inline void GPUNB_regf(
		const int ni,
		const double h2d[],
		const double dtr[],
		const double xid[][3],
		const double vid[][3],
		double acc[][3],
		double jrk[][3],
		double pot[],
		const int lmax,
		const int nbmax,
		int *listbase)
{
//printf("Enter GPUNB_regf, ni=%d, x1=%g...\n", ni, xid[0][3]); fflush(0);
	assert(nbmax<lmax);

	const int nj = ::nbody;

#pragma omp parallel for
	for(int i=0; i<ni; i++){
		const double xi  = xid[i][0];
		const double yi  = xid[i][1];
		const double zi  = xid[i][2];

		const double vxi = vid[i][0];
		const double vyi = vid[i][1];
		const double vzi = vid[i][2];

		const double h2i = h2d[i];
		const double dtri = dtr[i];

		double poti = 0.0;
		double Ax = 0.0;
		double Ay = 0.0;
		double Az = 0.0;
		double Jx = 0.0;
		double Jy = 0.0;
		double Jz = 0.0;

		int *nnbp = listbase + (lmax * i);
		int nnb = 0;
		
#pragma omp simd
		for(int j=0; j<nj; j++){
			const double xj = jparr[8*j+0];
			const double yj = jparr[8*j+1];
			const double zj = jparr[8*j+2];

			const double vxj = jparr[8*j+3];
			const double vyj = jparr[8*j+4];
			const double vzj = jparr[8*j+5];

			const double dx = xj - xi;
			const double dy = yj - yi;
			const double dz = zj - zi;

			const double r2 = dx*dx + dy*dy + dz*dz;
			if(r2<=0.0) continue;

			const double dvx = vxj - vxi;
			const double dvy = vyj - vyi;
			const double dvz = vzj - vzi;

			const double dxp = dx + dtri * dvx;
			const double dyp = dy + dtri * dvy;
			const double dzp = dz + dtri * dvz;
			const double r2p = dxp*dxp + dyp*dyp + dzp*dzp;

			const double r2min = std::min(r2, r2p);

			const double mj = jparr[8*j+6];

			const double mjh2i = mj * h2i;

			if(r2min<mjh2i){
				++nnb;
				if(nnb<=nbmax){
					nnbp[nnb] = j;
				}
				continue;
			}
			const double rinv = 1.0/sqrt(r2);
			const double rinv2 = rinv*rinv;
			const double rv = -3.0 * (dx*dvx + dy*dvy + dz*dvz) * rinv2;
			const double mrinv3 = mj * rinv * rinv2;

			poti += rinv;

			Jx += (dvx+rv*dx)*mrinv3;
			Jy += (dvy+rv*dy)*mrinv3;
			Jz += (dvz+rv*dz)*mrinv3;

			Ax += dx*mrinv3;
			Ay += dy*mrinv3;
			Az += dz*mrinv3;
		}

		acc[i][0] = Ax;
		acc[i][1] = Ay;
		acc[i][2] = Az;
		jrk[i][0] = Jx;
		jrk[i][1] = Jy;
		jrk[i][2] = Jz;
		pot[i]    = poti;

		nnbp[0] = (nnb > nbmax) ? -1 : nnb;
	}
}

extern "C" {
	void gpunb_open_(int *nbmax){
		GPUNB_open(*nbmax);
	}
	void gpunb_close_(){
		GPUNB_close();
	}
	void gpunb_send_(
			int *nj,
			double mj[],
			double xj[][3],
			double vj[][3]){
		GPUNB_send(*nj, mj, xj, vj);
	}
	void gpunb_regf_(
			int *ni,
			double h2[],
			double dtr[],
			double xi[][3],
			double vi[][3],
			double acc[][3],
			double jrk[][3],
			double pot[],
			int *lmax,
			int *nbmax,
			int *list){ // list[][lmax]
//printf(">gpunb_regf\n");fflush(stdout);
		GPUNB_regf(*ni, h2, dtr, xi, vi, acc, jrk, pot, *lmax, *nbmax, list);
//printf("<gpunb_regf\n");fflush(stdout);
	}
}
