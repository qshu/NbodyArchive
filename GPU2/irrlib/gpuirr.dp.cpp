#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <omp.h>
#include <stddef.h>
#include <new>
#include <cmath>

#define NNBMAX 768

#define _out_
//#define PROFILE

struct Particle{
	double time;
	double jrk6[3];
	double acc2[3];
	double vel[3];
	double pos[3];
	double mass;
	double pad[2];

	void make_v3(const double x[3], double y[3]) const {
		y[0] = x[0];
		y[1] = x[1];
		y[2] = x[2];
	}
	Particle(
		const double _pos [3],
		const double _vel [3],
		const double _acc2[3],
		const double _jrk6[3],
		const double _mass,
		const double _time)
	{
		time = _time;
		make_v3(_jrk6, jrk6);
		make_v3(_acc2, acc2);
		make_v3(_vel, vel);
		make_v3(_pos, pos);
		mass = _mass;
	}

	Particle(int) // constructor for a dummy particle
	{
		time = 0.0;
		jrk6[0] = jrk6[1] = jrk6[2] =
			acc2[0] = acc2[1] = acc2[2] =
			vel[0] = vel[1] = vel[2] = 0.0;
		pos[0] = pos[1] = pos[2] = 1e6;
		mass = 0.0;
	}
};

struct NBlist{
	int pad[15];
	int nnb;
	int nb[16*((NNBMAX+15)/16)];

	NBlist(){
		nnb = 0;
		for(int i=0; i<NNBMAX; i++){
			nb[i] = 0;
		}
	}

	void print(const int i, FILE *fp = stdout) const{
		fprintf(fp, "%6d%6d :", i, nnb);
		for(int k=0; k<nnb; k++){
			fprintf(fp, " %d", nb[k]);
		}
		fprintf(fp, "\n");
		fflush(fp);
	}
};

// The irregular force library
static Particle  *ptcl = NULL;
static NBlist    *list = NULL;
//__assume_aligned(ptcl, 64);
//__assume_aligned(list, 64);
static int        nmax;
static double     tnow;

static inline void gpuirr_open(
		const int nmax,
		const int lmax)
{
	assert(lmax <= 1 + NNBMAX);

	fprintf(stderr, "# Opening gpuirr lib. AVX ver. - nmax: %d, lmax: %d\n",nmax,lmax); 
	assert(0 == sizeof(NBlist)%16);

	void *ptr;
	int res;
	res = posix_memalign(&ptr, 64, (1+nmax) * sizeof(Particle));
	assert( 0 == res );
	memset(ptr, 0xff, (1+nmax) * sizeof(Particle));
	ptcl = (Particle *)ptr;

	// store dummy predictor to the last of the array
	for(int i=0; i<=nmax; ++i)
		new (&ptcl[i]) Particle(0);

	res = posix_memalign(&ptr, 64, nmax * sizeof(NBlist));
	assert( 0 == res );
	memset(ptr, 0xff, nmax * sizeof(NBlist));
	list = (NBlist *)ptr;

	::nmax = nmax;
}

static inline void gpuirr_close()
{
	free(ptcl); ptcl = NULL;
	free(list); list = NULL;

//	fprintf(stderr, "# Closing gpuirr lib. CPU ver.: grav time = %f sec\n", time_grav);
}

static inline void gpuirr_set_jp(
		const int addr,
		const double pos [3],
		const double vel [3],
		const double acc2[3],
		const double jrk6[3],
		const double mass,
		const double time)
{
//	__assume_aligned(ptcl, 64);
//	__assume_aligned(list, 64);

	new (&ptcl[addr]) Particle(pos, vel, acc2, jrk6, mass, time);
}

static inline void gpuirr_set_list(
		const int addr,
		const int nnb,
		const int nblist[])
{
//	__assume_aligned(ptcl, 64);
//	__assume_aligned(list, 64);

	assert(nnb <= NNBMAX);
	list[addr].nnb = nnb;

	int *dst = list[addr].nb;
	for(int k=0; k<nnb; k++){
		dst[k] = nblist[k];
	}
	// fill dummy
	const int kmax = 16 * ( (nnb+15)/16 );
	for(int k=nnb; k<kmax; k++){
		dst[k] = 0;
	}
}

static inline void gpuirr_pred_all(
		const int    js,
		const int    je,
		const double ti)
{
	::tnow = ti;
}

static inline void gpuirr_pred_act(
		const int ni,
		const int addr[],
		const double ti)
{
	::tnow = ti;
}

static inline void gpuirr_firr(
		const int addr,
		_out_ double accout[3],
		_out_ double jrkout[3])
{
//	__assume_aligned(ptcl, 64);
//	__assume_aligned(list, 64);

	const Particle *ptcl = ::ptcl;
	const Particle *pri = ptcl + addr;
	const NBlist   *nbi = list + addr;
	const int  nnb = nbi->nnb;

	// predict current position and velocity
	const double dt = tnow - pri->time;

	const double posxi = pri->pos[0] + dt*(pri->vel[0] + dt*(pri->acc2[0] + dt*pri->jrk6[0]));
	const double posyi = pri->pos[1] + dt*(pri->vel[1] + dt*(pri->acc2[1] + dt*pri->jrk6[1]));
	const double poszi = pri->pos[2] + dt*(pri->vel[2] + dt*(pri->acc2[2] + dt*pri->jrk6[2]));

	const double velxi = pri->vel[0] + dt*(2.0 * pri->acc2[0] + 3.0 * dt*pri->jrk6[0]);
	const double velyi = pri->vel[1] + dt*(2.0 * pri->acc2[1] + 3.0 * dt*pri->jrk6[1]);
	const double velzi = pri->vel[2] + dt*(2.0 * pri->acc2[2] + 3.0 * dt*pri->jrk6[2]);

	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;
	double jx = 0.0;
	double jy = 0.0;
	double jz = 0.0;

	// sum forces
#pragma omp simd
	for (int j = 0; j < nnb; j++) {
		const int offset = nbi->nb[j];
		const Particle *prj = ptcl + offset;

		const double dtj = tnow - prj->time;

		// get position and velocity

		const double jrkxj = prj->jrk6[0];
		const double accxj = prj->acc2[0];
		const double posxj = ((dtj*jrkxj + accxj)*dtj + prj->vel[0])*dtj +  prj->pos[0];
		const double velxj = (3.0*dtj*jrkxj + 2.0*accxj)*dtj + prj->vel[0];
		const double dx = posxj - posxi;
		const double dvx = velxj - velxi;

		const double jrkyj = prj->jrk6[1];
		const double accyj = prj->acc2[1];
		const double posyj = ((dtj*jrkyj + accyj)*dtj + prj->vel[1])*dtj +  prj->pos[1];
		const double velyj = (3.0*dtj*jrkyj + 2.0*accyj)*dtj + prj->vel[1];
		const double dy = posyj - posyi;
		const double dvy = velyj - velyi;

		const double jrkzj = prj->jrk6[2];
		const double acczj = prj->acc2[2];
		const double poszj = ((dtj*jrkzj + acczj)*dtj + prj->vel[2])*dtj +  prj->pos[2];
		const double velzj = (3.0*dtj*jrkzj + 2.0*acczj)*dtj + prj->vel[2];
		const double dz = poszj - poszi;
		const double dvz = velzj - velzi;

		// calc force
		const double r2 = dx*dx + dy*dy + dz*dz;

		if (r2<=1e-50) {
			printf("IRR0: %d %d\n",addr,offset);fflush(0);
		}		
		const double rinv = 1.0/sqrt(r2);
		const double rinv2 = rinv*rinv;
		const double rv = -3.0*(dx*dvx + dy*dvy + dz*dvz)*rinv2;
		const double mrinv3 = prj->mass * rinv * rinv2;

		ax += dx*mrinv3;
		ay += dy*mrinv3;
		az += dz*mrinv3;

		jx += (dvx+rv*dx)*mrinv3;
		jy += (dvy+rv*dy)*mrinv3;
		jz += (dvz+rv*dz)*mrinv3;
	}

	accout[0] = ax;
	accout[1] = ay;
	accout[2] = az;
	jrkout[0] = jx;
	jrkout[1] = jy;
	jrkout[2] = jz;
}

static inline void gpuirr_firr_vec(
		const int ni,
		const int addr[],
		_out_ double accout[][3],
		_out_ double jrkout[][3])
{
#pragma omp parallel for schedule(dynamic)
	for(int i=0; i<ni; i++){
		gpuirr_firr(addr[i], accout[i], jrkout[i]);
	}
}

// FORTRAN interface
extern "C"{
    void gpuirr_open_(int *nmax, int *lmax) {
        gpuirr_open(*nmax, *lmax);
	}
	void gpuirr_close_() {
		gpuirr_close();
	}
    
	void gpuirr_set_jp_(
		int    *addr,
		double  pos [3],
		double  vel [3],
		double  acc2[3],
		double  jrk6[3],
		double *mass,
		double *time)
	{
		gpuirr_set_jp(*addr, pos, vel, acc2, jrk6, *mass, *time);
	}

	void gpuirr_set_list_(
		int *addr,
		int *nblist)
	{
		gpuirr_set_list(*addr, *nblist, nblist+1);
	}

	void gpuirr_pred_all_(
		int    *js,
		int    *je,
		double *ti)
	{
		gpuirr_pred_all(*js, *je,  *ti);
	}

	void gpuirr_pred_act_(
		int    *ni,
		int     addr[],
		double *ti)
	{
		gpuirr_pred_act(*ni, addr, *ti);
	}

	void gpuirr_firr_vec_(
			int   *ni,
			int    addr[],
			double acc [][3],
			double jrk [][3])
	{
		gpuirr_firr_vec(*ni, addr, acc, jrk);
	}
}
