#include <algorithm>
#include <cmath>
#include <cassert>

static inline void cnbint(
		int i,
		const double pos[][3],
		const double vel[][3],
		const double mass[],
		int nnb,
		int list[],
		double f[3],
		double fdot[3]){

	const double xi = pos[i][0];
	const double yi = pos[i][1];
	const double zi = pos[i][2];
	const double vxi = vel[i][0];
	const double vyi = vel[i][1];
	const double vzi = vel[i][2];

	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;
	double jx = 0.0;
	double jy = 0.0;
	double jz = 0.0;

	for(int k=0; k<nnb; k++){
		const int j = list[k];
		assert(j != i);
		const double xj = pos[j][0];
		const double yj = pos[j][1];
		const double zj = pos[j][2];
		const double vxj = vel[j][0];
		const double vyj = vel[j][1];
		const double vzj = vel[j][2];
		const double mj = mass[j];
		const double dx = xj - xi;
		const double dy = yj - yi;
		const double dz = zj - zi;
		const double dvx = vxj - vxi;
		const double dvy = vyj - vyi;
		const double dvz = vzj - vzi;
		
		const double r2 = dx*dx + dy*dy + dz*dz;
		const double rinv = 1.0/sqrt(r2);
		const double rinv2 = rinv * rinv;
		const double rv = -3.0*rinv2*(dx*dvx + dy*dvy + dz*dvz);
		const double mrinv3 = mj * rinv * rinv2;
			 
		ax += dx*mrinv3;
		ay += dy*mrinv3;
		az += dz*mrinv3;
		jx += (dvx+rv*dx)*mrinv3;
		jy += (dvy+rv*dy)*mrinv3;
		jz += (dvz+rv*dz)*mrinv3;
	}
	f[0] = ax;
	f[1] = ay;
	f[2] = az;
	fdot[0] = jx;
	fdot[1] = jy;
	fdot[2] = jz;
	assert(f[0] == f[0]);
	assert(f[1] == f[1]);
	assert(f[2] == f[2]);
	assert(fdot[0] == fdot[0]);
	assert(fdot[1] == fdot[1]);
	assert(fdot[2] == fdot[2]);
}

extern "C" {
	void cnbint_(
		int *i,
		double pos[][3],
		double vel[][3],
		double mass[],
		int *nnb,
		int *nblist,
		double f[3],
		double fdot[3]){
		cnbint(*i, pos-1, vel-1, mass-1, *nnb-1, nblist, f, fdot);
	}
}
