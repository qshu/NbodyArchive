#include <math.h>

void gpupot(
		int n,
		double m[],
		double x[][3],
		double pot[]) {
#pragma omp parallel for
	for(int i=0; i<n; i++) {
		double poti = 0.0;
		double xi = x[i][0];
		double yi = x[i][1];
		double zi = x[i][2];
#pragma omp simd
		for(int j=0; j<n; j++) {
			double xj = x[j][0];
			double yj = x[j][1];
			double zj = x[j][2];
			double mj  = m[j];
			
			double dx = xj - xi;
			double dy = yj - yi;
			double dz = zj - zi;
			double r2 = dx*dx + dy*dy + dz*dz;
			if (r2 > 0.0) {
				poti += mj/sqrt(r2);
			}
		}
		pot[i] = poti;
	}
}

extern "C"{
	void gpupot_(
			int *n,
			double m[],
			double x[][3],
			double pot[]){
		gpupot(*n, m, x, pot);
	}
}
