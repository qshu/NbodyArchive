#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <vtools.h>
#include <visit.h>
#include <visit_types.h>
#include <lvisit.h>
#include <lvisit_nbody.h>
#define _FORTRANDOUBLEUNDERSCORE
#include <lvisit_nbody_fortran.h>



int flvisit_nbody_init() {
  int rc;
  rc=lvisit_nbody_init();
  return(rc);
}


int flvisit_nbody_check_connection(int *isactive) {
  int rc=0;
  rc=lvisit_nbody_check_connection();
  *isactive=lvisit_is_active();
  return(rc);
}


int flvisit_nbody_close() {
  int rc=0;
  rc=lvisit_nbody_close();
  return(rc);
}


int flvisit_nbody_split() {
  int rc;
  rc=0;
  return(rc);
}


int flvisit_comm_world(int *comm) {
  return(1);
}

int flvisit_nbody_particles_send (double *simtime, double *data_particles_1,double *data_particles_2,double *data_particles_3,double *data_particles_4,double *data_particles_5,double *data_particles_6,double *data_particles_7,int *data_particles_8,int *data_particles_9, int *n1) {
  MPRINTF(("\t\t%-40s: send data of size %d\n","lvisit_nbody_particles_send",*n1));
  lvisit_nbody_particles_send(*simtime,
				data_particles_1,
				data_particles_2,
				data_particles_3,
				data_particles_4,
				data_particles_5,
				data_particles_6,
				data_particles_7,
				data_particles_8,
				data_particles_9,

			   *n1
			   );
  return(1);
}


int flvisit_nbody_controlsim_send (int *VAction,int *nr,int *nr_min,int *nr_max) {
  int rc;
  MPRINTF(("\t\t%-40s: send data of size %d\n","lvisit_nbody_controlsim_send",1));
  {
    controlsim_t data_controlsim;
		data_controlsim.VAction=*VAction;
		data_controlsim.nr=*nr;
		data_controlsim.nr_min=*nr_min;
		data_controlsim.nr_max=*nr_max;

    rc=lvisit_nbody_controlsim_send(&data_controlsim);
  }

  return(rc);
}
int flvisit_nbody_info_send (double *posminx,double *posminy,double *posminz,double *posmaxx,double *posmaxy,double *posmaxz,double *TTOT,double *LOAT,double *BAR,double *MBAR,double *TIDE,double *IDAL,double *DENS1,double *DENS2,double *DENS3,double *TTOT_TCR,double *I6,double *NZERO,double *RC,double *NC,double *VC,double *RHOM,double *CMAX,double *RSCALE,double *RSMIN,double *DMIN1) {
  int rc;
  MPRINTF(("\t\t%-40s: send data of size %d\n","lvisit_nbody_info_send",1));
  {
    info_t data_info;
		data_info.posminx=*posminx;
		data_info.posminy=*posminy;
		data_info.posminz=*posminz;
		data_info.posmaxx=*posmaxx;
		data_info.posmaxy=*posmaxy;
		data_info.posmaxz=*posmaxz;
		data_info.TTOT=*TTOT;
		data_info.LOAT=*LOAT;
		data_info.BAR=*BAR;
		data_info.MBAR=*MBAR;
		data_info.TIDE=*TIDE;
		data_info.IDAL=*IDAL;
		data_info.DENS1=*DENS1;
		data_info.DENS2=*DENS2;
		data_info.DENS3=*DENS3;
		data_info.TTOT_TCR=*TTOT_TCR;
		data_info.I6=*I6;
		data_info.NZERO=*NZERO;
		data_info.RC=*RC;
		data_info.NC=*NC;
		data_info.VC=*VC;
		data_info.RHOM=*RHOM;
		data_info.CMAX=*CMAX;
		data_info.RSCALE=*RSCALE;
		data_info.RSMIN=*RSMIN;
		data_info.DMIN1=*DMIN1;

    rc=lvisit_nbody_info_send(&data_info);
  }

  return(rc);
}


int flvisit_nbody_steering_recv (double *dparm1,double *dparm2,double *dparm3,double *dparm4,int *iparm1,int *iparm2,int *iparm3,int *iparm4) {
  int rc=0;
  steering_t data_steering;

  rc=lvisit_nbody_steering_recv(&data_steering);
  		*dparm1=data_steering.dparm1;
		*dparm2=data_steering.dparm2;
		*dparm3=data_steering.dparm3;
		*dparm4=data_steering.dparm4;
		*iparm1=data_steering.iparm1;
		*iparm2=data_steering.iparm2;
		*iparm3=data_steering.iparm3;
		*iparm4=data_steering.iparm4;

  return(rc);
}
int flvisit_nbody_control_recv (int *VAction,int *nr) {
  int rc=0;
  control_t data_control;

  rc=lvisit_nbody_control_recv(&data_control);
  		*VAction=data_control.VAction;
		*nr=data_control.nr;

  return(rc);
}
int flvisit_nbody_proxyinfo_recv (int *connection_up) {
  int rc=0;
  proxyinfo_t data_proxyinfo;

  rc=lvisit_nbody_proxyinfo_recv(&data_proxyinfo);
  		*connection_up=data_proxyinfo.connection_up;

  return(rc);
}


