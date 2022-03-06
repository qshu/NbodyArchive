#ifndef LVISIT_NBODY_H
#define LVISIT_NBODY_H

#include <lvisit.h>

typedef struct {
   int start;
   int ack;
} init_t;
typedef struct {
   int VAction;
   int nr;
   int nr_min;
   int nr_max;
} controlsim_t;
typedef struct {
   double dparm1;
   double dparm2;
   double dparm3;
   double dparm4;
   int iparm1;
   int iparm2;
   int iparm3;
   int iparm4;
} steering_t;
typedef struct {
   int VAction;
   int nr;
} control_t;
typedef struct {
   int connection_up;
} proxyinfo_t;
typedef struct {
   double posminx;
   double posminy;
   double posminz;
   double posmaxx;
   double posmaxy;
   double posmaxz;
   double TTOT;
   double LOAT;
   double BAR;
   double MBAR;
   double TIDE;
   double IDAL;
   double DENS1;
   double DENS2;
   double DENS3;
   double TTOT_TCR;
   double I6;
   double NZERO;
   double RC;
   double NC;
   double VC;
   double RHOM;
   double CMAX;
   double RSCALE;
   double RSMIN;
   double DMIN1;
} info_t;


int lvisit_nbody_init_uts();

int lvisit_nbody_init();

int lvisit_nbody_check_drift();

int lvisit_nbody_check_connection();

int lvisit_nbody_init_send(init_t *data_init);
int lvisit_nbody_controlsim_send(controlsim_t *data_controlsim);
int lvisit_nbody_info_send(info_t *data_info);


int lvisit_nbody_particles_send(double simtime, double *data_particles_1,double *data_particles_2,double *data_particles_3,double *data_particles_4,double *data_particles_5,double *data_particles_6,double *data_particles_7,int *data_particles_8,int *data_particles_9, int n1);
int lvisit_nbody_particles_send_fld(double *data_particles, int n1, int n2);


int lvisit_nbody_steering_recv(steering_t *data_steering);
int lvisit_nbody_control_recv(control_t *data_control);
int lvisit_nbody_proxyinfo_recv(proxyinfo_t *data_proxyinfo);

int lvisit_nbody_close();

#endif
