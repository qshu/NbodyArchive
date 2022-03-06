#ifndef LVISIT_NBODY_FORTRAN_H
#define LVISIT_NBODY_FORTRAN_H
#include <lvisit_nbody.h>


#if defined(_FORTRANCAPS)
#define flvisit_nbody_init FLVISIT_NBODY_INIT
#elif defined(_FORTRANUNDERSCORE)
#define flvisit_nbody_init flvisit_nbody_init_
#elif defined(_FORTRANDOUBLEUNDERSCORE)
#define flvisit_nbody_init flvisit_nbody_init__
#elif defined(_FORTRANNOUNDERSCORE)
#else
#error nothing defined for fortran externals
#endif
int flvisit_nbody_init();


#if defined(_FORTRANCAPS)
#define flvisit_nbody_check_connection FLVISIT_NBODY_CHECK_CONNECTION
#elif defined(_FORTRANUNDERSCORE)
#define flvisit_nbody_check_connection flvisit_nbody_check_connection_
#elif defined(_FORTRANDOUBLEUNDERSCORE)
#define flvisit_nbody_check_connection flvisit_nbody_check_connection__
#elif defined(_FORTRANNOUNDERSCORE)
#else
#error nothing defined for fortran externals
#endif
int flvisit_nbody_check_connection(int *isactive);


#if defined(_FORTRANCAPS)
#define flvisit_nbody_close FLVISIT_NBODY_CLOSE
#elif defined(_FORTRANUNDERSCORE)
#define flvisit_nbody_close flvisit_nbody_close_
#elif defined(_FORTRANDOUBLEUNDERSCORE)
#define flvisit_nbody_close flvisit_nbody_close__
#elif defined(_FORTRANNOUNDERSCORE)
#else
#error nothing defined for fortran externals
#endif
int flvisit_nbody_close();


#if defined(_FORTRANCAPS)
#define flvisit_nbody_split FLVISIT_NBODY_SPLIT
#elif defined(_FORTRANUNDERSCORE)
#define flvisit_nbody_split flvisit_nbody_split_
#elif defined(_FORTRANDOUBLEUNDERSCORE)
#define flvisit_nbody_split flvisit_nbody_split__
#elif defined(_FORTRANNOUNDERSCORE)
#else
#error nothing defined for fortran externals
#endif
int flvisit_nbody_split();


#if defined(_FORTRANCAPS)
#define flvisit_comm_world FLVISIT_COMM_WORLD
#elif defined(_FORTRANUNDERSCORE)
#define flvisit_comm_world flvisit_comm_world_
#elif defined(_FORTRANDOUBLEUNDERSCORE)
#define flvisit_comm_world flvisit_comm_world__
#elif defined(_FORTRANNOUNDERSCORE)
#else
#error nothing defined for fortran externals
#endif
int flvisit_comm_world(int *comm);


#if defined(_FORTRANCAPS)
#define flvisit_nbody_particles_send FLVISIT_NBODY_PARTICLES_SEND
#elif defined(_FORTRANUNDERSCORE)
#define flvisit_nbody_particles_send flvisit_nbody_particles_send_
#elif defined(_FORTRANDOUBLEUNDERSCORE)
#define flvisit_nbody_particles_send flvisit_nbody_particles_send__
#elif defined(_FORTRANNOUNDERSCORE)
#else
#error nothing defined for fortran externals
#endif
int flvisit_nbody_particles_send(double *simtime, double *data_particles_1,double *data_particles_2,double *data_particles_3,double *data_particles_4,double *data_particles_5,double *data_particles_6,double *data_particles_7,int *data_particles_8,int *data_particles_9, int *n1);

#if defined(_FORTRANCAPS)
#define flvisit_nbody_controlsim_send FLVISIT_NBODY_CONTROLSIM_SEND
#elif defined(_FORTRANUNDERSCORE)
#define flvisit_nbody_controlsim_send flvisit_nbody_controlsim_send_
#elif defined(_FORTRANDOUBLEUNDERSCORE)
#define flvisit_nbody_controlsim_send flvisit_nbody_controlsim_send__
#elif defined(_FORTRANNOUNDERSCORE)
#else
#error nothing defined for fortran externals
#endif
int flvisit_nbody_controlsim_send(int *VAction,int *nr,int *nr_min,int *nr_max);
#if defined(_FORTRANCAPS)
#define flvisit_nbody_info_send FLVISIT_NBODY_INFO_SEND
#elif defined(_FORTRANUNDERSCORE)
#define flvisit_nbody_info_send flvisit_nbody_info_send_
#elif defined(_FORTRANDOUBLEUNDERSCORE)
#define flvisit_nbody_info_send flvisit_nbody_info_send__
#elif defined(_FORTRANNOUNDERSCORE)
#else
#error nothing defined for fortran externals
#endif
int flvisit_nbody_info_send(double *posminx,double *posminy,double *posminz,double *posmaxx,double *posmaxy,double *posmaxz,double *TTOT,double *LOAT,double *BAR,double *MBAR,double *TIDE,double *IDAL,double *DENS1,double *DENS2,double *DENS3,double *TTOT_TCR,double *I6,double *NZERO,double *RC,double *NC,double *VC,double *RHOM,double *CMAX,double *RSCALE,double *RSMIN,double *DMIN1);


#if defined(_FORTRANCAPS)
#define flvisit_nbody_steering_recv FLVISIT_NBODY_STEERING_RECV
#elif defined(_FORTRANUNDERSCORE)
#define flvisit_nbody_steering_recv flvisit_nbody_steering_recv_
#elif defined(_FORTRANDOUBLEUNDERSCORE)
#define flvisit_nbody_steering_recv flvisit_nbody_steering_recv__
#elif defined(_FORTRANNOUNDERSCORE)
#else
#error nothing defined for fortran externals
#endif
int flvisit_nbody_steering_recv(double *dparm1,double *dparm2,double *dparm3,double *dparm4,int *iparm1,int *iparm2,int *iparm3,int *iparm4);
#if defined(_FORTRANCAPS)
#define flvisit_nbody_control_recv FLVISIT_NBODY_CONTROL_RECV
#elif defined(_FORTRANUNDERSCORE)
#define flvisit_nbody_control_recv flvisit_nbody_control_recv_
#elif defined(_FORTRANDOUBLEUNDERSCORE)
#define flvisit_nbody_control_recv flvisit_nbody_control_recv__
#elif defined(_FORTRANNOUNDERSCORE)
#else
#error nothing defined for fortran externals
#endif
int flvisit_nbody_control_recv(int *VAction,int *nr);
#if defined(_FORTRANCAPS)
#define flvisit_nbody_proxyinfo_recv FLVISIT_NBODY_PROXYINFO_RECV
#elif defined(_FORTRANUNDERSCORE)
#define flvisit_nbody_proxyinfo_recv flvisit_nbody_proxyinfo_recv_
#elif defined(_FORTRANDOUBLEUNDERSCORE)
#define flvisit_nbody_proxyinfo_recv flvisit_nbody_proxyinfo_recv__
#elif defined(_FORTRANNOUNDERSCORE)
#else
#error nothing defined for fortran externals
#endif
int flvisit_nbody_proxyinfo_recv(int *connection_up);


#endif
