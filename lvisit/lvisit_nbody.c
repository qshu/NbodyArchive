#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vtools.h>
#include <visit.h>
#include <visit_types.h>
#include <lvisit.h>
#include <lvisit_nbody.h>
static lvisit_conn_t *conn = NULL;


static visit_type uty_init;
static init_t stdata_init= {
            -1,
            -1,
   };

static visit_type uty_controlsim;

static visit_type uty_steering;

static visit_type uty_control;

static visit_type uty_proxyinfo;
static proxyinfo_t stdata_proxyinfo= {
            -1,
   };

static visit_type uty_info;




int lvisit_nbody_init_uts() {
   {
      init_t init_v;
      int  len_init[2];
      int  off_init[2];
      int subt_init[2];
      len_init[0] = 1;
      len_init[1] = 1;
      off_init[0] = (char *)&init_v.start - (char *)&init_v;
      off_init[1] = (char *)&init_v.ack - (char *)&init_v;
      subt_init[0] = VISIT_INT;
      subt_init[1] = VISIT_INT;
      uty_init = visit_new_type("init",
                                     2, sizeof(init_t),
                                     len_init,
                                     off_init,
                                     subt_init);
   }
   {
      controlsim_t controlsim_v;
      int  len_controlsim[4];
      int  off_controlsim[4];
      int subt_controlsim[4];
      len_controlsim[0] = 1;
      len_controlsim[1] = 1;
      len_controlsim[2] = 1;
      len_controlsim[3] = 1;
      off_controlsim[0] = (char *)&controlsim_v.VAction - (char *)&controlsim_v;
      off_controlsim[1] = (char *)&controlsim_v.nr - (char *)&controlsim_v;
      off_controlsim[2] = (char *)&controlsim_v.nr_min - (char *)&controlsim_v;
      off_controlsim[3] = (char *)&controlsim_v.nr_max - (char *)&controlsim_v;
      subt_controlsim[0] = VISIT_INT;
      subt_controlsim[1] = VISIT_INT;
      subt_controlsim[2] = VISIT_INT;
      subt_controlsim[3] = VISIT_INT;
      uty_controlsim = visit_new_type("controlsim",
                                     4, sizeof(controlsim_t),
                                     len_controlsim,
                                     off_controlsim,
                                     subt_controlsim);
   }
   {
      steering_t steering_v;
      int  len_steering[8];
      int  off_steering[8];
      int subt_steering[8];
      len_steering[0] = 1;
      len_steering[1] = 1;
      len_steering[2] = 1;
      len_steering[3] = 1;
      len_steering[4] = 1;
      len_steering[5] = 1;
      len_steering[6] = 1;
      len_steering[7] = 1;
      off_steering[0] = (char *)&steering_v.dparm1 - (char *)&steering_v;
      off_steering[1] = (char *)&steering_v.dparm2 - (char *)&steering_v;
      off_steering[2] = (char *)&steering_v.dparm3 - (char *)&steering_v;
      off_steering[3] = (char *)&steering_v.dparm4 - (char *)&steering_v;
      off_steering[4] = (char *)&steering_v.iparm1 - (char *)&steering_v;
      off_steering[5] = (char *)&steering_v.iparm2 - (char *)&steering_v;
      off_steering[6] = (char *)&steering_v.iparm3 - (char *)&steering_v;
      off_steering[7] = (char *)&steering_v.iparm4 - (char *)&steering_v;
      subt_steering[0] = VISIT_DOUBLE;
      subt_steering[1] = VISIT_DOUBLE;
      subt_steering[2] = VISIT_DOUBLE;
      subt_steering[3] = VISIT_DOUBLE;
      subt_steering[4] = VISIT_INT;
      subt_steering[5] = VISIT_INT;
      subt_steering[6] = VISIT_INT;
      subt_steering[7] = VISIT_INT;
      uty_steering = visit_new_type("steering",
                                     8, sizeof(steering_t),
                                     len_steering,
                                     off_steering,
                                     subt_steering);
   }
   {
      control_t control_v;
      int  len_control[2];
      int  off_control[2];
      int subt_control[2];
      len_control[0] = 1;
      len_control[1] = 1;
      off_control[0] = (char *)&control_v.VAction - (char *)&control_v;
      off_control[1] = (char *)&control_v.nr - (char *)&control_v;
      subt_control[0] = VISIT_INT;
      subt_control[1] = VISIT_INT;
      uty_control = visit_new_type("control",
                                     2, sizeof(control_t),
                                     len_control,
                                     off_control,
                                     subt_control);
   }
   {
      proxyinfo_t proxyinfo_v;
      int  len_proxyinfo[1];
      int  off_proxyinfo[1];
      int subt_proxyinfo[1];
      len_proxyinfo[0] = 1;
      off_proxyinfo[0] = (char *)&proxyinfo_v.connection_up - (char *)&proxyinfo_v;
      subt_proxyinfo[0] = VISIT_INT;
      uty_proxyinfo = visit_new_type("proxyinfo",
                                     1, sizeof(proxyinfo_t),
                                     len_proxyinfo,
                                     off_proxyinfo,
                                     subt_proxyinfo);
   }
   {
      info_t info_v;
      int  len_info[26];
      int  off_info[26];
      int subt_info[26];
      len_info[0] = 1;
      len_info[1] = 1;
      len_info[2] = 1;
      len_info[3] = 1;
      len_info[4] = 1;
      len_info[5] = 1;
      len_info[6] = 1;
      len_info[7] = 1;
      len_info[8] = 1;
      len_info[9] = 1;
      len_info[10] = 1;
      len_info[11] = 1;
      len_info[12] = 1;
      len_info[13] = 1;
      len_info[14] = 1;
      len_info[15] = 1;
      len_info[16] = 1;
      len_info[17] = 1;
      len_info[18] = 1;
      len_info[19] = 1;
      len_info[20] = 1;
      len_info[21] = 1;
      len_info[22] = 1;
      len_info[23] = 1;
      len_info[24] = 1;
      len_info[25] = 1;
      off_info[0] = (char *)&info_v.posminx - (char *)&info_v;
      off_info[1] = (char *)&info_v.posminy - (char *)&info_v;
      off_info[2] = (char *)&info_v.posminz - (char *)&info_v;
      off_info[3] = (char *)&info_v.posmaxx - (char *)&info_v;
      off_info[4] = (char *)&info_v.posmaxy - (char *)&info_v;
      off_info[5] = (char *)&info_v.posmaxz - (char *)&info_v;
      off_info[6] = (char *)&info_v.TTOT - (char *)&info_v;
      off_info[7] = (char *)&info_v.LOAT - (char *)&info_v;
      off_info[8] = (char *)&info_v.BAR - (char *)&info_v;
      off_info[9] = (char *)&info_v.MBAR - (char *)&info_v;
      off_info[10] = (char *)&info_v.TIDE - (char *)&info_v;
      off_info[11] = (char *)&info_v.IDAL - (char *)&info_v;
      off_info[12] = (char *)&info_v.DENS1 - (char *)&info_v;
      off_info[13] = (char *)&info_v.DENS2 - (char *)&info_v;
      off_info[14] = (char *)&info_v.DENS3 - (char *)&info_v;
      off_info[15] = (char *)&info_v.TTOT_TCR - (char *)&info_v;
      off_info[16] = (char *)&info_v.I6 - (char *)&info_v;
      off_info[17] = (char *)&info_v.NZERO - (char *)&info_v;
      off_info[18] = (char *)&info_v.RC - (char *)&info_v;
      off_info[19] = (char *)&info_v.NC - (char *)&info_v;
      off_info[20] = (char *)&info_v.VC - (char *)&info_v;
      off_info[21] = (char *)&info_v.RHOM - (char *)&info_v;
      off_info[22] = (char *)&info_v.CMAX - (char *)&info_v;
      off_info[23] = (char *)&info_v.RSCALE - (char *)&info_v;
      off_info[24] = (char *)&info_v.RSMIN - (char *)&info_v;
      off_info[25] = (char *)&info_v.DMIN1 - (char *)&info_v;
      subt_info[0] = VISIT_DOUBLE;
      subt_info[1] = VISIT_DOUBLE;
      subt_info[2] = VISIT_DOUBLE;
      subt_info[3] = VISIT_DOUBLE;
      subt_info[4] = VISIT_DOUBLE;
      subt_info[5] = VISIT_DOUBLE;
      subt_info[6] = VISIT_DOUBLE;
      subt_info[7] = VISIT_DOUBLE;
      subt_info[8] = VISIT_DOUBLE;
      subt_info[9] = VISIT_DOUBLE;
      subt_info[10] = VISIT_DOUBLE;
      subt_info[11] = VISIT_DOUBLE;
      subt_info[12] = VISIT_DOUBLE;
      subt_info[13] = VISIT_DOUBLE;
      subt_info[14] = VISIT_DOUBLE;
      subt_info[15] = VISIT_DOUBLE;
      subt_info[16] = VISIT_DOUBLE;
      subt_info[17] = VISIT_DOUBLE;
      subt_info[18] = VISIT_DOUBLE;
      subt_info[19] = VISIT_DOUBLE;
      subt_info[20] = VISIT_DOUBLE;
      subt_info[21] = VISIT_DOUBLE;
      subt_info[22] = VISIT_DOUBLE;
      subt_info[23] = VISIT_DOUBLE;
      subt_info[24] = VISIT_DOUBLE;
      subt_info[25] = VISIT_DOUBLE;
      uty_info = visit_new_type("info",
                                     26, sizeof(info_t),
                                     len_info,
                                     off_info,
                                     subt_info);
   }

  return(0);
}

int lvisit_nbody_init() {
  lvisit_nbody_init_uts();
  conn=lvisit_init("nbody");
  return(0);
}

int lvisit_nbody_check_connection() {
  int changed,lastproxytovisok;
/*    MPRINTF(("\t\t%-40s:\n","lvisit_nbody_check_connection")); */
  changed=lvisit_check_connection();
  lastproxytovisok=conn->proxytovisok;

  if(conn->active) {
    if(lvisit_do_Proxy() && !lvisit_iam_Proxy()) {
      MPRINTF(("\t\t%-40s: there is a proxy between sim. and vis.\n","lvisit_nbody_check_connection"));
      if(lvisit_nbody_proxyinfo_recv(&stdata_proxyinfo)) {
	conn->proxytovisok=stdata_proxyinfo.connection_up;
	MPRINTF(("\t\t%-40s: connection_up=%d between proxy and vis.\n","lvisit_nbody_check_connection",conn->proxytovisok));
      }
    }
  }

  if(lastproxytovisok!=conn->proxytovisok) changed=1;

  if(!conn->active) {return(changed);}
  if((conn->doproxy) && (!conn->iamproxy) && (!conn->proxytovisok)) {return(changed);}

  if(changed) {
    stdata_init.start=1;
    stdata_init.ack=0;
    lvisit_nbody_init_send(&stdata_init);
  }

  return(changed);
}

int lvisit_nbody_check_drift() {
  lvisit_check_drift();
  return(0);
}

int lvisit_nbody_init_send (init_t *data_init) {
  int ret;
  if(!conn->active) {return(-1);}
  if((conn->doproxy) && (!conn->iamproxy) && (!conn->proxytovisok)) {return(-1);}
  if(!conn->SendData) {return(-1);}
  lvisit_timing_setid(1);
  ret=lvisit_sendproto(1,3.1415,data_init,uty_init,1,1,1,1,1,"init");
  if(!lvisit_checkreturncode(ret)) return 0;
  return(1);
}
int lvisit_nbody_controlsim_send (controlsim_t *data_controlsim) {
  int ret;
  if(!conn->active) {return(-1);}
  if((conn->doproxy) && (!conn->iamproxy) && (!conn->proxytovisok)) {return(-1);}
  if(!conn->SendData) {return(-1);}
  lvisit_timing_setid(3);
  ret=lvisit_sendproto(3,3.1415,data_controlsim,uty_controlsim,1,1,1,1,1,"controlsim");
  if(!lvisit_checkreturncode(ret)) return 0;
  return(1);
}
int lvisit_nbody_info_send (info_t *data_info) {
  int ret;
  if(!conn->active) {return(-1);}
  if((conn->doproxy) && (!conn->iamproxy) && (!conn->proxytovisok)) {return(-1);}
  if(!conn->SendData) {return(-1);}
  lvisit_timing_setid(5);
  ret=lvisit_sendproto(5,3.1415,data_info,uty_info,1,1,1,1,1,"info");
  if(!lvisit_checkreturncode(ret)) return 0;
  return(1);
}


int lvisit_nbody_particles_send (double simtime, double *data_particles_1,double *data_particles_2,double *data_particles_3,double *data_particles_4,double *data_particles_5,double *data_particles_6,double *data_particles_7,int *data_particles_8,int *data_particles_9, int n1) {
  int ret,i;
  double before_time;
  static int vecsize=0;
  static double *lvec;
  if(!conn->active) {return(-1);}
  if((conn->doproxy) && (!conn->iamproxy) && (!conn->proxytovisok)) {return(-1);}
  if(!conn->SendData) {return(-1);}

  if(vecsize < n1) {
    if(vecsize > 0) { FREEN(lvec,"lvec");}
    lvec=MALLOCN(n1*9*sizeof(double),"lvec");
    vecsize=n1;
  }

  for(i=0;i<n1;i++) {
    lvec[i*9+1-1]=data_particles_1[i];
    lvec[i*9+2-1]=data_particles_2[i];
    lvec[i*9+3-1]=data_particles_3[i];
    lvec[i*9+4-1]=data_particles_4[i];
    lvec[i*9+5-1]=data_particles_5[i];
    lvec[i*9+6-1]=data_particles_6[i];
    lvec[i*9+7-1]=data_particles_7[i];
    lvec[i*9+8-1]=data_particles_8[i];
    lvec[i*9+9-1]=data_particles_9[i];
    
  }

  lvisit_timing_setid(4);
  ret=lvisit_sendproto(4,simtime,lvec,VISIT_DOUBLE,2,n1,9,1,1,"particles");
  if(!lvisit_checkreturncode(ret)) return 0;
  return(1);
}
int lvisit_nbody_particles_send_fld (double *data_particles, int n1, int n2) {
  int ret;
  unsigned long pnn[4]={1,1,1,1};
  double before_time,time;
  if(!conn->active) {return(-1);}
  if((conn->doproxy) && (!conn->iamproxy) && (!conn->proxytovisok)) {return(-1);}
  if(!conn->SendData) {return(-1);}
  pnn[1-1]=n1;
  pnn[2-1]=n2;
  
  lvisit_timing_setid(4);
  ret=lvisit_sendproto(4,time,data_particles,VISIT_DOUBLE,2,pnn[0],pnn[1],pnn[2],pnn[3],"particles");
  if(!lvisit_checkreturncode(ret)) return 0;
  return(1);
}



int lvisit_nbody_steering_recv (steering_t *data_steering) {
  int id,n1,n2,n3,n4,ndim,ret;
  double time;
  visit_type vtype=uty_steering;
  id=6; n1=1;n2=1,n3=1,n4=1; ndim=1;time=3.14;
  if(!conn->active) {return(-1);}
  if((conn->doproxy) && (!conn->iamproxy) && (!conn->proxytovisok)) {return(-1);}
  if(!conn->RecvData) {return(-1);}
  lvisit_timing_setid(6);
  ret=lvisit_recvproto(id,&time,data_steering,&vtype,&ndim,&n1,&n2,&n3,&n4,"steering");
  if(!lvisit_checkreturncode(ret)) return 0;
  return(1);
}
int lvisit_nbody_control_recv (control_t *data_control) {
  int id,n1,n2,n3,n4,ndim,ret;
  double time;
  visit_type vtype=uty_control;
  id=2; n1=1;n2=1,n3=1,n4=1; ndim=1;time=3.14;
  if(!conn->active) {return(-1);}
  if((conn->doproxy) && (!conn->iamproxy) && (!conn->proxytovisok)) {return(-1);}
  if(!conn->RecvData) {return(-1);}
  lvisit_timing_setid(2);
  ret=lvisit_recvproto(id,&time,data_control,&vtype,&ndim,&n1,&n2,&n3,&n4,"control");
  if(!lvisit_checkreturncode(ret)) return 0;
  return(1);
}

int lvisit_nbody_proxyinfo_recv (proxyinfo_t *data_proxyinfo) {
  int id,n1,n2,n3,n4,ndim,ret;
  double time;
  visit_type vtype=uty_proxyinfo;
  id=29; n1=1;n2=1,n3=1,n4=1; ndim=1;time=3.14;
  if(!conn->active) {return(-1);}
/*    if((conn->doproxy) && (!conn->iamproxy) && (!conn->proxytovisok)) {return(-1);} */
  if(!conn->RecvData) {return(-1);}
  lvisit_timing_setid(29);
  ret=lvisit_recvproto(id,&time,data_proxyinfo,&vtype,&ndim,&n1,&n2,&n3,&n4,"proxyinfo");
  if(!lvisit_checkreturncode(ret)) return 0;
  return(1);
}

int lvisit_nbody_close() {
  MPRINTF(("\t\tlvisit_nbody_close()\n"));
  lvisit_close();
  return(0);
}


