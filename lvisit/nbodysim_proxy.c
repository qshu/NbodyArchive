#include <stdlib.h>    
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include "vtools.h"
#include "visit_srv.h"
#include "visit_types.h"
#include <lvisit_nbody.h>
#include <lvisit_srv.h>







#define CONN_TIMEOUT -1
#define VISIT_MAXDIM 4
/* program options */
int parm_iter;
int parm_startnr;
int parm_endnr;
int parm_factor; /* scaling of fields, not necessary for server */
int loop=0;
char helpstr1[255],helpstr2[255],helpstr3[255];

void usage();

char* infilename;
void scan_program_options(int argc,char **argv);


/* --------------------------------------------------------------- */
int main(int argc,char **argv) {
  int vscd,id,size,ret1,ret2;
  visit_request *req;
  double local_starttime,remote_starttime,StartTime,RemoteStartTime,LastTime;
  double delta,local_t0, local_t1, remote_t0;

   init_t data_init = {
  -1,-1}
  ;
   controlsim_t data_controlsim = {
  -1,-1,-1,-1}
  ;
   steering_t data_steering = {
  -1.0,-1.0,-1.0,-1.0,-1,-1,-1,-1}
  ;
   control_t data_control = {
  -1,-1}
  ;
   proxyinfo_t data_proxyinfo = {
  -1}
  ;
   info_t data_info = {
  -1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0}
  ;
  
   int esize_particles;
   int data_particles_size=0;
   double *data_particles=NULL;
  


  req = visit_new_request(VISIT_MAXDIM); /* max. dim */
 
  parm_iter=10;
  infilename="nbodytest";
  parm_startnr=1;
  parm_endnr=  100;
  
  
  scan_program_options(argc,argv);

  lvisit_nbody_init();
  lvisit_set_iam_Proxy();

 LOOP:
  
  if( (vscd = visit_srv_init_socket(lvisit_get_Proxy_service(), lvisit_get_Proxy_passwd(), "*", 0,
                                    VISIT_SRV_SEAP_TOGGLE,
                                    CONN_TIMEOUT,
                                    NULL, NULL)) < 0 ) {
    fprintf(stderr, "vserv: init failed\n");
    exit(1);
  }
  
  printf("----- starting nbodysim_proxy\n");
  
  if( visit_srv_connect(vscd, NULL, NULL) < 0) {
    fprintf(stderr, "vserv: connect failed\n");
    goto VERROR;
  }
  

  lvisit_timing_init(visit_srv_basetime(vscd),"PRO");
  lvisit_nbody_init_uts();

  /* start connection to visualization */
  lvisit_nbody_check_connection();

  /* forever loop */
  for(;;) {
    switch( visit_srv_get_id(vscd,&id) ) {
    case -1:
      printf("visit_srv_get_id: error\n");
      goto VERROR;
      break;
    case 0:
      /* nothing to do (e.g. drift-exchange) */
      printf("visit_srv_get_id: drift change\n");
      continue;
      break;
    default:
      break;
    }
    printf("+=======================================================+\n");
    printf("| ID %10d ...                                     |\n",id);
    printf("+=======================================================+\n");
    lvisit_timing_setid(id);
    visit_srv_timings(vscd, &local_t0, &local_t1, &remote_t0);
    lvisit_timing_start_time("get_id",local_t0+0.000001);
    lvisit_timing_end("get_id");

    lvisit_timing_start("get_request");
    if( (size=visit_srv_get_request(vscd,req)) < 0) goto VERROR;
    lvisit_timing_end("get_request");
    printf("\tgot_request id = %d size = %d dir = %d\n",req->id,size,req->dir);
    
     if(id==6) { /* steering */
    
       /* send data to visualization */
       lvisit_nbody_steering_recv(&data_steering);
    
       lvisit_timing_start("write_data");
       if(! (ret1=visit_srv_write_data(vscd, &data_steering, req->vtype, req)) ) goto VERROR;
       lvisit_timing_end_time_abs("write_data",vtools_dtime(0.0)-0.000003);
       visit_srv_timings(vscd, &local_t0, &local_t1, &remote_t0);
       lvisit_timing_start("send_ack2");
       if(! (ret2=visit_srv_ack2(vscd))) goto VERROR;
       lvisit_timing_end("send_ack2");
    
       lvisit_timing_start_time("send",local_t0);
       printf("\t\tSEND -> from=%02d id=%02d group=%-15s size=%-8d Bytes n=%-14s ret=%d,%d bw=%7.4f MB/s (%s)\n",
    	  visit_srv_rank(vscd),id,"steering",size,
    	  printcoord(req->n[0],req->n[1],req->n[2],req->n[3],helpstr1),ret1,ret2,size/(delta*1024*1024),
    	  visit_vtype2name(req->vtype));
       lvisit_timing_end_time("send",local_t1);
    
     }
     if(id==2) { /* control */
    
       /* send data to visualization */
       lvisit_nbody_control_recv(&data_control);
    
       lvisit_timing_start("write_data");
       if(! (ret1=visit_srv_write_data(vscd, &data_control, req->vtype, req)) ) goto VERROR;
       lvisit_timing_end_time_abs("write_data",vtools_dtime(0.0)-0.000003);
       visit_srv_timings(vscd, &local_t0, &local_t1, &remote_t0);
       lvisit_timing_start("send_ack2");
       if(! (ret2=visit_srv_ack2(vscd))) goto VERROR;
       lvisit_timing_end("send_ack2");
    
       lvisit_timing_start_time("send",local_t0);
       printf("\t\tSEND -> from=%02d id=%02d group=%-15s size=%-8d Bytes n=%-14s ret=%d,%d bw=%7.4f MB/s (%s)\n",
    	  visit_srv_rank(vscd),id,"control",size,
    	  printcoord(req->n[0],req->n[1],req->n[2],req->n[3],helpstr1),ret1,ret2,size/(delta*1024*1024),
    	  visit_vtype2name(req->vtype));
       lvisit_timing_end_time("send",local_t1);
    
     }
    
     if(id==29) { /* proxyinfo */
    
       /* send data to visualization */
       lvisit_nbody_check_connection();
       data_proxyinfo.connection_up=lvisit_is_active();
    
       lvisit_timing_start("write_data");
       if(! (ret1=visit_srv_write_data(vscd, &data_proxyinfo, req->vtype, req)) ) goto VERROR;
       lvisit_timing_end_time_abs("write_data",vtools_dtime(0.0)-0.000003);
       visit_srv_timings(vscd, &local_t0, &local_t1, &remote_t0);
       lvisit_timing_start("send_ack2");
       if(! (ret2=visit_srv_ack2(vscd))) goto VERROR;
       lvisit_timing_end("send_ack2");
    
       lvisit_timing_start_time("send",local_t0);
       printf("\t\tSEND -> from=%02d id=%02d group=%-15s size=%-8d Bytes n=%-14s ret=%d,%d bw=%7.4f MB/s (%s)\n",
    	  visit_srv_rank(vscd),id,"proxyinfo",size,
    	  printcoord(req->n[0],req->n[1],req->n[2],req->n[3],helpstr1),ret1,ret2,size/(delta*1024*1024),
    	  visit_vtype2name(req->vtype));
       lvisit_timing_end_time("send",local_t1);
    
     }
    
    
     if(id==1) { /* init */
       lvisit_timing_start("read_data");
       if(! (ret1=visit_srv_read_data(vscd, &data_init, req->vtype, req)) ) goto VERROR;
       lvisit_timing_end_time_abs("read_data",vtools_dtime(0.0)-0.000003);
       visit_srv_timings(vscd, &local_t0, &local_t1, &remote_t0);
       lvisit_timing_start("send_ack2");
       if(! (ret2=visit_srv_ack2(vscd))) goto VERROR;
       lvisit_timing_end("send_ack2");
       printf("\t\tinit: start=%d, ack=%d\n", data_init.start, data_init.ack);
       
       
    
       delta=local_t1-local_t0;if (delta==0) delta=-1.0;
      
       lvisit_timing_start_time("receive",local_t0);
       printf("\t\tRECV -> from=%02d id=%02d group=%-15s size=%-8d Bytes n=%-14s ret=%d,%d bw=%7.4f MB/s (%s)\n",
    	  visit_srv_rank(vscd),id,"init",size,
    	  printcoord(req->n[0],req->n[1],req->n[2],req->n[3],helpstr1),ret1,ret2,size/(delta*1024*1024),
    	  visit_vtype2name(req->vtype));
       lvisit_timing_end_time("receive",local_t1);
    
       /* send data to visualization */
       lvisit_nbody_init_send(&data_init);
     }
     if(id==3) { /* controlsim */
       lvisit_timing_start("read_data");
       if(! (ret1=visit_srv_read_data(vscd, &data_controlsim, req->vtype, req)) ) goto VERROR;
       lvisit_timing_end_time_abs("read_data",vtools_dtime(0.0)-0.000003);
       visit_srv_timings(vscd, &local_t0, &local_t1, &remote_t0);
       lvisit_timing_start("send_ack2");
       if(! (ret2=visit_srv_ack2(vscd))) goto VERROR;
       lvisit_timing_end("send_ack2");
       
       
    
       delta=local_t1-local_t0;if (delta==0) delta=-1.0;
      
       lvisit_timing_start_time("receive",local_t0);
       printf("\t\tRECV -> from=%02d id=%02d group=%-15s size=%-8d Bytes n=%-14s ret=%d,%d bw=%7.4f MB/s (%s)\n",
    	  visit_srv_rank(vscd),id,"controlsim",size,
    	  printcoord(req->n[0],req->n[1],req->n[2],req->n[3],helpstr1),ret1,ret2,size/(delta*1024*1024),
    	  visit_vtype2name(req->vtype));
       lvisit_timing_end_time("receive",local_t1);
    
       /* send data to visualization */
       lvisit_nbody_controlsim_send(&data_controlsim);
     }
     if(id==5) { /* info */
       lvisit_timing_start("read_data");
       if(! (ret1=visit_srv_read_data(vscd, &data_info, req->vtype, req)) ) goto VERROR;
       lvisit_timing_end_time_abs("read_data",vtools_dtime(0.0)-0.000003);
       visit_srv_timings(vscd, &local_t0, &local_t1, &remote_t0);
       lvisit_timing_start("send_ack2");
       if(! (ret2=visit_srv_ack2(vscd))) goto VERROR;
       lvisit_timing_end("send_ack2");
       
       
    
       delta=local_t1-local_t0;if (delta==0) delta=-1.0;
      
       lvisit_timing_start_time("receive",local_t0);
       printf("\t\tRECV -> from=%02d id=%02d group=%-15s size=%-8d Bytes n=%-14s ret=%d,%d bw=%7.4f MB/s (%s)\n",
    	  visit_srv_rank(vscd),id,"info",size,
    	  printcoord(req->n[0],req->n[1],req->n[2],req->n[3],helpstr1),ret1,ret2,size/(delta*1024*1024),
    	  visit_vtype2name(req->vtype));
       lvisit_timing_end_time("receive",local_t1);
    
       /* send data to visualization */
       lvisit_nbody_info_send(&data_info);
     }
    
     if(id==4) { /* particles */
       esize_particles = size / visit_sizeof(req->vtype);
       if(data_particles_size < esize_particles) {
          if(data_particles_size > 0) { FREEN(data_particles,"data_particles");}
          data_particles=MALLOCN(esize_particles*sizeof(double),"data_particles");
          data_particles_size=esize_particles;
       }
       
       lvisit_timing_start("read_data");
       if(! (ret1=visit_srv_read_data(vscd, data_particles, req->vtype, req)) ) goto VERROR;
       lvisit_timing_end_time_abs("read_data",vtools_dtime(0.0)-0.000003);
       visit_srv_timings(vscd, &local_t0, &local_t1, &remote_t0);
       lvisit_timing_start("send_ack2");
       if(! (ret2=visit_srv_ack2(vscd))) goto VERROR;
       lvisit_timing_end("send_ack2");
    
       delta=local_t1-local_t0;if (delta==0) delta=-1.0;
       
       lvisit_timing_start_time("receive",local_t0);
       printf("\t\tRECV -> from=%02d id=%02d group=%-15s size=%-8d Bytes n=%-14s ret=%d,%d bw=%7.4f MB/s (%s)\n",
    	  visit_srv_rank(vscd),id,"particles",size,
    	  printcoord(req->n[0],req->n[1],req->n[2],req->n[3],helpstr1),ret1,ret2,size/(delta*1024*1024),
    	  visit_vtype2name(req->vtype));
       lvisit_timing_end_time("receive",local_t1);
    
       /* send data to visualization */
       lvisit_nbody_particles_send_fld(data_particles
    				,req->n[0]
    ,req->n[1]
    
    				);
    
     }
    
    

  }
    

 VERROR:
  visit_srv_disconnect(vscd);
  lvisit_nbody_close();
  printf("----- ending  nbodysim_proxy, return code rc=0\n");

  if(loop) goto LOOP;
  return(0);
}
/* --------------------------------------------------------------- */

void scan_program_options(int argc,char **argv) {
  int count;
    for(count=0;count<argc;count++) {
      if(argv[count][0]=='-') {
          switch(argv[count][1]) {
              case 'i':
                  parm_iter=atoi(argv[++count]);
                  break;
              case 's':
                  parm_startnr=atoi(argv[++count]);
                  break;
              case 'l':
                  loop=1;
                  break;
              case 'e':
                  parm_endnr=atoi(argv[++count]);
                  break;
              case 'f':
                  infilename=argv[count+1];
                  break;
              case '?':
                  usage();
                  exit(0);
                  break;
          }
      }
  }

  printf("----- PARAMETERS --------------------------------------\n");
  printf("\t Number of Iterations: %d\n",parm_iter);
  printf("\t Start Nr:             %d\n",parm_startnr);
  printf("\t End Nr:               %d\n",parm_endnr);
  printf("-------------------------------------------------------\n");
  
}

void usage() {
    fprintf(stderr,"Usage: nbodysim                                 \n");
    fprintf(stderr,"    -i <iter>                                  \n");
    fprintf(stderr,"    -s <startnumber>                           \n");
    fprintf(stderr,"    -e <endnumber>                             \n");
    fprintf(stderr,"    -?  : print this usage help                \n");
}
