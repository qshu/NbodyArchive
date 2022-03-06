#include <stdio.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <errno.h>

#define BUFFERED_IO 1
#if BUFFERED_IO
FILE * out_dump_file;


write_array_(array, length)
    int * array;
    int * length;
{
    if(fwrite((char *) array,sizeof(*array),*length,out_dump_file) != *length){
	printf("(write_array) failed (length %d) \n",*length);
	return 1;
    }
/*    sleep(1);*/
#ifdef HARP3TESTRUN
    harp3_testrun_();
#endif    
    fflush(out_dump_file);
    return 0;
}

write_array_by_putc(array, length)
    int * array;
    int * length;
{
    unsigned char * p = (unsigned char*) array;
    int i, nbytes;
    nbytes = sizeof(int)*(*length);
    for(i=0;i<nbytes;i++){
	fputc(*p, out_dump_file);
	p++;
    }
    return 0;
}

read_array_(array, length)
    int * array;
    int * length;
{
    if(fread((char *) array,sizeof(*array),*length,out_dump_file) != *length){
	printf("(read_array)  failed (length %d) \n",*length);
	return 1;
    }
    return 0;
}


open_dump_(unit_number, mode)
    int * unit_number;
    int * mode; /* mode =0 : read, mode = 1 : write */
{
    char fname[100];
    sprintf(fname, "NBODY_DUMP%d", *unit_number);
    fprintf(stderr,"(open_dump) name =%s mode = %d\n", fname, *mode);
    if(*mode ==0){
	out_dump_file = fopen(fname, "r");
    }else{
#ifdef HARP3TESTRUN
	harp3_testrun_();
#endif	
	unlink(fname);
	out_dump_file = fopen(fname, "w+");
#ifdef HARP3TESTRUN
        harp3_testrun_();
#endif	
    }
    if (out_dump_file == NULL){
	printf("Open failed for file %s mode %d\n", fname, mode);
    }
}

close_dump_()
{
  fprintf(stderr,"(close_dump) called\n:");
    if(fclose(out_dump_file)){
	printf("Close failed for dump file\n");
    }
}

#else /* not BUFFERED_IO */

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

int out_dump_file;


write_array_(array, length)
    int * array;
    int * length;
{
    if(write(out_dump_file, (char *) array,sizeof(*array)*(*length))
	     != (*length) * sizeof(*array)){
	printf("(write_array) failed (length %d) \n",*length);
	return 1;
    }
    return 0;
}


read_array_(array, length)
    int * array;
    int * length;
{
    if(read(out_dump_file, (char *) array,sizeof(*array)*(*length))
	     != (*length) * sizeof(*array)){
	printf("(write_array) failed (length %d) \n",*length);
	return 1;
    }
    return 0;
}


open_dump_(unit_number, mode)
    int * unit_number;
    int * mode; /* mode =0 : read, mode = 1 : write */
{
    char fname[100];
    sprintf(fname, "NBODY_DUMP%d", *unit_number);
    fprintf(stderr,"(open_dump) name =%s mode = %d\n", fname, *mode);
    if(*mode ==0){
	out_dump_file = open(fname, O_RDONLY);
    }else{
	unlink(fname);
	out_dump_file = open(fname, O_WRONLY | O_CREAT|O_TRUNC,
			     S_IRUSR|S_IWUSR);
    }
    if (out_dump_file == -1){
	printf("Open failed for file %s mode %d\n", fname, mode);
    }
}

close_dump_()
{
    if(close(out_dump_file)){
	printf("Close failed for dump file\n");
    }
}

#endif
#ifdef MEMORY_LOCKING


void lock_memory_()
{
    int retval;
    if(getuid() == 0){
	/* ROOT - can mlock */
	retval = mlockall (MCL_FUTURE);
	if(retval){
	    int err = errno;
	    printf("lock_memory:failed err = %d\n", err);
	}
    }else{
	fprintf(stderr,"lock_memory:not root\n");
    }
}


void lock_memory_future_()
{
    int retval;
    if(getuid() == 0){
	/* ROOT - can mlock */
	retval = mlockall (MCL_FUTURE);
	if(retval){
	    int err = errno;
	    printf("lock_memory:failed err = %d\n", err);
	}
    }else{
	fprintf(stderr,"lock_memory_future:not root\n");
    }
}

void lock_memory_current_()
{
    int retval;
    if(getuid() == 0){
	/* ROOT - can mlock */
	retval = mlockall (MCL_CURRENT);
	if(retval){
	    int err = errno;
	    printf("lock_memory_current:failed err = %d\n", err);
	}
    }else{
	fprintf(stderr,"lock_memory:not root\n");
    }
}


void unlock_memory_()
{
    int retval;
    if(getuid() == 0){
	/* ROOT - can mlock */
	retval = munlockall();
	if(retval){
	    int err = errno;
	    printf("unlock_memory:failed err = %d\n", err);
	}
    }else{
	fprintf(stderr,"unlock_memory:not root\n");
    }
	
}
#else

void lock_memory_()
{
}

void unlock_memory_()
{
}
#endif

    
    
