#
#  Simple python script to automate the building procedure
#  create by Francesco Rizzuto on May 2018
#
#
#
from glob import glob               
import fileinput
import sys
import os


# %%
# --- STEP 1: Choose working directory and initial conditionfile name 
Main_dir = 'peter'
#Out_dir  = '/p/home/jusers/rizzutp1/juwels'
Out_dir  = '/ptmp/mpa/frizzuto/New_models/Test'
#Out_dir = '/ptmp/mpa/frizzuto/New_models/2k_b0.2k/Model1/run0'
#Out_dir = '/ptmp/mpa/frizzuto/mcluster-master_mod/CollidingClusters/run0'

run_dir  = '/u/frizzuto/Nbody6ppGPU/nbody6++gpu-Maria-2020.6'

Configure = True
Make = True
Debug = False
Stack = False

# --------

os.getcwd()
os.chdir(run_dir)
os.system('ls')

#os.system('cp -r ../Main_all/Main_' + Main_dir + ' src')

# %%
#os.system('rm -f -r src/Main')
#os.system('mv src/Main_' + Main_dir + ' src/Main')


def Replace(file_in, old, new):
    for line in fileinput.input(file_in, inplace = 1):
        if(old in line):
            new_line = line.replace(old, new)
            sys.stdout.write(new_line)
        else:
            sys.stdout.write(line)


# %%

# Create The Make file with the appropriate setting
if(Configure):
    #rm = './configure --enable-simd=sse --disable-hdf5 --enable-mpi --enable-tools --enable-gpu --with-cuda --with-cuda_sdk --enable-mcmodel=large --enable-option-checking  --with-par=1m --disable-debug'  
    #rm = './configure --enable-simd=avx --enable-mpi --disable-tools --enable-gpu --enable-mcmodel=large --enable-option-checking --with-par=b512k --with-kmax=200000' 
    rm  = './configure --enable-mcmodel=large --enable-simd=sse --with-par=1m --enable-openmp'

    os.system(rm)
    Replace('build/Makefile', 'FC = gfortran', 'FC = mpif77')
    
    HDF5_new = 'HDF5_FLAGS = -D H5OUTPUT -I$(HDF5_HOME)/include/ -L$(HDF5_HOME)/lib/ -lhdf5_fortran'
    HDF5_old = 'HDF5_FLAGS = -D H5OUTPUT'
    Replace('build/Makefile', HDF5_old, HDF5_new)
    #
    inc_old = '     &          MLD=22,MLR=22,MLV=200,MCL=10,NCMAX=10,NTMAX=100)'
    inc_new = '     &          MLD=22,MLR=22,MLV=200,MCL=10,NCMAX=10,NTMAX=100)'
#    inc_new = '     &          MLD=22,MLR=600,MLV=200,MCL=10,NCMAX=10,NTMAX=100)'  
    Replace('include/params.h', inc_old, inc_new)

if(Stack):
    db_old = 'FFLAGS =  -O3 -fPIC -mcmodel=large -fopenmp -I../include'
    db_new = 'FFLAGS =  -O3 -fPIC -fmax-stack-var-size=700000 -mcmodel=large -fopenmp -I../include'

if(Debug):
    db_old = 'FFLAGS =  -O3 -fPIC -mcmodel=large -fopenmp -I../include'
    #db_new = 'FFLAGS =  -O3 -fPIC -g -fcheck=bounds -fbacktrace -Wall -ffpe-trap=invalid,zero,overflow,underflow -mcmodel=large -fopenmp -I../include'
    db_new = 'FFLAGS =  -O3 -fPIC -g -fcheck=bounds -fbacktrace -Wall  -mcmodel=large -fopenmp -I../include'
    Replace('build/Makefile', db_old, db_new)
 

if(Make):
    os.system("make clean")    
    os.system("make")
    #
    
    
    print('Moving executables in ' , Out_dir)
    Test_dir = run_dir + '/build'
    FilesNb = glob(Test_dir + '/nbody6++.*')
    if(Debug):
        os.system('mkdir ' + Out_dir + '/deb')
        os.system('cp /ptmp/mpa/frizzuto/Nbody6_File/debug.input ' + Out_dir + '/deb')
        os.system('cp /ptmp/mpa/frizzuto/Nbody6_File/debug.job ' + Out_dir + '/deb')
    #
    for f_old in FilesNb:
        f_new = f_old.replace(Test_dir + '/nbody6++', Main_dir)
        if(Debug):
            os.system('mv ' + f_old + ' ' + Out_dir + '/deb/' + f_new + '.deb')
        else:
            os.system('mv ' + f_old + ' ' + Out_dir + '/' + f_new )
    

'''
  --with-nmax=size        Set Maximum number of particles (will overwrite
                          value set in --with-par)
  --with-kmax=size        Set Maximum number of K.S. particles (will overwrite
                          value set in --with-par)
  --with-lmax=size        Set Maximum number of neighbor list particles (will
                          overwrite value set in --with-par)
  --with-mmax=size        Set Maximum number of merger (stable triple)
                          particles (will overwrite value set in --with-par)
'''

