# 3D visualization

[Go to top](../README.md)  

This document is used for CfCA hydro school 2022. 

## Make HDF files and xdmf text

To learn about VisIt and paraview, you need to make the data for it. First login the analysis server.

    ssh <your account>@an**.cfca.nao.ac.jp
    
To treat HDF5 format, you need library. Add the follwoing command in your`~/.bashrc`. You can find sample `.bash_profile` and `.bashrc` in `/home/skel/`. (`tty erase '^H'` in `.basrhrc` may make a trouble. Then delete the line.)
    
    module load intel
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hydro17/hdf5/lib
    module load visit
    
Then copy the source code if you have not copy it yet.

    cd /cfca-work/<your account>
    cp -r /cfca-work/hydro17/HydroSchool2022 .
To run the code, you need to compile `main.f90`.
    
    cd HydroSchool2022/vis3D
    make makedata.x
    
Then `makedata.x`is made in this directory.
    
    ./makedata.x
    
You can find the data in `hdfdata`.

## details of xdmf
See the links.
- https://qiita.com/hsimyu/items/6e163428477d19429576
- https://www.visitusers.org/index.php?title=Using_XDMF_to_read_HDF5

## Visualization

To use VisIt, go to the [instruction](https://www.notion.so/VisIt-535a672b1e024c018c04e38015d2a249?pvs=4).

## Preparation
Abobe, you used a library in `~/hydro17/hdf5`. After the school you cannot use it. You need to install `hdf5` in your home directory in analysis server. 

### Download
Download `hdf5-1.10.6.tar.gz` from the following URL: https://www.hdfgroup.org/downloads/hdf5/source-code/
Here we use the old version of HDF5. Then put it in your home directory of analysis server. 

### Setup
The follwing command enable you use the HDF5. Now we assume the use of the intel compiller(`module load intel` has been already written in ``~/.basrhrc ). If you prefer to use another compiler, change the argument.

    
    tar xzvf hdf5-1.10.6.tar.gz
    cd hdf5-1.10.6
    ./configure --enable-fortran FC=ifort --prefix=/home/<your account>/hdf5
    make
    make install
     
Finally you have to add the follwoing command in your`~/.bashrc`.
     
     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/<your account>/hdf5/lib
    
