# 3D visualization

[Go to top](../README.md)  

This is the instruction for CfCA hydro school 2022. 

## Make HDF files and xdmf text

To learn about VisIt and paraview, you need to make the data for it. First login the analysis server.

    ssh <your account>@an**.cfca.nao.ac.jp
    
Then copy the source code if you have not copy it yet.

    cd /cfca-work/<your account>
    cp -r /cfca-work/hydro17/HydroSchool2022 .
To run the code, you need to compile `main.f90`.
    
    cd HydroSchool2022/vis3D
    module load intel
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/takiwkkz/hdf5/lib
    make makedata.x
    
Then `makedata.x`is made in this directory.
    
    ./makedata.x
    
You can find the data in `hdfdata`.

## Visualization

## Preparation
Abobe, you used a library in `~/hydro17/hdf5`. After the school you cannot use it. You need to install `hdf5` in your local dir. 
