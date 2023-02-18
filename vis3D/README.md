# 3D visualization

[Go to top](../README.md)  

## How to run

### compile 
This is the instruction for CfCA hydro school 2022. First login the analysis server.

    ssh <your account>@an**.cfca.nao.ac.jp
    
Then copy the source code if you have not copy it yet.

    cd /cfca-work/<your account>
    cp -r /cfca-work/hydro17/HydroSchool2022 .
To run the code, you need to compile `main.f90`.
    
    cd HydroSchool2022/vis3D
    module load intel
    make makedata.x
    
Then `makedata.x`is made in this directory.

