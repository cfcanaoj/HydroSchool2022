dird := hdfdata

exe=makedata.x

all: ${exe}

.PHONY: all clean allclean

fc=ifort
fopt=-g -traceback -O2

HDF5 = /home/takiwkkz/hdf5
FLAGHDF =  -I${HDF5}/include -L${HDF5}/lib -lhdf5 -lhdf5_fortran 

${exe}: main.f90
	${fc} ${fopt} ${FLAGHDF} $< -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -f ${dird}/*
