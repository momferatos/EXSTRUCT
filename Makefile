FC=nvfortran
ARCH=haswell
DEBUGFLAGS=-g
OPTFLAGS=-fast -fastsse -tp=$(ARCH)
HDF5ROOT=/home/giorgos/libs/hdf5/
LDFLAGS=-lhdf5_fortran -lhdf5
FCFLAGS=-cpp $(OPTFLAGS) -I $(HDF5ROOT)/include -L $(HDF5ROOT)/lib


all: exstruct.f90
	$(FC) $(FCFLAGS) -o exstruct.exe exstruct.f90 $(LDFLAGS)

clean:
	rm -rf *.o *.mod
