# directory specified by --prefix=DIR in ./configure of HDF5
# written in config.log
# hdf5-x.y.z/hdf5 by default
hdf5 = /home/hiroaki/hdf5-1.12.1/hdf5

# compile option can be checked by
# >$(hdf5)/bin/h5c++ (filename) -show
# >$(hdf5)/bin/h5c++ (filename) -c -show

# common option
commonOption = -DOLD_HEADER_FILENAME -DHDF_NO_NAMESPACE -DNO_STATIC_CAST 

# compile option: $(hdf5)/include contains the header file "H5Cpp.h"
compileOption = $(commonOption) -I $(hdf5)/include

# link option 1: -L addes the path where libraries are searched
linkOption1 = -L $(hdf5)/lib

# link option 2: necessary libraries in HDF5
# -l XX addes the library libXX.a (static) or libXX.so (shared)
# -Wl,-rpath,$(hdf5)/lib tells "ld" command to search $(hdf5)/lib directory 
linkOption2 = -l hdf5_hl_cpp -l hdf5_cpp -l hdf5_hl -l hdf5 -Wl,-rpath,$(hdf5)/lib

# configuration for intel MKL
MKLROOT = /opt/intel/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64_lin
option = -L${MKLROOT} -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -liomp5 -lpthread -lm -ldl -lgfortran -fopenmp $(linkOption1) $(linkOption2)

target = calcPSF_main.o setup.o load_input.o parse_inputs.o log.o validation_Thomas_Fermi.o calculation_Thomas_Fermi.o validation_atomic_wfn.o calculation_atomic_wfn.o validation_scf_atom.o calculation_scf_atom.o validation_phase_shift.o calculation_phase_shift.o validation_psf.o calculation_psf.o HDF5_tools.o physical_tools.o

all: calcPSF.o

calcPSF.o: $(target)
	$(CXX) $(target) $(option) -o calcPSF.o

.cpp.o:
	$(CXX) -fopenmp $(compileOption) -c $<


calcPSF_main.o: variables.hpp setup.hpp load_input.hpp validation_Thomas_Fermi.hpp calculation_Thomas_Fermi.hpp validation_atomic_wfn.hpp log.hpp calculation_atomic_wfn.hpp validation_scf_atom.hpp calculation_scf_atom.hpp validation_phase_shift.hpp calculation_phase_shift.hpp validation_psf.hpp calculation_psf.hpp HDF5_tools.hpp
load_input.o: parse_inputs.hpp variables_ext.hpp
log.o: variables_ext.hpp
validation_Thomas_Fermi.o: log.hpp variables_ext.hpp setup.hpp
calculation_Thomas_Fermi.o: log.hpp variables_ext.hpp
validation_atomic_wfn.o: log.hpp variables_ext.hpp setup.hpp
calculation_atomic_wfn.o: log.hpp variables_ext.hpp setup.hpp
setup.o: variables_ext.hpp
validation_psf.o: log.hpp variables_ext.hpp
calculation_psf.o: log.hpp variables_ext.hpp HDF5_tools.hpp physical_tools.hpp
validation_scf_atom.o: log.hpp variables_ext.hpp setup.hpp
calculation_scf_atom.o: log.hpp variables_ext.hpp setup.hpp calculation_atomic_wfn.hpp HDF5_tools.hpp
validation_phase_shift.o: log.hpp variables_ext.hpp setup.hpp
calculation_phase_shift.o: log.hpp variables_ext.hpp setup.hpp calculation_atomic_wfn.hpp HDF5_tools.hpp
HDF5_tools.o: log.hpp
