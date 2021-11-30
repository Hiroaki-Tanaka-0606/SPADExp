MKLROOT = /opt/intel/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64_lin
option = -L${MKLROOT} -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -liomp5 -lpthread -lm -ldl -lgfortran -fopenmp

target = calcPSF_main.o setup.o load_input.o parse_inputs.o log.o validation_Thomas_Fermi.o calculation_Thomas_Fermi.o validation_atomic_wfn.o calculation_atomic_wfn.o validation_psf.o calculation_psf.o

all: calcPSF.o

calcPSF.o: $(target)
	$(CXX) $(target) $(option) -o calcPSF.o

.cpp.o:
	$(CXX) -fopenmp -c $<


calcPSF_main.o: variables.hpp setup.hpp load_input.hpp validation_Thomas_Fermi.hpp calculation_Thomas_Fermi.hpp validation_atomic_wfn.hpp log.hpp calculation_atomic_wfn.hpp validation_psf.hpp calculation_psf.hpp
load_input.o: parse_inputs.hpp variables_ext.hpp
log.o: variables_ext.hpp
validation_Thomas_Fermi.o: log.hpp variables_ext.hpp setup.hpp
calculation_Thomas_Fermi.o: log.hpp variables_ext.hpp
validation_atomic_wfn.o: log.hpp variables_ext.hpp setup.hpp
calculation_atomic_wfn.o: log.hpp variables_ext.hpp setup.hpp
setup.o: variables_ext.hpp
