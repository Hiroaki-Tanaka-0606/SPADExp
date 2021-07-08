option = -lm
target = calcPSF_main.o setup.o load_input.o parse_inputs.o log.o validation_Thomas_Fermi.o calculation_Thomas_Fermi.o

all: calcPSF.o

calcPSF.o: $(target)
	$(CXX) $(target) $(option) -o calcPSF.o

.cpp.o:
	$(CXX) -c $<

calcPSF_main.o: variables.hpp setup.hpp load_input.hpp validation_Thomas_Fermi.hpp calculation_Thomas_Fermi.hpp log.hpp
load_input.o: parse_inputs.hpp variables_ext.hpp
log.o: variables_ext.hpp
validation_Thomas_Fermi.o: log.hpp variables_ext.hpp setup.hpp
calculation_Thomas_Fermi.o: log.hpp variables_ext.hpp
setup.o: variables_ext.hpp
