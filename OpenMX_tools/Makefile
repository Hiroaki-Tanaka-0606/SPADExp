option = -lm
library = inputTools.o
h5cpp = ~/hdf5-1.12.1/hdf5/bin/h5c++

all: preproc.o postproc.o

preproc.o: preproc_main.o $(library)
	$(CXX) preproc_main.o $(library) -o preproc.o

postproc.o: postproc_main.o HDF5_tools.o $(library)
	$(h5cpp) postproc_main.o HDF5_tools.o $(library) -o postproc.o

.cpp.o:
	$(CXX) -c $<

postproc_main.o: postproc_main.cpp
	$(h5cpp) -c postproc_main.cpp

HDF5_tools.o: HDF5_tools.cpp
	$(h5cpp) -c HDF5_tools.cpp

preproc.o: inputTools.hpp
postproc.o: inputTools.hpp HDF5_tools.hpp
