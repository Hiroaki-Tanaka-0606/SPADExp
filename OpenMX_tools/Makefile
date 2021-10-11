option = -lm
library = inputTools.o

all: preproc.o

preproc.o: preproc_main.o $(library)
	$(CXX) preproc_main.o $(library) -o preproc.o

.cpp.o:
	$(CXX) -c $<

preproc.o: inputTools.hpp
