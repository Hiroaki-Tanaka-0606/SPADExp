// calcPSF tools for OpenMX
// HDF5_tools

#include <H5Cpp.h>
#include <string>

using namespace std;
using namespace H5;

// for attributes
void w_att_str(Group g, const char* key, string value);
void w_att_int(Group g, const char* key, int value);
void w_att_double(Group g, const char* key, double value);
void w_att_bool(Group g, const char* key, bool value);
void w_att_1d(Group g, const char* key, int size, double* value);
void w_att_2d(Group g, const char* key, int size1, int size2, double** value);

// for datasets
void w_data_1c(Group g, const char* key, int size, int length, char** value);
void w_data_2d(Group g, const char* key, int size1, int size2, double** value);
void w_data_4d(Group g, const char* key,
							 int size1, int size2, int size3, int size4, double**** value);
