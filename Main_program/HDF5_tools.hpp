// calcPSF main program
// HDF5_tools
// write functions are the same as OpenMX_toole/HDF5_tools

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
void w_data_1d(Group g, const char* key, int size, double* value);
void w_data_2d(Group g, const char* key, int size1, int size2, double** value);
void w_data_3d(Group g, const char* key, int size1, int size2,
							 int size3, double*** value);
void w_data_4d(Group g, const char* key,
							 int size1, int size2, int size3, int size4, double**** value);

// read size and load data
void s_data_2d(Group g, const char* key, int* size1, int* size2);
void r_data_2d(Group g, const char* key, int size1, int size2, double** value);

void s_data_4d(Group g, const char* key, int* size1, int* size2, int* size3, int* size4);
void r_data_4d(Group g, const char* key, int size1, int size2, int size3, int size4,  double**** value);

void s_data_1c(Group g, const char* key, int* size, int* length);
void r_data_1c(Group g, const char* key, int size, int length, char** value);

string r_att_str(Group g, const char* key);
int r_att_int(Group g, const char* key);
bool r_att_bool(Group g, const char* key);
void r_att_1d(Group g, const char* key, int size, double* value);
double r_att_double(Group g, const char* key);
