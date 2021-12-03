// calcPSF main program
// HDF5_tools
// write functions are the same as OpenMX_toole/HDF5_tools

#include <H5Cpp.h>
#include <string>
#include "log.hpp"

using namespace std;
using namespace H5;

void w_att_str(Group g, const char* key, string value){
	StrType str_var(PredType::C_S1, H5T_VARIABLE);
	str_var.setCset(H5T_CSET_UTF8);
	Attribute at=g.createAttribute(key, str_var, H5S_SCALAR);
	string data_s[1]={value};
	at.write(str_var, data_s);
}
		
void w_att_int(Group g, const char* key, int value){
	Attribute at=g.createAttribute(key, PredType::NATIVE_INT, H5S_SCALAR);
	int data_i[1]={value};
	at.write(PredType::NATIVE_INT, data_i);
}
		
void w_att_double(Group g, const char* key, double value){
	Attribute at=g.createAttribute(key, PredType::NATIVE_DOUBLE, H5S_SCALAR);
	double data_i[1]={value};
	at.write(PredType::NATIVE_DOUBLE, data_i);
}

void w_att_bool(Group g, const char* key, bool value){
	Attribute at=g.createAttribute(key, PredType::NATIVE_HBOOL, H5S_SCALAR);
	bool data_b[1]={value};
	at.write(PredType::NATIVE_HBOOL, data_b);
}

void w_att_1d(Group g, const char* key, int size, double* value){
	hsize_t dims[1]={(hsize_t)size};
	DataSpace ds(1, dims);
	Attribute at=g.createAttribute(key, PredType::NATIVE_DOUBLE, ds);
	at.write(PredType::NATIVE_DOUBLE, value);
}

void w_att_2d(Group g, const char* key, int size1, int size2, double** value){
	hsize_t dims[2]={(hsize_t)size1, (hsize_t)size2};
	DataSpace ds(2, dims);
	Attribute at=g.createAttribute(key, PredType::NATIVE_DOUBLE, ds);
	at.write(PredType::NATIVE_DOUBLE, value);
}

void w_data_1c(Group g, const char* key, int size, int length, char** value){
	StrType char_arr(PredType::C_S1, length);
	char_arr.setCset(H5T_CSET_UTF8);
	hsize_t dims[1]={(hsize_t)size};
	DataSpace dsp(1, dims);
	DataSet dse=g.createDataSet(key, char_arr, dsp);
	dse.write(value, char_arr);
}

void w_data_2d(Group g, const char* key, int size1, int size2, double** value){
	hsize_t dims[2]={(hsize_t)size1, (hsize_t)size2};
	DataSpace dsp(2, dims);
	DataSet dse=g.createDataSet(key, PredType::NATIVE_DOUBLE, dsp);
	dse.write(value, PredType::NATIVE_DOUBLE);
}

void w_data_3d(Group g, const char* key, int size1, int size2, int size3, double*** value){
	hsize_t dims[3]={(hsize_t)size1, (hsize_t)size2, (hsize_t)size3};
	DataSpace dsp(3, dims);
	DataSet dse=g.createDataSet(key, PredType::NATIVE_DOUBLE, dsp);
	dse.write(value, PredType::NATIVE_DOUBLE);
}


void w_data_4d(Group g, const char* key,
							 int size1, int size2, int size3, int size4, double**** value){
	hsize_t dims[4]={(hsize_t)size1, (hsize_t)size2, (hsize_t)size3, (hsize_t)size4};
	DataSpace dsp(4, dims);
	DataSet dse=g.createDataSet(key, PredType::NATIVE_DOUBLE, dsp);
	dse.write(value, PredType::NATIVE_DOUBLE);
}

void s_data_2d(Group g, const char* key, int* size1, int* size2){
	DataSet dse=g.openDataSet(key);
	DataSpace dsp=dse.getSpace();

	int rank=dsp.getSimpleExtentNdims();
	if(rank!=2){
		write_log((char*)"s_data_2d error: different rank");
		return;
	}
	hsize_t dims[2];
	dsp.getSimpleExtentDims(dims, NULL);
	*size1=dims[0];
	*size2=dims[1];
}

void r_data_2d(Group g, const char* key, int size1, int size2, double** value){
	DataSet dse=g.openDataSet(key);
	DataSpace dsp=dse.getSpace();

	int rank=dsp.getSimpleExtentNdims();
	if(rank!=2){
		write_log((char*)"r_data_2d error: different rank");
		return;
	}
	hsize_t dims[2];
	dsp.getSimpleExtentDims(dims, NULL);
	if(dims[0]!=size1 || dims[1]!=size2){
		write_log((char*)"r_data_2d error: size mismatch");
		return;
	}
	dse.read(value, PredType::NATIVE_DOUBLE);
}

void s_data_1c(Group g, const char* key, int* size, int* length){
	DataSet dse=g.openDataSet(key);
	DataSpace dsp=dse.getSpace();

	int rank=dsp.getSimpleExtentNdims();
	if(rank!=1){
		write_log((char*)"s_data_1c error: different rank");
		return;
	}
	hsize_t dims[1];
	dsp.getSimpleExtentDims(dims, NULL);
	*size=dims[0];

	StrType st=dse.getStrType();
	*length=st.getSize();
}

void r_data_1c(Group g, const char* key, int size, int length, char** value){
	DataSet dse=g.openDataSet(key);
	DataSpace dsp=dse.getSpace();

	int rank=dsp.getSimpleExtentNdims();
	if(rank!=1){
		write_log((char*)"r_data_1c error: different rank");
		return;
	}
	hsize_t dims[1];
	dsp.getSimpleExtentDims(dims, NULL);
	if(size!=(int)dims[0]){
		write_log((char*)"r_data_1c error: size mismatch");
		return;
	}
	StrType st=dse.getStrType();
	if(length!=(int)st.getSize()){
		write_log((char*)"r_data_1c error: length mismatch");
		return;
	}
	dse.read(value, st);
}

string r_att_str(Group g, const char* key){
	Attribute at=g.openAttribute(key);
	StrType st=at.getStrType();
	string data;
	at.read(st, data);
	return data;
}
		
int r_att_int(Group g, const char* key){
	Attribute at=g.openAttribute(key);
	int data;
	at.read(PredType::NATIVE_INT, &data);
	return data;
}
		
bool r_att_bool(Group g, const char* key){
	Attribute at=g.openAttribute(key);
	bool data;
	at.read(PredType::NATIVE_HBOOL, &data);
	return data;
}


void r_att_1d(Group g, const char* key, int size, double* value){
	Attribute at=g.openAttribute(key);
	DataSpace dsp=at.getSpace();

	int rank=dsp.getSimpleExtentNdims();
	if(rank!=1){
		write_log((char*)"r_att_1d error: different rank");
		return;
	}
	hsize_t dims[1];
	dsp.getSimpleExtentDims(dims, NULL);
	if(dims[0]!=size){
		write_log((char*)"r_att_1d error: size mismatch");
		return;
	}
	at.read(PredType::NATIVE_DOUBLE, value);
}

double r_att_double(Group g, const char* key){
	Attribute at=g.openAttribute(key);
	double data;
	at.read(PredType::NATIVE_DOUBLE, &data);
	return data;
}


void s_data_4d(Group g, const char* key, int* size1, int* size2, int* size3, int* size4){
	DataSet dse=g.openDataSet(key);
	DataSpace dsp=dse.getSpace();

	int rank=dsp.getSimpleExtentNdims();
	if(rank!=4){
		write_log((char*)"s_data_4d error: different rank");
		return;
	}
	hsize_t dims[4];
	dsp.getSimpleExtentDims(dims, NULL);
	*size1=dims[0];
	*size2=dims[1];
	*size3=dims[2];
	*size4=dims[3];
}

void r_data_4d(Group g, const char* key, int size1, int size2, int size3, int size4, double**** value){
	DataSet dse=g.openDataSet(key);
	DataSpace dsp=dse.getSpace();

	int rank=dsp.getSimpleExtentNdims();
	if(rank!=4){
		write_log((char*)"r_data_4d error: different rank");
		return;
	}
	hsize_t dims[4];
	dsp.getSimpleExtentDims(dims, NULL);
	if(dims[0]!=size1 || dims[1]!=size2 || dims[2]!=size3 || dims[3]!=size4){
		write_log((char*)"r_data_4d error: size mismatch");
		return;
	}
	dse.read(value, PredType::NATIVE_DOUBLE);
}
