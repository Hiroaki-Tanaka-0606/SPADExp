// calcPSF tools for OpenMX
// inputTools
// see also OpenMX/source/Inputtools.c

#include <fstream>
using namespace std;

void load_file(ifstream* fp);

int load_logical(const char* key, bool* value);
int load_int(const char* key, int* value);
int load_doublev(const char* key, const int count, double* values);
int load_str(const char* key, char* value, int length);

int find_str(const char* key);
char* get_line(int line_number);

void copyInput(ofstream* fp);
void copyInput_Weyl(ofstream* fp);

