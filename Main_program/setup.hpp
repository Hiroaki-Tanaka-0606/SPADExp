#include <H5Cpp.h>
using namespace std;
using namespace H5;

void initialize();
void finalize();
void initialize_radial_grid();
void generate_radial_grid();
void setup_radial_grid();
int setup_potential(int Z, double mu);
void modify_potential(int Z, double mu);
double load_potential_H5(Group atomG, double mu);
