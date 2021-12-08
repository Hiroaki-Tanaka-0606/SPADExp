// calcPSF calculation of atomic wavefuncion

void sequence_atomic_wfn();

double calc_atomic_wfn(double mu, int l, int nodes);
double calc_atomic_wfn2(double mu, int l, int nodes, double eigen_guess);

void Atomic_wfn_evolution(double mu, int l, double E, bool reverse);
