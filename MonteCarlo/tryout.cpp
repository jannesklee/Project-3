#include <fstream>
#include <iomanip>
#include <armadillo>
#include "singleparticle.h"
#include "lib.h"

using namespace arma;

int main()
{  
    int j;
    double step_length = 2.0;
    vec r_new, r_old;
    double wf;
    long idum = -1.;
    int dimension = 2;

    r_new = zeros<vec>(dimension);
    r_old = zeros<vec>(dimension);

    for (j = 0; j < dimension; j++) {
        r_old(j) = step_length*(ran1(&idum)-0.5);
    }
    //r_old(0) = 0.3;
    //r_old(1) = -0.1;

    SingleParticle particle_1(r_old, 0, 0, dimension, 1.0);
    wf = particle_1.Wavefunction();

    cout << r_old << setw(15) << wf << endl;

    //r_new(0) = -0.2;
    //r_new(1) = -0.4;

    for (j = 0; j < dimension; j++) {
        r_new(j) = r_old(j) + step_length*(ran1(&idum)-0.5);
    }
    particle_1.SetPosition(r_new);
    wf = particle_1.Wavefunction();

    cout << r_new << setw(15) << wf << endl;

    return 0;
}
