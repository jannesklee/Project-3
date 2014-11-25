#include <fstream>
#include <iomanip>
#include <armadillo>
#include "singleparticle.h"
#include "lib.h"
#include "manybody.h"

using namespace arma;

int main()
{  
    int j, i;
    double step_length = 2.0;
    mat r_new, r_old;
    double wf1, wf2;
    long idum = -1.;
    int dimension = 2;
    int number_particles = 6;
    double alpha = 0.6;
    double beta = 0.1;
    double omega = 1.0;

    r_new = zeros<mat>(number_particles, dimension);
    r_old = zeros<mat>(number_particles, dimension);

    for (i = 0; i < number_particles; i++) {
        for (j = 0; j < dimension; j++) {
            r_old(i,j) = step_length*(ran1(&idum)-0.5);
        }
    }
    //r_old(0) = 0.3;
    //r_old(1) = -0.1;

    //SingleParticle particle_1(r_old, 0, 0, dimension, 1.0);
    //wf1 = particle_1.Wavefunction();

    cout << "Start Test:" << endl;

    ManyBody System(r_old, alpha, beta, dimension, number_particles, omega);
    wf1 = System.SixElectronSystem();
    //wf1 = 1.9;

    //SingleParticle particle[4];
    //cout << r_old << setw(15) << wf1 << setw(15) << wf2 << endl;

    //r_new(0) = -0.2;
    //r_new(1) = -0.4;
    //
    for (i = 0; i < number_particles; i++) {
        for (j = 0; j < dimension; j++) {
            r_new(i,j) = r_old(i,j) + step_length*(ran1(&idum)-0.5);
        }
    }

    //particle_1.SetPosition(r_new);
    //wf1 = particle_1.Wavefunction();

    //cout << r_new << setw(15) << wf1 << endl;

    return 0;
}
