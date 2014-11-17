#ifndef SINGLEPARTICLE_H
#define SINGLEPARTICLE_H
#include <armadillo>

class SingleParticle
{
public:
    arma::mat m_r;
    int m_n, m_dimension, m_number_particles;
    double m_omega;
    SingleParticle(arma::mat r, int nx, int dimension, \
                   int number_particles, double omega);
    double Wavefunction();
};

#endif // SINGLEPARTICLE_H
