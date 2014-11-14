#ifndef SINGLEPARTICLE_H
#define SINGLEPARTICLE_H
#include <armadillo>

class SingleParticle
{
public:
    arma::mat m_r;
    double m_alpha;
    int m_dimension, m_number_particles;
    SingleParticle(arma::mat r, double alpha, int dimension, \
                   int number_particles);

    double Wavefunction();
};

#endif // SINGLEPARTICLE_H
