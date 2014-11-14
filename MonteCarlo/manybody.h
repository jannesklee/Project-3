#ifndef MANYBODY_H
#define MANYBODY_H
#include <armadillo>

class ManyBody
{
public:
    arma::mat m_r;
    double m_alpha, m_beta;
    int m_dimension, m_number_particles;
    ManyBody(arma::mat r, double alpha, double beta, int dimension, \
                   int number_particles);
    double Wavefunction();
};

#endif // MANYBODY_H
