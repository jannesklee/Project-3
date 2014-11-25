#ifndef MANYBODY_H
#define MANYBODY_H
#include <armadillo>

using namespace arma;

class ManyBody
{
    protected:
        mat m_r;
        double m_alpha, m_beta, m_omega;
        int m_dimension, m_number_particles;
    public:
        ManyBody(mat r, double alpha, double beta, int dimension, \
                       int number_particles, double omega);
        double PerturbedWavefunction();
        double UnperturbedWavefunction();
        double SixElectronSystem();
        void SetPosition(mat r);
};

#endif // MANYBODY_H
