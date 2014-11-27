#ifndef MANYBODY_H
#define MANYBODY_H
#include <armadillo>

using namespace arma;
//! The ManyBody class

//! Sets up a system of the multiply particles. Yet implemented are perturbed
//! and unperturbed wave function for two electron systems and a purturbed
//! version for a one electron system 
class ManyBody
{
    protected:
        mat m_r;
        double m_alpha, m_beta, m_omega;
        int m_dimension, m_number_particles;
    public:
        ManyBody();
        ManyBody(mat r, double alpha, double beta, int dimension, \
                       int number_particles, double omega);
        ManyBody(double alpha, double beta, int dimension, \
                       int number_particles, double omega);
        ManyBody(int dimension, int number_particles, double omega);
        double PerturbedWavefunction();
        double UnperturbedWavefunction();
        double SixElectronSystem();
        void SetPosition(mat r);
        void SetVariables(double alpha, double beta);
};

#endif // MANYBODY_H
