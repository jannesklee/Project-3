#ifndef MANYBODY_H
#define MANYBODY_H
#include <armadillo>

using namespace arma;
//! \brief{Class for the calculation of multiple particles.}

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

//!\brief trial wave function for six electron system
 
//! SixElectronSystem calculates the trial wave function in two steps. The first
//! one consists of the calculation of the Jastrow-factor, the second the unpur-
//! turbed part of the wavefunction. The latter thereby fills in two slater 
//! determinants for either spin up or spin down particles and evaluates the
//! determinant applying Sarrus' rule.
        double SixElectronSystem();
        void SetPosition(mat r);
        void SetVariables(double alpha, double beta);
};

#endif // MANYBODY_H
