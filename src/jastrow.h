#ifndef JASTROW_H
#define JASTROW_H
#include <armadillo>

using namespace arma;
//! The Jastro class

//! This class delivers either closed form or brute force solutions of the 
//! derivations of the Jastro-factor.
class Jastrow
{
    private:
        mat m_jastrow_grad;
        double m_jastrow_lap;

        mat m_r;
        double m_alpha, m_beta, m_omega;
        int m_dimension, m_number_particles;
    public:
        // constructors
        Jastrow(mat r, double alpha, double beta, int dimension, \
        int number_particles, double omega);

        // calculates derivatives
        vec Gradient(int);
        double Laplacian(int);
        double Factor(); 

        // getters, setters
        vec GetGradient();
        double GetLaplacian();
};

#endif // JASTROW_H
