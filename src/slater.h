#ifndef SLATER_H
#define SLATER_H
#include <armadillo>

using namespace arma;


class Slater 
{
    private:
        mat m_slater_grad;
        double m_slater_lap;

        mat m_r;
        double m_alpha, m_beta, m_omega;
        int m_dimension, m_number_particles;
    public:
        // constructors
        Slater(mat r, double alpha, double beta, int dimension, \
        int number_particles, double omega);
        
        // calculates derivatives
        void UpdateDerivatives();
        mat JastrowFirstDerivative();
        
        // getters, setters
        mat GetGradient();
        double GetLaplacianSum();
}

#endif // SLATER_H
