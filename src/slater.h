#ifndef SLATER_H
#define SLATER_H
#include <armadillo>

using namespace arma;

class Slater
{
    private:
        mat m_slater_up;
        mat m_slater_down;
        mat m_slater_invup;
        mat m_slater_invdown;

        mat m_r;
        double m_alpha, m_beta, m_omega;
        int m_dimension, m_number_particles;
    public:
        // constructors
        Slater();
        Slater(mat r, double alpha, double beta, int dimension, \
        int number_particles, double omega);

        void SetupSixElectron();
        void SetPosition(mat r); 
        vec Gradient(int i);
        double Laplacian(int i);
        void InitInverseSlaterMatrices();
        void UpdateInverseSlaterMatrices();
};

#endif // SLATER_H
