#ifndef SINGLEPARTICLE_H
#define SINGLEPARTICLE_H
#include <armadillo>

using namespace arma;

class SingleParticle
{
    protected:
        vec m_r;
        int m_nx, m_ny, m_dimension;
        double m_omega;
    public:
        SingleParticle(vec r, int nx, int ny, int dimension, double omega);
        SingleParticle();
        void SetPosition(vec r);
        void SetAll(vec r, int nx, int ny, int dimension, double omega);
        double Wavefunction();
};

#endif // SINGLEPARTICLE_H
