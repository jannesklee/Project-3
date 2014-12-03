#ifndef SINGLEPARTICLE_H
#define SINGLEPARTICLE_H
#include <armadillo>

using namespace arma;

//! Single Particle class
class SingleParticle
{
    protected:
        vec m_r;
        int m_nx, m_ny, m_dimension;
        double m_omega, m_alpha;
    public:
        SingleParticle(vec r, int nx, int ny, int dimension, double omega,\
                double alpha);
        SingleParticle();
        void SetPosition(vec r);
        void SetAll(vec r, int nx, int ny, int dimension, double omega,\
                double alpha);
        double Wavefunction();
        vec GetGradient();
        double GetLaplacian();
};

#endif // SINGLEPARTICLE_H
