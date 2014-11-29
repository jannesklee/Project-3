#include "singleparticle.h"
#include "manybody.h"
#include "jastrow.h"
#include <armadillo>
#include <iomanip>
#include <math.h>

using namespace arma;
using namespace std;

Slater::Slater(mat r, double alpha, double beta, int dimension, \
        int number_particles, double omega)
{
    m_r = r;
    m_alpha = alpha;
    m_beta = beta;
    m_dimension = dimension;
    m_number_particles = number_particles;
    m_omega = omega;
}

void Slater::UpdateDerivatives()
{
    m_jastrow_grad = Jastrow::GetGradient();
    m_jastrow_lap  = Jastrow::GetLaplacianSum();
}

mat Slater::GetGradient()
{
    return m_jastrow_grad;
}

double Slater::GetLaplacianSum()
{
    return m_jastrow_lap;
}


