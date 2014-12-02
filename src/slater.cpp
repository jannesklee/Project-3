#include "singleparticle.h"
#include "manybody.h"
#include "jastrow.h"
#include <armadillo>
#include <iomanip>
#include <math.h>

using namespace arma;
using namespace std;

// ---------------------------- constructors -------------------------------- // 
Slater::Slater()
{
}

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


// ------------------------ setters, getters -------------------------------- //
double GetDetSlaterUp()
{
    return det(m_slater_up);
}

double GetDetSlaterDown()
{
    return det(m_slater_down);
}

// TODO: GO ON HERE
mat GetSlaterGrad()
{
    for (k = 0; k < m_dimension; k++) {
        grad
    }
    return grad_slater;
}

// ------------------------ functions --------------------------------------- //
void Slater::SetupTwoElectron(){
    SingleParticle particle[2];

    particle[0].SetAll(conv_to<vec>::from(m_r.row(0)), 0, 0, m_dimension,\
            m_omega, m_alpha);
    particle[1].SetAll(conv_to<vec>::from(m_r.row(1)), 0, 0, m_dimension,\
            m_omega, m_alpha);

    // spin up
    particle[0].SetPosition(conv_to<vec>::from(m_r.row(0))); 
    m_slater_up(0,0) = particle[0].Wavefunction();
    // spin down
    particle[1].SetPosition(conv_to<vec>::from(m_r.row(1))); 
    m_slater_down(0,0) = particle[1].Wavefunction();
}

void Slater::SetupSixElectron(){
    SingleParticle particle[6];

    // spin-up particles
    particle[0].SetAll(conv_to<vec>::from(m_r.row(0)), 0, 0, m_dimension,\
            m_omega, m_alpha);
    particle[2].SetAll(conv_to<vec>::from(m_r.row(0)), 1, 0, m_dimension,\
            m_omega, m_alpha);
    particle[4].SetAll(conv_to<vec>::from(m_r.row(0)), 0, 1, m_dimension,\
            m_omega, m_alpha);

    // spin-down particles
    particle[1].SetAll(conv_to<vec>::from(m_r.row(3)), 0, 0, m_dimension,\
            m_omega, m_alpha);
    particle[3].SetAll(conv_to<vec>::from(m_r.row(3)), 1, 0, m_dimension,\
            m_omega, m_alpha);
    particle[5].SetAll(conv_to<vec>::from(m_r.row(3)), 0, 1, m_dimension,\
            m_omega, m_alpha);

    // filling of slater matrix
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            // spin up
            particle[i*2].SetPosition(conv_to<vec>::from(m_r.row(j))); 
            m_slater_up(i,j) = particle[i*2].Wavefunction();
            // spin down
            particle[i*2+1].SetPosition(conv_to<vec>::from(m_r.row(j+3))); 
            m_slater_down(i,j) = particle[i*2+1].Wavefunction();
        }
    }
}

void Slater::InitInverseSlaterMatrices() {
    m_slater_invup = inv(m_slater_up);
    m_slater_invdown = inv(m_slater_down);
}

//! \todo{make this when the rest is working}
//void Slater::UpdateInverseSlaterMatrices() {
//        
//}

void Slater::UpdateSlaterGrad() {
    SingleParticle particle[m_number_particles];
    
    
    for (k = 0; k < m_number_particles; k++) {
            
    }    
}

