#include <armadillo>
#include <iomanip>
#include <math.h>
#include "singleparticle.h"
#include "manybody.h"
#include "jastrow.h"
#include "slater.h"

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
void Slater::SetPosition(mat r) {
    m_r = r; 
}

void Slater::SetupSixElectron(){
    int i, j;

    SingleParticle particle[6];

    m_slater_up = mat(m_number_particles/2, m_number_particles/2, fill::zeros);
    m_slater_down = mat(m_number_particles/2, m_number_particles/2, fill::zeros);

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
    for (i = 0; i < m_number_particles/2; i++) {
        for (j = 0; j < m_number_particles/2; j++) {
            // spin up
            particle[i*2].SetPosition(conv_to<vec>::from(m_r.row(j))); 
            m_slater_up(i,j) = particle[i*2].Wavefunction();
            // spin down
            particle[i*2+1].SetPosition(conv_to<vec>::from(m_r.row(j+3))); 
            m_slater_down(i,j) = particle[i*2+1].Wavefunction();
        }
    }
}

vec Slater::Gradient(int i) {
    int k, j, offset; 
    vec slater_grad, particle_grad;
    mat m_slater, m_slater_inv;
    SingleParticle particle[m_number_particles/2];

    slater_grad = vec(m_dimension);
    m_slater = mat(m_number_particles/2, m_number_particles/2); 
    m_slater_inv = mat(m_number_particles/2, m_number_particles/2); 


    if (i < m_number_particles/2) {
        m_slater = m_slater_up;
        offset = 0;
    }
    else {
        m_slater = m_slater_down;
        offset = m_number_particles/2; 
    }

    particle[0].SetAll(conv_to<vec>::from(m_r.row(0)), 0, 0, m_dimension,\
            m_omega, m_alpha);
    particle[1].SetAll(conv_to<vec>::from(m_r.row(0)), 1, 0, m_dimension,\
            m_omega, m_alpha);
    particle[2].SetAll(conv_to<vec>::from(m_r.row(0)), 0, 1, m_dimension,\
            m_omega, m_alpha);

    m_slater_inv = inv(m_slater);

    for (k = 0; k < m_number_particles/2; k++) {
        particle[k].SetPosition(conv_to<vec>::from(m_r.row(i)));
        particle_grad = particle[k].GetGradient();

        for (j = 0; j < m_dimension; j++) {
            slater_grad(j) += particle_grad(j)*m_slater_inv(k,i-offset);    
        }
    }    

    return slater_grad; 
}

double Slater::Laplacian(int i) {
    int k, offset; 
    double slater_lap, particle_lap;
    mat m_slater, m_slater_inv;
    SingleParticle particle[m_number_particles/2];

    m_slater = mat(m_number_particles/2, m_number_particles/2); 
    m_slater_inv = mat(m_number_particles/2, m_number_particles/2); 

    if (i < m_number_particles/2) {
        m_slater = m_slater_up;
        offset = 0;
    }
    else {
        m_slater = m_slater_down;
        offset = 3; 
    }

    particle[0].SetAll(conv_to<vec>::from(m_r.row(0)), 0, 0, m_dimension,\
            m_omega, m_alpha);
    particle[1].SetAll(conv_to<vec>::from(m_r.row(0)), 1, 0, m_dimension,\
            m_omega, m_alpha);
    particle[2].SetAll(conv_to<vec>::from(m_r.row(0)), 0, 1, m_dimension,\
            m_omega, m_alpha);

    m_slater_inv = inv(m_slater);

    slater_lap = 0;
    for (k = 0; k < m_number_particles/2; k++) {
        particle[k].SetPosition(conv_to<vec>::from(m_r.row(i)));
        particle_lap = particle[k].GetLaplacian();
        slater_lap += particle_lap*m_slater_inv(k,i-offset);
    }    

    return slater_lap;
}

