#include <iomanip>
#include <armadillo>
#include "singleparticle.h"
#include "manybody.h"
#include <omp.h>

using namespace arma;
using namespace std; 

double ManyBody::UnperturbedWavefunction()
{
  int i, j;
  double wf, argument, r_single_particle;

  argument = wf = 0;
  for (i = 0; i < m_number_particles; i++) {
    r_single_particle = 0;
    for (j = 0; j < m_dimension; j++) {
      r_single_particle  += m_r(i,j)*m_r(i,j);
    }
    argument += sqrt(r_single_particle);
  }
  wf = exp(-argument*m_alpha) ;
  return wf;
}

double ManyBody::PerturbedWavefunction()
{
  int i, j, k;
  double wf, argument, r_single_particle, r_12;
  double a = 1.;

  argument = wf = 0;
  for (i = 0; i < m_number_particles; i++) {
    r_single_particle = 0;
    for (j = 0; j < m_dimension; j++) {
      r_single_particle  += m_r(i,j)*m_r(i,j);
    }
    argument += r_single_particle; 
  }

  // calculates the relative distance // TODO: At that point too much calculation (loop)
  r_12 = 0;
  for (i = 0; i < m_number_particles-1; i++) { 
    for (j = i+1; j < m_number_particles; j++) {
      r_12 = 0;
      for (k = 0; k < m_dimension; k++) {
        r_12 += (m_r(i,k)-m_r(j,k))*(m_r(i,k)-m_r(j,k));
      }
    }
  }
  r_12 = sqrt(r_12);

  wf = exp(-m_alpha*m_omega*argument*0.5)*exp(a*r_12/(1.+ m_beta*r_12));
  return wf;
}

double ManyBody::SixElectronSystem()
{
    int i, j, k, n2;
    double a, wf, r_12, psi_c;
    double det_slater_up, det_slater_down;
    mat slater_up, slater_down;
    vec m_r1, m_r2, m_r3, m_r4, m_r5, m_r6;
    SingleParticle particle[6];

    n2 = m_number_particles/2;  
    slater_up = mat(n2, n2);  
    slater_down = mat(n2, n2);

    // ----------------------- Pade-Jastrow factor -------------------------- //
    wf = 0;
    // calculates the relative distance 
    psi_c = 0.0;
    for (i = 0; i < m_number_particles-1; i++) { 
        for (j = i+1; j < m_number_particles; j++) {
            r_12 = 0;
            for (k = 0; k < m_dimension; k++) {
                r_12 += (m_r(i,k)-m_r(j,k))*(m_r(i,k)-m_r(j,k));
            }
            r_12 = sqrt(r_12);
                // evaluate a for parallel or antiparallel spin
                if ((i < n2 && j < n2) || (i >= n2 && j >= n2)) { 
                    a = 1./3.;
                }
                else {
                    a = 1.;
                }
            psi_c += a*r_12/(1. + m_beta*r_12);
        }
    }
    psi_c = exp(psi_c); 

    // ------------------------ Unpurturbed part ---------------------------- //
    // spin-up particles
    particle[0].SetAll(conv_to<vec>::from(m_r.row(0)), 0, 0, m_dimension,\
            m_omega, m_alpha);
    particle[1].SetAll(conv_to<vec>::from(m_r.row(0)), 1, 0, m_dimension,\
            m_omega, m_alpha);
    particle[2].SetAll(conv_to<vec>::from(m_r.row(0)), 0, 1, m_dimension,\
            m_omega, m_alpha);

    // filling of slater matrix
    for (i = 0; i < n2; i++) {
        for (j = 0; j < n2; j++) {
            // spin up
            particle[i].SetPosition(conv_to<vec>::from(m_r.row(j))); 
            slater_up(i,j) = particle[i].Wavefunction();
            // spin down
            particle[i].SetPosition(conv_to<vec>::from(m_r.row(j+n2))); 
            slater_down(i,j) = particle[i].Wavefunction();
        }
    }

    det_slater_up = det(slater_up);
    det_slater_down = det(slater_down);

    // ---------------------- calculate wavefunction ------------------------ //
    
    wf = det_slater_up*det_slater_down*psi_c;
    return wf;
}

// ------------------------ setters, getters -------------------------------- //
void ManyBody::SetPosition(mat r) {
    m_r = r;
}

void ManyBody::SetVariables(double alpha, double beta) {
    m_alpha = alpha;
    m_beta = beta;
}

// ------------------------- constructors ----------------------------------- //
ManyBody::ManyBody()
{
}

ManyBody::ManyBody(int dimension, int number_particles, double omega)
{
    m_dimension = dimension;
    m_number_particles = number_particles;
    m_omega = omega;
}

ManyBody::ManyBody(double alpha, double beta, int dimension, \
        int number_particles, double omega)
{
    m_alpha = alpha;
    m_beta = beta;
    m_dimension = dimension;
    m_number_particles = number_particles;
    m_omega = omega;
}

ManyBody::ManyBody(mat r, double alpha, double beta, int dimension, \
        int number_particles, double omega)
{
    m_r = r;
    m_alpha = alpha;
    m_beta = beta;
    m_dimension = dimension;
    m_number_particles = number_particles;
    m_omega = omega;
}
