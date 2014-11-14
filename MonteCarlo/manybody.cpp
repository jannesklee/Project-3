#include "manybody.h"
#include <armadillo>

/* -------------------------------------------------------------------------- *
 *        Function to compute the squared wave function, simplest form        *
 * -------------------------------------------------------------------------- */
double ManyBody::UnperturbedWavefunction()
{
  int i, j, k;
  double wf, argument, r_single_particle, r_12;

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
  double a, C;
  double wf, argument, omega, r_single_particle, r_12;
  //TODO: implement a and C in a proper way
  a = 1.0; // antiparallel spin
  C = 1.0; // normalization
  omega = 1.0;

  argument = wf = 0;
  for (i = 0; i < m_number_particles; i++) {
    r_single_particle = 0;
    for (j = 0; j < m_dimension; j++) {
      r_single_particle  += m_r(i,j)*m_r(i,j);
    }
    argument += r_single_particle; //TODO: check if removed squareroot is reasonable
  }

  // TODO: copied from below
  for (i = 0; i < m_number_particles-1; i++) { // for 2 electrons no loop
    for (j = i+1; j < m_number_particles; j++) {
      r_12 = 0;
      for (k = 0; k < m_dimension; k++) {
        r_12 += (m_r(i,k)-m_r(j,k))*(m_r(i,k)-m_r(j,k));
      }
    }
  }

  wf = C*exp(-m_alpha*omega*argument*0.5)*exp(a*r_12/(1.+ m_beta*r_12));
  return wf;
}


ManyBody::ManyBody(arma::mat r, double alpha, double beta, int dimension, int number_particles)
{
    m_r = r;
    m_alpha = alpha;
    m_beta = beta;
    m_dimension = dimension;
    m_number_particles = number_particles;
}
