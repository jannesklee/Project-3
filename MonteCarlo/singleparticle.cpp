#include "singleparticle.h"
#include <armadillo>

double SingleParticle::Wavefunction()
{
  int i, j;
  double C;
  double wf, argument, r_single_particle, hermite;

  wf = 0;
  r_single_particle = 0;
  for (i = 0; i < m_dimension; i++) {  // TODO: This has to be changed. Not reasonable
    for (j = 0; j < m_dimension; j++) {
      r_single_particle  += m_r(i,j)*m_r(i,j);
    }
  }
  argument = sqrt(r_single_particle);

  // Hermite polynomials and normalization factor C for n = 0,1,2
  hermite = 0;
  C = 1;
  if (m_n==0){
      hermite = 1;
      C = 1*sqrt(sqrt(m_omega/M_PI));
  }
  if (m_n==1){
      hermite = 2*sqrt(m_omega)*r_single_particle;
      C = 0.5*sqrt(sqrt(m_omega/M_PI));
  }
  if (m_n==2){
      hermite = 4*m_omega*r_single_particle*r_single_particle-2;
      C = 1/8*sqrt(sqrt(m_omega/M_PI));
  }


  //single particle wave function
  wf = C*exp(-m_omega*argument*0.5)*hermite;
  return wf;
}


SingleParticle::SingleParticle(arma::mat r,int nx, int dimension, \
        int number_particles, double omega)
{
    m_r = r;
    m_n = nx;
    m_dimension = dimension;
    m_number_particles = number_particles;
    m_omega = omega;
}
