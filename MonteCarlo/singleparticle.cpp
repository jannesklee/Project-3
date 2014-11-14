#include "singleparticle.h"
#include <armadillo>

double SingleParticle::Wavefunction()
{
  int i, j, k;
  double C;
  double wf, argument, omega, r_single_particle, hermite;
  C = 1.0;
  //TODO: implement C in a proper way
  omega = 1.0;


  wf = 0;
  r_single_particle = 0;
  for (j = 0; j < m_dimension; j++) {
      r_single_particle  += m_r(i,j)*m_r(i,j);
    }
  argument = sqrt(r_single_particle);

  // Hermite polynomials for n = 0,1,2
  hermite = 0;
  if (m_n=0){
      hermite = 1;
  }
  if (m_n=1){
      hermite = 2*sqrt(omega)*r_single_particle;
  }
  if (m_n=2){
      hermite = 4*omega*r_single_particle*r_single_particle-2;
  }

  //single particle wave function
  wf = C*exp(-omega*argument*0.5)*hermite;
  return wf;
}


SingleParticle::SingleParticle(arma::mat r,int nx, int dimension, int number_particles)
{
    m_r = r;
    m_n = nx;
    m_dimension = dimension;
    m_number_particles = number_particles;
}
