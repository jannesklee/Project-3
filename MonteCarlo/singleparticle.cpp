#include "singleparticle.h"
#include <armadillo>

double SingleParticle::Wavefunction()
// Function to compute the squared wave function, simplest form
//double  wave_function(mat r, double alpha,int dimension, int number_particles) //this function is limited to two electrons
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


SingleParticle::SingleParticle(arma::mat r, double alpha,int dimension, int number_particles)
{
    m_r = r;
    m_alpha = alpha;
    m_dimension = dimension;
    m_number_particles = number_particles;
}
