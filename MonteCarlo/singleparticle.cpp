#include "singleparticle.h"
#include <armadillo>

using namespace arma;

double SingleParticle::Wavefunction()
{
  int i, j;
  double C1, C2, wf, r_single_particle2, hermite1, hermite2;

  r_single_particle2 = 0;
  for (j = 0; j < m_dimension; j++) {
    r_single_particle2 += m_r(j) * m_r(j);
  }

  // ------------------------ hermite's polynomials ------------------------- //
  // for x-direction
  if (m_nx==0){
      hermite1 = 1.;
      C1 = 1.*sqrt(sqrt(m_omega/M_PI));
  }
  if (m_nx==1){
      hermite1 = 2.*sqrt(m_omega)*m_r(0);
      C1 = 0.5*sqrt(sqrt(m_omega/M_PI));
  }
  if (m_nx==2){
      hermite1 = 4.*m_omega*m_r(0)*m_r(0) - 2.;
      C1 = 1./8.*sqrt(sqrt(m_omega/M_PI));
  } 
  // for y-direction
  if (m_ny==0){
      hermite2 = 1.;
      C2 = 1.*sqrt(sqrt(m_omega/M_PI));
  }
  if (m_ny==1){
      hermite2 = 2.*sqrt(m_omega)*m_r(1);
      C2 = 0.5*sqrt(sqrt(m_omega/M_PI));
  }
  if (m_ny==2){
      hermite2 = 4.*m_omega*m_r(1)*m_r(1) - 2.;
      C2 = 1./8.*sqrt(sqrt(m_omega/M_PI));
  }

  // -------------------- single particle wave function --------------------- //
  wf = C1*C2*exp(-m_omega*r_single_particle2*0.5)*hermite1*hermite2;
  return wf;
}

void SingleParticle::SetPosition(vec r) {
    m_r = r;
}

SingleParticle::SingleParticle(vec r,int nx, int ny, int dimension,\
        double omega)
{
    m_r = r;
    m_nx = nx;
    m_ny = ny;
    m_dimension = dimension;
    m_omega = omega;
}
