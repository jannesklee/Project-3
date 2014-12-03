#include "singleparticle.h"
#include <armadillo>

using namespace arma;

// -------------------------- constructors ---------------------------------- //
SingleParticle::SingleParticle()
{
    m_r = 1.0;
    m_nx = 0;
    m_ny = 0;
    m_dimension = 2;
    m_omega = 1.0;
    m_alpha = 1.0; 
}

SingleParticle::SingleParticle(vec r,int nx, int ny, int dimension,\
        double omega, double alpha)
{
    m_r = r;
    m_nx = nx;
    m_ny = ny;
    m_dimension = dimension;
    m_omega = omega;
    m_alpha = alpha; 
}


// -------------------- setters, getters ------------------------------------ //
void SingleParticle::SetPosition(vec r) {
    m_r = r;
}

void SingleParticle::SetAll(vec r, int nx, int ny, int dimension, double omega,\
        double alpha)
{
    m_r = r;
    m_nx = nx;
    m_ny = ny;
    m_dimension = dimension;
    m_omega = omega;
    m_alpha = alpha; 
}


// -------------------------- wavefunctions --------------------------------- //
double SingleParticle::Wavefunction()
{
  int i;
  double wf, r_single_particle2, hermite1, hermite2;

  r_single_particle2 = 0;
  for (i = 0; i < m_dimension; i++) {
    r_single_particle2 += m_r(i) * m_r(i);
  }

  // ------------------------ hermite's polynomials ------------------------- //
  // for x-direction
  // TODO: delete normalizations
  if (m_nx==0){
      hermite1 = 1.;
  }
  if (m_nx==1){
      hermite1 = 2.*sqrt(m_alpha*m_omega)*m_r(0);
  }
  if (m_nx==2){
      hermite1 = 4.*m_alpha*m_omega*m_r(0)*m_r(0) - 2.;
  } 
  // for y-direction
  if (m_ny==0){
      hermite2 = 1.;
  }
  if (m_ny==1){
      hermite2 = 2.*sqrt(m_alpha*m_omega)*m_r(1);
  }
  if (m_ny==2){
      hermite2 = 4.*m_alpha*m_omega*m_r(1)*m_r(1) - 2.;
  }

  // -------------------- single particle wave function --------------------- //
  wf = exp(-m_omega*m_alpha*r_single_particle2*0.5)*hermite1*hermite2;
  return wf;
}

vec SingleParticle::GetGradient() {
    vec slater_grad;
    double c;

    slater_grad = vec(m_dimension); 

    c = m_omega*m_alpha;
    if (m_nx == 0 && m_ny == 0) {
        slater_grad(0) = -c*m_r(0);
        slater_grad(1) = -c*m_r(1);
    }
    else if (m_nx == 1 && m_ny == 0) {
        slater_grad(0) = -(sqrt(c)*m_r(0)-1.)*(sqrt(c)*m_r(0)+1.);
        slater_grad(1) = -c*m_r(0)*m_r(1);
    }
    else if (m_nx == 0 && m_ny == 1) {
        slater_grad(0) = -c*m_r(0)*m_r(1);
        slater_grad(1) = -(sqrt(c)*m_r(1)-1.)*(sqrt(c)*m_r(1)+1.);
    }
    slater_grad = slater_grad*exp(-0.5*c*(m_r(0)*m_r(0)+m_r(1)*m_r(1)));

    return slater_grad;
}

double SingleParticle::GetLaplacian() {
    double slater_lap; 
    double r_single_particle2;
    double c;

    c = m_omega*m_alpha;
    r_single_particle2 = m_r(0)*m_r(0) + m_r(1)*m_r(1);
    slater_lap = 0;
  
    if (m_nx == 0 && m_ny == 0){
        slater_lap = c*(c*r_single_particle2 - 2.);
    }
    else if (m_nx == 1 && m_ny == 0) {
        slater_lap = c*m_r(0)*(c*r_single_particle2 - 4.);
    }
    else if(m_nx == 0 && m_ny == 1) {
        slater_lap = c*m_r(1)*(c*r_single_particle2 - 4.);
    }
    else{
    }

    slater_lap = slater_lap*exp(-0.5*c*(m_r(0)*m_r(0)+m_r(1)*m_r(1)));
    return slater_lap;
}
