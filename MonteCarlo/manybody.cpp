#include "singleparticle.h"
#include "manybody.h"
#include <armadillo>
#include <iomanip>

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

  // calculates the absolute position
  argument = wf = 0;
  for (i = 0; i < m_number_particles; i++) {
    r_single_particle = 0;
    for (j = 0; j < m_dimension; j++) {
      r_single_particle  += m_r(i,j)*m_r(i,j);
    }
    argument += r_single_particle; 
  }

  // calculates the relative distance // TODO: At that point too much calculation (loop)
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
    // TODO: the matrix m_r is not coupled to the vecotrs m_r in singleparticle
    // this might get into trouble if only one of these things is done
    int i, j, k, red_slater_size;
    double a, wf, argument, r_single_particle, r_12;
    double psi_c;
    mat slater_up, slater_down;
    red_slater_size = m_number_particles/2;  // needs an even amount of particles
    vec m_r1, m_r2, m_r3, m_r4, m_r5, m_r6;
    SingleParticle particle[6];
    double det_slater_up, det_slater_down;

    slater_up = mat(red_slater_size, red_slater_size);  
    slater_down = mat(red_slater_size, red_slater_size);

    // calculates the absolute position
    argument = wf = 0;
    for (i = 0; i < m_number_particles; i++) {
        r_single_particle = 0;
        for (j = 0; j < m_dimension; j++) {
          r_single_particle  += m_r(i,j)*m_r(i,j);
        }
        argument += r_single_particle; 
    }
    //cout << "-- calculated absolute position" << endl;
    //cout << argument << endl;

    // calculates the relative distance 
    psi_c = 1.0;
    for (i = 0; i < m_number_particles-1; i++) { 
        for (j = i+1; j < m_number_particles; j++) {
            r_12 = 0;
            for (k = 0; k < m_dimension; k++) {
                r_12 += (m_r(i,k)-m_r(j,k))*(m_r(i,k)-m_r(j,k));
            }
                // evaluate a for parallel or antiparallel spin
                if ((i+j)%2 == 0) { 
                    a = 1./3.;
                }
                else {
                    a = 1.;
                }
            r_12 = sqrt(r_12);
            psi_c = psi_c * exp(a*r_12/(1. + m_beta*r_12));
        }
    }
    //cout << "-- calculated relative distances" << endl;
    //cout << r_12 << setw(15) << psi_c << endl;

    // spin-up particles
    particle[0].SetAll(conv_to<vec>::from(m_r.row(0)), 0, 0, m_dimension,\
            m_omega);
    particle[2].SetAll(conv_to<vec>::from(m_r.row(0)), 1, 0, m_dimension,\
            m_omega);
    particle[4].SetAll(conv_to<vec>::from(m_r.row(0)), 0, 1, m_dimension,\
            m_omega);
    // spin-down particles
    particle[1].SetAll(conv_to<vec>::from(m_r.row(3)), 0, 0, m_dimension,\
            m_omega);
    particle[3].SetAll(conv_to<vec>::from(m_r.row(3)), 1, 0, m_dimension,\
            m_omega);
    particle[5].SetAll(conv_to<vec>::from(m_r.row(3)), 0, 1, m_dimension,\
           m_omega);

    //cout << "-- setup particles" << endl;

    // filling of slater matrix
    for (i = 0; i < red_slater_size; i++) {
        for (j = 0; j < red_slater_size; j++) {
            particle[i*2].SetPosition(conv_to<vec>::from(m_r.row(j))); // spin_up
            slater_up(i,j) = particle[i*2].Wavefunction();
            particle[i*2+1].SetPosition(conv_to<vec>::from(m_r.row(j+3))); // spin_down
            slater_down(i,j) = particle[i*2+1].Wavefunction();
        }
    }
    //cout << "-- filled slater matrix" << endl;
    //cout << slater_up << endl << slater_down << endl; 

    // calculate slater determinant applying sarrus' rule
    det_slater_up = det_slater_down = 0;
    for (j = 0; j < red_slater_size; j++) {
        det_slater_up += slater_up(0,j)*slater_up(1,(j+1)%3)*\
                         slater_up(2,(j+2)%3);
        det_slater_up -= slater_up(2,j)*slater_up(1,(j+1)%3)*\
                         slater_up(0,(j+2)%3);

        det_slater_down += slater_down(0,j)*slater_down(1,(j+1)%3)*\
                         slater_down(2,(j+2)%3);
        det_slater_down -= slater_down(2,j)*slater_down(1,(j+1)%3)*\
                         slater_down(0,(j+2)%3);
    }
    //cout << "-- calculated determinant" << endl;
    //cout << det_slater_up << setw(15) << det_slater_down << endl;

    // ---------------------- calculate wavefunction ----------------------- //
    wf = det_slater_up*det_slater_down*psi_c;
    //cout << "wavefunction = " << wf << endl;
    return wf;
}

void ManyBody::SetPosition(mat r) {
    m_r = r;
}


ManyBody::ManyBody(mat r, double alpha, double beta, int dimension, int number_particles, double omega)
{
    m_r = r;
    m_alpha = alpha;
    m_beta = beta;
    m_dimension = dimension;
    m_number_particles = number_particles;
    m_omega = omega;
}
