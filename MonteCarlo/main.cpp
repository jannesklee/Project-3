// Variational Monte Carlo for atoms with up to two electrons

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "lib.h"
#include "singleparticle.h"  //class for single particles
#include "manybody.h" //class for many-body problems
#include <omp.h>
//#include <mpi.h>

using namespace  std;
using namespace arma;


// output file as global variable
ofstream ofile;
// the step length and its squared inverse for the second derivative
#define h 0.001
#define h2 1000000
#define abegin 0.98
#define bbegin 0.38
#define astep 0.02
#define bstep 0.02

/* -------------------------------------------------------------------------- *
 *                        Declaration of functions                            *
 * -------------------------------------------------------------------------- */
// The Mc sampling for the variational Monte Carlo
void  mc_sampling(int, int, int, int, int, int, double, mat &, mat &, double,\
        int);
// The local energy
double  local_energy(mat, double, double, double, int, int, int, double, int);
// prints to screen the results of the calculations
void  output(int, int, int, mat, mat);
// quantum force for importance sampling
void quantum_force(int, int, double, double, double, double, mat, mat &);
// function for gaussian random numbers
double gaussian_deviate(long *);

/* -------------------------------------------------------------------------- *
 *                        Begin of main program                               *
 * -------------------------------------------------------------------------- */
int main()
{
  int number_cycles = 100000;                 // number of Monte-Carlo steps  //
  int max_variations = 5;                     // max. var. params             //
  int thermalization = 0;                 // Thermalization steps         //
  int charge = 1;                             // nucleus' charge              //
  int dimension = 2;                          // dimensionality               //
  int number_particles = 6;                   // number of particles          //
  double step_length= 0.1;                    // step length                  //
  mat cumulative_e, cumulative_e2;            // energy-matrices              // 
  mat cumulative_e_temp, cumulative_e2_temp;  // energy-matrix (squared)      // 
  double omega = 1.0;                          // freq. harm. osc.             //
  int nx = 0;
  int num_threads;                            // number of threads            // 

  cumulative_e = mat(max_variations+1, max_variations+1);
  cumulative_e2 = mat(max_variations+1, max_variations+1);

  // ----------------------- MC sampling ------------------------------------ //
omp_set_num_threads(1);
#pragma omp parallel shared(cumulative_e_temp, cumulative_e2_temp)
  {
  cumulative_e_temp = mat(max_variations+1, max_variations+1);
  cumulative_e2_temp = mat(max_variations+1, max_variations+1);
  mc_sampling(dimension, number_particles, charge, \
              max_variations, thermalization, number_cycles, \
              step_length, cumulative_e_temp, cumulative_e2_temp, omega, nx);
#pragma omp barrier
#pragma omp critical
  {
    cumulative_e += cumulative_e_temp;
    cumulative_e2 += cumulative_e2_temp;
  }
    num_threads = omp_get_num_threads();
  }
  cout << num_threads << endl;
  cumulative_e = cumulative_e/num_threads;
  cumulative_e2 = cumulative_e2/num_threads;
  

  // ------------------------- Output --------------------------------------- // 
  ofile.open("vmc.dat");
  output(max_variations, number_cycles, charge, cumulative_e, cumulative_e2);
  ofile.close();  // close output file
  return 0;
}


/* -------------------------------------------------------------------------- *
 *             Monte Carlo sampling with the Metropolis algorithm             *
 * -------------------------------------------------------------------------- */
void mc_sampling(int dimension, int number_particles, int charge,
                 int max_variations,
                 int thermalization, int number_cycles, double step_length,
                 mat &cumulative_e, mat &cumulative_e2, double omega, int nx){

  int cycles, variate, variate2, accept, i, j, thread;
  long idum;
  double alpha, beta, energy, energy2, delta_e, wfold, wfnew;
  double D, greensfunction;
  D = 0.5; // 
  alpha = abegin*charge;
  idum=-1;

  // initial positions, initial forces
  mat r_old = zeros<mat>(number_particles, dimension);
  mat r_new = zeros<mat>(number_particles, dimension);
  mat qforce_old = zeros<mat>(number_particles, dimension);
  mat qforce_new = zeros<mat>(number_particles, dimension);
   
  // -------------- Loop over different values of alpha, beta --------------- //
  for (variate = 1; variate <= max_variations; variate++){
      alpha += astep;
      beta = bbegin*charge;
      for (variate2 = 1; variate2 <= max_variations; variate2++){ 
          beta += bstep;
          energy = energy2 = 0; accept = 0; delta_e = 0;

          //  initial trial position
          for (i = 0; i < number_particles; i++) {
            for (j = 0; j < dimension; j++) {
//              r_old(i,j) = step_length*(ran2(&idum)-0.5);
              r_old(i,j) = gaussian_deviate(&idum);//*sqrt(step_length);
            }
          }


          /*SingleParticle particle_old(r_old, nx, dimension,\
                                          number_particles, omega);*/
          ManyBody system_old(r_old, alpha, beta, dimension,\
                                      number_particles, omega);
         
          // clearify which wavefunction shall be used: perturbed or unperturbed 
          //wfold = particle_old.PerturbedWavefunction();
          //wfold= particle_old.UnperturbedWavefunction();
          wfold = system_old.SixElectronSystem();

          quantum_force(number_particles, dimension, alpha, beta, omega, \
                  wfold, r_old, qforce_old);
          
          // -------------- loop over monte carlo cycles -------------------- //
          for (cycles = 1; cycles <= number_cycles+thermalization; cycles++){
            // new position
            for (i = 0; i < number_particles; i++) {
              for (j = 0; j < dimension; j++) {
//                r_new(i,j) = r_old(i,j) + step_length*(ran1(&idum)-0.5);
                r_new(i,j) = r_old(i,j) + gaussian_deviate(&idum)*
                             sqrt(step_length) + step_length*D*qforce_old(i,j);// 
              }
            }

            /* SingleParticle particle_new(r_new, nx, dimension,\
                                            number_particles, omega);*/
            ManyBody system_new(r_new, alpha, beta, dimension,\
                                            number_particles, omega);

            // clearify which wavefunction shall be used: perturbed or unpert.
            //wfnew = system_new.PerturbedWavefunction();
            wfnew = system_new.SixElectronSystem();
            //cout << wfnew << endl; 
            //wfnew= particle_new.UnperturbedWavefunction();
           
            quantum_force(number_particles, dimension, alpha, beta, omega,\
                    wfnew, r_new, qforce_new);
            
            // ------------------ greensfunction ---------------------------- //
            greensfunction = 0.0;
            for (i = 0; i < number_particles; i++) {
                for (j = 0; j < dimension; j++) {
                    greensfunction += 0.5*(qforce_old(i,j) + qforce_new(i,j))* \
                       (D*step_length*0.5*(qforce_old(i,j) - qforce_new(i,j))- \
                        r_new(i,j) + r_old(i,j));
                }
            }
            greensfunction = exp(greensfunction);

            // ----------------- metropolis test ---------------------------- //
            if (ran1(&idum) <= greensfunction*wfnew*wfnew/wfold/wfold){
                for (i = 0; i < number_particles; i++) {
                    for (j = 0; j < dimension; j++){
                        r_old(i,j) = r_new(i,j);
                        qforce_old(i,j) = qforce_new(i,j); 
                    }
                }
                wfold = wfnew;
                accept = accept + 1;
            }


            // ----------------- local energy ------------------------------- //
            if (cycles > thermalization) {
              delta_e = local_energy(r_old, alpha, beta, wfold, dimension,
                                     number_particles, charge, omega, nx);
              // update energies
              energy += delta_e;
              energy2 += delta_e*delta_e;
            }
          } 
          // ------------- end loop over monte carlo cycles ----------------- //    
           
          // update the energy average and its squared
          cumulative_e(variate, variate2) = energy/number_cycles;
          cumulative_e2(variate, variate2) = energy2/number_cycles;

          thread = omp_get_thread_num();
#pragma omp critical
          {
          cout << "alpha = " << alpha << setw(15)
               << "beta = " << beta << setw(20)
               << "accepted steps = " << accept << setw(20)
               << "energy = " << cumulative_e(variate,variate2) << setw(20)
               << "thread = " << thread << endl;
          }
      }
  }    // end of loop over variational  steps
}   // end mc_sampling function


void quantum_force(int number_particles, int dimension, double alpha, \
        double beta, double omega, double wf, mat r, mat &qforce)
{
    int i, j;
    double wfminus , wfplus;
    mat r_plus, r_minus;

    r_plus = zeros<mat>(number_particles,dimension);
    r_minus = zeros<mat>(number_particles,dimension);

    for(i = 0 ; i < number_particles; i++) {
        for (j = 0; j < dimension ; j++) {
            r_plus(i,j) = r_minus(i,j) = r(i,j);
        }
    }

    // quantum the first derivative
    for (i = 0; i < number_particles; i++){
        for (j = 0; j < dimension; j++){
            r_plus(i,j) = r(i,j) + h;
            r_minus(i,j) = r(i,j) - h; 
            ManyBody particle_minus(r_minus, alpha, beta, dimension,\
                                            number_particles, omega);
            wfminus = particle_minus.SixElectronSystem();
            ManyBody particle_plus(r_plus, alpha, beta, dimension,\
                                            number_particles, omega);
            wfplus = particle_plus.SixElectronSystem();
            qforce(i,j) = (wfplus - wfminus)/(wf*h);
            r_plus(i,j) = r(i,j);
            r_minus(i,j) = r(i,j);
        }
    }
}

/* -------------------------------------------------------------------------- *
 *        Function to calculate the local energy with num derivative          *
 * -------------------------------------------------------------------------- */
double  local_energy(mat r, double alpha, double beta, double wfold,\
            int dimension, int number_particles, int charge, double omega, int nx)
{
  int i, j , k;
  double e_local, e_kinetic, e_potential, r_12, r_single_particle;
  double wfminus, wfplus;
  mat r_plus,r_minus;

  // allocate matrices which contain the position of the particles
  // the function matrix is defined in the progam library
  r_plus = zeros<mat>(number_particles,dimension);
  r_minus = zeros<mat>(number_particles,dimension);
  for (i = 0; i < number_particles; i++) {
    for (j = 0; j < dimension; j++) {
      r_plus(i,j) = r_minus(i,j) = r(i,j);
    }
  }

  // ---------------------- kinetic energy ---------------------------------- //
  e_kinetic = 0;
  for (i = 0; i < number_particles; i++) {
    for (j = 0; j < dimension; j++) {
      r_plus(i,j) = r(i,j) + h;
      r_minus(i,j) = r(i,j) - h;
      // SingleParticle particle_minus(r_minus, nx, dimension,
      //                                number_particles, omega);
      ManyBody particle_minus(r_minus, alpha, beta, dimension,\
                              number_particles, omega);
      wfminus= particle_minus.SixElectronSystem();
      // double wfminus= particle_minus.UnperturbedWavefunction();
      // SingleParticle particle_plus(r_plus, nx, dimension,
      //                                number_particles, omega);
      ManyBody particle_plus(r_plus, alpha, beta, dimension,\
                                           number_particles, omega);
      wfplus= particle_plus.SixElectronSystem();
      // double wfplus= particle_plus.UnperturbedWavefunction();
      e_kinetic -= (wfminus + wfplus - 2*wfold);
      r_plus(i,j) = r(i,j);
      r_minus(i,j) = r(i,j);
    }
  }
  // include electron mass and hbar squared and divide by wave function
  e_kinetic = 0.5*h2*e_kinetic/wfold;

  // ---------------------- potential energy -------------------------------- //
  e_potential = 0;
  // contribution from electron-proton potential
  for (i = 0; i < number_particles; i++) {
    r_single_particle = 0;
    for (j = 0; j < dimension; j++) {
      r_single_particle += r(i,j)*r(i,j);
    }
    e_potential += 0.5*omega*omega*r_single_particle; 
  }
  // contribution from electron-electron potential  
  for (i = 0; i < number_particles-1; i++) {
    for (j = i+1; j < number_particles; j++) {
      r_12 = 0;
      for (k = 0; k < dimension; k++) {
        r_12 += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
      }
      e_potential += 1/sqrt(r_12);
    }
  }

  // -------------------------- total energy -------------------------------- //
  e_local = e_potential + e_kinetic;
  return e_local;
}


// output function
void output(int max_variations, int number_cycles, int charge,
            mat cumulative_e, mat cumulative_e2)
{
  int i, j;
  double alpha, beta, variance, error;
  alpha = abegin*charge;   
  for(i = 1; i <= max_variations; i++){
      alpha += astep;   
      beta = bbegin;           
      for (j = 1; j <= max_variations; j++){
          beta += bstep;
          variance = cumulative_e2(i,j)-cumulative_e(i,j)*cumulative_e(i,j);
          error=sqrt(variance/number_cycles);
          ofile << setiosflags(ios::showpoint | ios::uppercase);
          ofile << setw(15) << setprecision(8) << alpha;
          ofile << setw(15) << setprecision(8) << beta;
          ofile << setw(15) << setprecision(8) << cumulative_e(i,j);
          ofile << setw(15) << setprecision(8) << variance;
          ofile << setw(15) << setprecision(8) << error << endl;
      }
      ofile << endl;
  }
} 


///* -------------------------------------------------------------------------- *
// *                    Random number generator                                 *
// * -------------------------------------------------------------------------- */
//#define IA 16807
//#define IM 2147483647
//#define AM (1.0/IM)
//#define IQ 127773
//#define IR 2836
//#define NTAB 32
//#define NDIV (1+(IM-1)/NTAB)
//#define EPS 1.2e-7
//#define RNMX (1.0-EPS)
//
//double ran1(long *idum)
//{
//   int             j;
//   long            k;
//   static long     iy=0;
//   static long     iv[NTAB];
//   double          temp;
//
//   if (*idum <= 0 || !iy) {
//      if (-(*idum) < 1) *idum=1;
//      else *idum = -(*idum);
//      for(j = NTAB + 7; j >= 0; j--) {
//         k     = (*idum)/IQ;
//         *idum = IA*(*idum - k*IQ) - IR*k;
//         if(*idum < 0) *idum += IM;
//         if(j < NTAB) iv[j] = *idum;
//      }
//      iy = iv[0];
//   }
//   k     = (*idum)/IQ;
//   *idum = IA*(*idum - k*IQ) - IR*k;
//   if(*idum < 0) *idum += IM;
//   j     = iy/NDIV;
//   iy    = iv[j];
//   iv[j] = *idum;
//   if((temp=AM*iy) > RNMX) return RNMX;
//   else return temp;
//}

