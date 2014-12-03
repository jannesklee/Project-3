// Variational Monte Carlo for atoms with up to six electrons
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <omp.h>
#include <unittest++/UnitTest++.h>
#include "lib.h"
#include "singleparticle.h"  //class for single particles
#include "manybody.h" //class for many-body problems
#include "jastrow.h" //class for derivations of jastrow factor
#include "slater.h"
//#include <mpi.h>

using namespace  std;
using namespace arma;

// output file as global variable
ofstream ofile;
// the step length and its squared inverse for the second derivative
#define h 0.001
#define h2 1000000
#define abegin 0.9
#define bbegin 0.3
#define astep 0.04
#define bstep 0.04

/* -------------------------------------------------------------------------- *
 *                        Declaration of functions                            *
 * -------------------------------------------------------------------------- */
// The Mc sampling for the variational Monte Carlo
void mc_sampling(int, int, int, int, int, int, double, mat &, mat &, double,\
        mat &, mat &);
// The local energy
double local_energy(mat, double, double, double, int, int, int, double,\
        double &, double &, int);
// prints to screen the results of the calculations
void output(int, int, int, mat, mat, mat, mat);
// quantum force for importance sampling
void quantum_force(int, int, double, double, double, double, mat, mat &, int);
// qunatum force init
void quantum_force_init(int, int, double, double, double, double, mat, mat &);

/* -------------------------------------------------------------------------- *
 *                        Begin of main program                               *
 * -------------------------------------------------------------------------- */
int main()
{
  int number_cycles = 100000;                 // number of Monte-Carlo steps  //
  int max_variations = 5;                     // max. var. params             //
  int thermalization = 0; 
  int charge = 1;                             // nucleus' charge              //
  int dimension = 2;                          // dimensionality               //
  int number_particles = 2;                   // number of particles          //
  double step_length= 0.1;                    // either f. br.for. or imp.samp//
  mat cumulative_e, cumulative_e2;            // energy-matrices              //
  mat cumulative_e_temp, cumulative_e2_temp;  // energy-matrix (squared)      //
  mat kin_e, pot_e;
  mat kin_e_temp, pot_e_temp;
  double omega = 1.;                         // freq. harm. osc.             //
  int num_threads;                            // number of threads            //

  cumulative_e = mat(max_variations+1, max_variations+1, fill::zeros);
  cumulative_e2 = mat(max_variations+1, max_variations+1, fill::zeros);
  kin_e = mat(max_variations+1, max_variations+1, fill::zeros);
  pot_e = mat(max_variations+1, max_variations+1, fill::zeros);


  // ----------------------- MC sampling ------------------------------------ //
//omp_set_num_threads(1);
#pragma omp parallel shared(cumulative_e_temp, cumulative_e2_temp)
  {
  cumulative_e_temp = mat(max_variations+1, max_variations+1, fill::zeros);
  cumulative_e2_temp = mat(max_variations+1, max_variations+1, fill::zeros);
  kin_e_temp = mat(max_variations+1, max_variations+1, fill::zeros);
  pot_e_temp = mat(max_variations+1, max_variations+1, fill::zeros);
  mc_sampling(dimension, number_particles, charge, \
              max_variations, number_cycles, thermalization, \
              step_length, cumulative_e_temp, cumulative_e2_temp, omega, \
              kin_e_temp, pot_e_temp);
#pragma omp barrier
#pragma omp critical
  {
    cumulative_e += cumulative_e_temp;
    cumulative_e2 += cumulative_e2_temp;
    kin_e += kin_e_temp;
    pot_e += pot_e_temp;
  }
    num_threads = omp_get_num_threads();
  }

  cout << num_threads << endl;
  cumulative_e = cumulative_e/num_threads;
  cumulative_e2 = cumulative_e2/num_threads;
  kin_e = kin_e/num_threads;
  pot_e = pot_e/num_threads;


  // ------------------------- Output --------------------------------------- //
  ofile.open("vmc.dat");
  output(max_variations, number_cycles, charge, cumulative_e, cumulative_e2,\
          kin_e, pot_e);
  ofile.close();  // close output file


  // ------------------------ unittests ------------------------------------- //
  //return UnitTest::RunAllTests();
  return 0;
}


/* -------------------------------------------------------------------------- *
 *             Monte Carlo sampling with the Metropolis algorithm             *
 * -------------------------------------------------------------------------- */
void mc_sampling(int dimension, int number_particles, int charge,
                 int max_variations, int number_cycles, int thermalization, double step_length,
                 mat &cumulative_e, mat &cumulative_e2, double omega,
                 mat &kin_e, mat &pot_e){

  int cycles, variate, variate2, accept, i, j, k, thread;
  long idum;
  double alpha, beta, energy, energy2, delta_e, wfold, wfnew;
  double D, greensfunction, del_kin_e, del_pot_e;
  double kinetic_energy, potential_energy;

  D = 0.5; 
  idum=-1;

  // initial positions, initial forces
  mat r_old = zeros<mat>(number_particles, dimension);
  mat r_new = zeros<mat>(number_particles, dimension);
  mat qforce_old = zeros<mat>(number_particles, dimension);
  mat qforce_new = zeros<mat>(number_particles, dimension);
  ManyBody system(dimension, number_particles, omega);

  alpha = abegin*charge;
  // -------------- Loop over different values of alpha, beta --------------- //
  for (variate = 1; variate <= max_variations; variate++){
      beta = bbegin*charge;
      for (variate2 = 1; variate2 <= max_variations; variate2++){
          energy = energy2 = kinetic_energy = potential_energy = 0; 
          accept = 0; delta_e = 0;
          system.SetVariables(alpha, beta);

          //  initial trial position
          for (i = 0; i < number_particles; i++) {
            for (j = 0; j < dimension; j++) {
//             r_old(i,j) = step_length*(ran2(&idum)-0.5);
              r_old(i,j) = gaussian_deviate(&idum);//*sqrt(step_length);
            }
          }

          system.SetPosition(r_old);
          wfold = system.SixElectronSystem();
          //wfold = particle_old.UnperturbedWavefunction();

          i = 0; 
          quantum_force_init(number_particles, dimension, alpha, beta, omega, \
                  wfold, r_old, qforce_old);

          // -------------- loop over monte carlo cycles -------------------- //
          for (cycles = 1; cycles <= number_cycles + thermalization; cycles++){
            // new position
            for (i = 0; i < number_particles; i++) {
              for (j = 0; j < dimension; j++) {
//              r_new(i,j) = r_old(i,j) + step_length*(ran1(&idum)-0.5);
                r_new(i,j) = r_old(i,j) + gaussian_deviate(&idum)*sqrt(step_length)
                           + step_length*D*qforce_old(i,j);
              }
              
              // we move only one particle at the time
              for (k = 0; k < number_particles; k++) {
                  if (k != i) { 
                      for (j = 0; j < dimension; j++) {  // resets all elements to old
                          r_new(k,j) = r_old(k,j);
                      }
                  }
              }

              system.SetPosition(r_new);
              wfnew = system.SixElectronSystem();
              
              quantum_force(number_particles, dimension, alpha, beta, omega,\
                      wfnew, r_new, qforce_new, i);
              
              // ------------------ greensfunction ---------------------------- //
              greensfunction = 0.0;
              for (j = 0; j < dimension; j++) {
                  greensfunction += 0.5*(qforce_old(i,j) + qforce_new(i,j))* \
                     (D*step_length*0.5*(qforce_old(i,j) - qforce_new(i,j))- \
                      r_new(i,j) + r_old(i,j));
              }
              greensfunction = exp(greensfunction);

//              greensfunction = 1.;
              // ----------------- metropolis test ---------------------------- //
              if (ran2(&idum) <= greensfunction*wfnew*wfnew/wfold/wfold){
                  for (j = 0; j < dimension; j++){
                      r_old(i,j) = r_new(i,j);
                      qforce_old(i,j) = qforce_new(i,j); 
                  }
              wfold = wfnew;
              accept = accept + 1;
              }
            }

            // ----------------- local energy ------------------------------- //
            if (cycles > thermalization){
                delta_e = local_energy(r_old, alpha, beta, wfold, dimension,
                                       number_particles, charge, omega, \
                                       del_kin_e, del_pot_e, i);
                // update energies
                energy += delta_e;
                energy2 += delta_e*delta_e;
                kinetic_energy += del_kin_e;
                potential_energy += del_pot_e;
            }
          }
          // ------------- end loop over monte carlo cycles ----------------- //

          // update the energy average and its squared
          cumulative_e(variate, variate2) = energy/number_cycles;
          cumulative_e2(variate, variate2) = energy2/number_cycles;
          kin_e(variate, variate2) = kinetic_energy/number_cycles;
          pot_e(variate, variate2) = potential_energy/number_cycles;

          thread = omp_get_thread_num();
#pragma omp critical
          {
          cout << "alpha = " << alpha << setw(15)
               << "beta = " << beta << setw(20)
               << "accepted steps = " << accept << setw(20)
               << "energy = " << cumulative_e(variate,variate2) << setw(20)
               << "kin. energy = " << kin_e(variate,variate2) << setw(20)
               << "pot. energy = " << pot_e(variate,variate2) << setw(20)
               << "thread = " << thread << endl;
          }
          beta += bstep;
      }
      alpha += astep;
  }    // end of loop over variational  steps
}   // end mc_sampling function


void quantum_force(int number_particles, int dimension, double alpha, \
        double beta, double omega, double wf, mat r, mat &qforce, \
        int particle_ind)
{
    int k;
    vec grad_slater, grad_jastrow;
    (void) wf; 

    // Setup Slater Determinant
    Slater slater_obj(r, alpha, beta, dimension, number_particles, omega);
    slater_obj.SetupSixElectron();
    grad_slater = slater_obj.Gradient(particle_ind); 

    // Setup Jastrow
    Jastrow jastrow_obj(r, alpha, beta, dimension, number_particles, omega);
    grad_jastrow = jastrow_obj.Gradient(particle_ind);

    for(k = 0; k < dimension; k++) {
        qforce(particle_ind,k) = 2.*(grad_slater(k) + grad_jastrow(k));
    }
}

void quantum_force_init(int number_particles, int dimension, double alpha, \
        double beta, double omega, double wf, mat r, mat &qforce)
{
    int k, particle_ind;
    vec grad_slater, grad_jastrow;
    (void) wf; 

    // Setup Slater
    Slater slater_obj(r, alpha, beta, dimension, number_particles, omega);
    slater_obj.SetupSixElectron();

    // Setup Jastrow
    Jastrow jastrow_obj(r, alpha, beta, dimension, number_particles, omega);

    for (particle_ind = 0; particle_ind < number_particles; particle_ind++) {
        grad_slater = slater_obj.Gradient(particle_ind); 
        grad_jastrow = jastrow_obj.Gradient(particle_ind);

        for(k = 0; k < dimension; k++) {
            qforce(particle_ind,k) = 2.*(grad_slater(k) + grad_jastrow(k));
        }
    }
}

/* -------------------------------------------------------------------------- *
 *        Function to calculate the local energy with num derivative          *
 * -------------------------------------------------------------------------- */
double local_energy(mat r, double alpha, double beta, double wfold,\
          int dimension, int number_particles, int charge, double omega,\
          double &kin_e, double &pot_e, int particle)
{
  int i, j , k;
  double e_local, e_kinetic, e_potential, r_12, r_single_particle;
  vec grad_slater, grad_jastrow;
  double lap_slater, lap_jastrow;
  (void) charge;
  (void) particle; 
  (void) wfold;


  // ---------------------- kinetic energy ---------------------------------- //
  // closed-form solution 
   
  // Setup Slater
  Slater slater_obj(r, alpha, beta, dimension, number_particles, omega);
  slater_obj.SetupSixElectron();

  // Setup Jastrow
  Jastrow jastrow_obj(r, alpha, beta, dimension, number_particles, omega);


  e_kinetic = 0.0;
  for (i = 0; i < number_particles; i++) {
      grad_slater = slater_obj.Gradient(i); 
      lap_slater = slater_obj.Laplacian(i);
      grad_jastrow = jastrow_obj.Gradient(i); // gradient has to be calculated befor laplacian
      lap_jastrow = jastrow_obj.Laplacian(i);

      e_kinetic += -0.5*(lap_slater + lap_jastrow + 2*dot(grad_slater,grad_jastrow));
  }


//  // allocate matrices which contain the position of the particles
//  // the function matrix is defined in the progam library
//  r_plus = zeros<mat>(number_particles,dimension);
//  r_minus = zeros<mat>(number_particles,dimension);
//  for (i = 0; i < number_particles; i++) {
//    for (j = 0; j < dimension; j++) {
//      r_plus(i,j) = r_minus(i,j) = r(i,j);
//    }
//  }

//  // brute force derivation
//  e_kinetic = 0;
//  for (i = 0; i < number_particles; i++) {
//    for (j = 0; j < dimension; j++) {
//      r_plus(i,j) = r(i,j) + h;
//      r_minus(i,j) = r(i,j) - h;
//
//      system.SetPosition(r_minus);
//      //wfminus = system.PerturbedWavefunction();
//      wfminus = system.SixElectronSystem();
//      system.SetPosition(r_plus);
//      //wfplus = system.PerturbedWavefunction();
//      wfplus = system.SixElectronSystem();
//
//      e_kinetic -= (wfminus + wfplus - 2.*wfold);
//      r_plus(i,j) = r(i,j);
//      r_minus(i,j) = r(i,j);
//    }
//  }
//  // include electron mass and hbar squared and divide by wave function
//  e_kinetic = 0.5*h2*e_kinetic/wfold;

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
  kin_e = e_kinetic;
  pot_e = e_potential;
  e_local = e_potential + e_kinetic;
  return e_local;
}


// output function
void output(int max_variations, int number_cycles, int charge, \
        mat cumulative_e, mat cumulative_e2, mat kin_e, mat pot_e)
{
  int i, j;
  double alpha, beta, variance, error;
  alpha = abegin*charge;
  ofile << setw(15) << "alpha";
  ofile << setw(15) << "beta";
  ofile << setw(15) << "cumulative_e(i,j)";
  ofile << setw(15) << "kinetic energy";
  ofile << setw(15) << "potential energy"; 
  ofile << setw(15) << "variance (cum_e)";
  ofile << setw(15) << "error (cum_e)" << endl;
  for(i = 1; i <= max_variations; i++){
      beta = bbegin;
      for (j = 1; j <= max_variations; j++){
          variance = cumulative_e2(i,j)-cumulative_e(i,j)*cumulative_e(i,j);
          error=sqrt(variance/number_cycles);
          ofile << setiosflags(ios::showpoint | ios::uppercase);
          ofile << setw(15) << setprecision(8) << alpha;
          ofile << setw(15) << setprecision(8) << beta;
          ofile << setw(15) << setprecision(8) << cumulative_e(i,j);
          ofile << setw(15) << setprecision(8) << kin_e(i,j);
          ofile << setw(15) << setprecision(8) << pot_e(i,j);
          ofile << setw(15) << setprecision(8) << variance;
          ofile << setw(15) << setprecision(8) << error << endl;
          beta += bstep;
      }

      alpha += astep;
  }
}
