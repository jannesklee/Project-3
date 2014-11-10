// Variational Monte Carlo for atoms with up to two electrons

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
//#include <mpi.h>

using namespace  std;
using namespace arma;


// output file as global variable
ofstream ofile;
// the step length and its squared inverse for the second derivative
#define h 0.001
#define h2 1000000
#define abegin 0.7
#define bbegin 0.1
#define astep 0.01
#define bstep 0.01

/* -------------------------------------------------------------------------- *
 *                        Declaration of functions                            *
 * -------------------------------------------------------------------------- */
// The Mc sampling for the variational Monte Carlo
void  mc_sampling(int, int, int, int, int, int, double, mat &, mat &);
// The variational wave function
// TODO: implementation with classes 
//double  wave_function(mat, double, int, int);  // with repulsion
double  wave_function(mat, double, double, int, int);
// The local energy
double  local_energy(mat, double, double, double, int, int, int);
// prints to screen the results of the calculations
void  output(int, int, int, mat, mat);
// pseudo-random numbers generator
double ran1(long *);

/* -------------------------------------------------------------------------- *
 *                        Begin of main program                               *
 * -------------------------------------------------------------------------- */
int main()
{
  int number_cycles = 1000000;                // number of Monte-Carlo steps  //
  int max_variations = 50;                     // max. var. params             //
  int thermalization = 1000000;                // Thermalization steps         //
  int charge = 1;                             // nucleus' charge              //
  int dimension = 2;                          // dimensionality               //
  int number_particles = 2;                   // number of particles          //
  double step_length= 1.0;                    // step length                  //
  mat cumulative_e;                           // energy-matrix                // 
  mat cumulative_e2;                          // energy-matrix (squared)      // 

  cumulative_e = mat(max_variations+1, max_variations+1);
  cumulative_e2 = mat(max_variations+1, max_variations+1);

  ofile.open("vmc.dat");
  // ----------------------- MC sampling ------------------------------------ //
  mc_sampling(dimension, number_particles, charge, \
              max_variations, thermalization, number_cycles, \
              step_length, cumulative_e, cumulative_e2);
  // ------------------------- Output --------------------------------------- // 
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
                 mat &cumulative_e, mat &cumulative_e2){
  int cycles, variate, variate2, accept, i, j;
  __attribute__((unused)) int dim;
  long idum;
  double wfnew, wfold, alpha, beta, energy, energy2, delta_e;
  alpha = abegin*charge;
  idum=-1;

  // initial positions
  mat r_old = zeros<mat>(number_particles, dimension);
  mat r_new = zeros<mat>(number_particles, dimension);
   
  // -------------- Loop over different values of alpha, beta --------------- //
  for (variate = 1; variate <= max_variations; variate++){
      alpha += astep;
      beta = bbegin*charge;
      for (variate2 = 1; variate2 <= max_variations; variate2++){ 
          beta += bstep;
          energy = energy2 = 0; accept = 0; delta_e = 0;

          //  initial trial position
          for (i = 0; i < number_particles; i++) {
            for ( j = 0; j < dimension; j++) {
              r_old(i,j) = step_length*(ran1(&idum)-0.5);
            }
          }

          wfold = wave_function(r_old, alpha, beta, dimension,\
                  number_particles);
          
          // -------------- loop over monte carlo cycles -------------------- //
          for (cycles = 1; cycles <= number_cycles+thermalization; cycles++){
            // new position
            for (i = 0; i < number_particles; i++) {
              for (j = 0; j < dimension; j++) {
                r_new(i,j) = r_old(i,j) + step_length*(ran1(&idum)-0.5);
              }
            }
            wfnew = wave_function(r_new, alpha, beta, dimension,\
                    number_particles);
            
            // metropolis test
            if(ran1(&idum) <= wfnew*wfnew/wfold/wfold ) {
              for (i = 0; i < number_particles; i++) {
                for ( j = 0; j < dimension; j++) {
                  r_old(i,j)=r_new(i,j);
                }
              }
              wfold = wfnew;
              accept = accept + 1;
            }

            // compute local energy
            if (cycles > thermalization) {
              delta_e = local_energy(r_old, alpha, beta, wfold, dimension,
                                     number_particles, charge);
              // update energies
              energy += delta_e;
              energy2 += delta_e*delta_e;
            }
          }   

          cout << "alpha = " << alpha << setw(20)
               << "beta = " << beta << setw(20)
               << "accepted steps = " << accept << endl;
          // update the energy average and its squared
          cumulative_e(variate, variate2) = energy/number_cycles;
          cumulative_e2(variate, variate2) = energy2/number_cycles;
      }
  }    // end of loop over variational  steps
}   // end mc_sampling function


// Function to compute the squared wave function, simplest form
//double  wave_function(mat r, double alpha,int dimension, int number_particles) //this function is limited to two electrons
//{
//  int i, j, k;
//  double wf, argument, r_single_particle, r_12;
//
//  argument = wf = 0;
//  for (i = 0; i < number_particles; i++) {
//    r_single_particle = 0;
//    for (j = 0; j < dimension; j++) {
//      r_single_particle  += r(i,j)*r(i,j);
//    }
//    argument += sqrt(r_single_particle);
//  }
//  wf = exp(-argument*alpha) ;
//  return wf;
//}

/* -------------------------------------------------------------------------- *
 *        Function to compute the squared wave function, simplest form        *
 * -------------------------------------------------------------------------- */
double  wave_function(mat r, double alpha, double beta, int dimension, \
        int number_particles) 
{
  int i, j, k;
  double a, C;
  double wf, argument, omega, r_single_particle, r_12;
  //TODO: implement a and C in a proper way
  a = 1.0; // antiparallel spin
  C = 1.0; // normalization
  omega = 1.0;

  argument = wf = 0;
  for (i = 0; i < number_particles; i++) {
    r_single_particle = 0;
    for (j = 0; j < dimension; j++) {
      r_single_particle  += r(i,j)*r(i,j);
    }
    argument += r_single_particle; //TODO: check if removed squarroot is reasonable
  }

  // TODO: copied from below 
  for (i = 0; i < number_particles-1; i++) { // for 2 electrons no loop
    for (j = i+1; j < number_particles; j++) {
      r_12 = 0;
      for (k = 0; k < dimension; k++) {
        r_12 += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
      }
    }
  }

  wf = C*exp(-alpha*omega*argument*0.5)*exp(a*r_12/(1.+ beta*r_12));
  return wf;
}

/* -------------------------------------------------------------------------- *
 *        Function to calculate the local energy with num derivative          *
 * -------------------------------------------------------------------------- */
double  local_energy(mat r, double alpha, double beta, double wfold,\
            int dimension, int number_particles, int charge)
{
  int i, j , k;
  double e_local, wfminus, wfplus, e_kinetic, e_potential, r_12, \
    r_single_particle;
  mat r_plus,r_minus;

  // allocate matrices which contain the position of the particles
  // the function matrix is defined in the progam library
  r_plus = zeros<mat>(number_particles,dimension);
  r_minus = zeros<mat>(number_particles,dimension);
  for (i = 0; i < number_particles; i++) {
    for (j=0; j < dimension; j++) {
      r_plus(i,j) = r_minus(i,j) = r(i,j);
    }
  }

  // compute the kinetic energy
  e_kinetic = 0;
  for (i = 0; i < number_particles; i++) {
    for (j = 0; j < dimension; j++) {
      r_plus(i,j) = r(i,j) + h;
      r_minus(i,j) = r(i,j) - h;
      wfminus = wave_function(r_minus, alpha, beta, dimension, number_particles);
      wfplus  = wave_function(r_plus, alpha, beta, dimension, number_particles);
      e_kinetic -= (wfminus + wfplus - 2*wfold);
      r_plus(i,j) = r(i,j);
      r_minus(i,j) = r(i,j);
    }
  }
  // include electron mass and hbar squared and divide by wave function
  e_kinetic = 0.5*h2*e_kinetic/wfold;

  // compute the potential energy
  e_potential = 0;
  // contribution from electron-proton potential
  for (i = 0; i < number_particles; i++) {
    r_single_particle = 0;
    for (j = 0; j < dimension; j++) {
      r_single_particle += r(i,j)*r(i,j);
    }
//    e_potential -= charge/sqrt(r_single_particle);
    e_potential += 0.5*r_single_particle; // TODO: OMEGA IS MISSING HERE! 
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

/* -------------------------------------------------------------------------- *
 *                    Random number generator                                 *
 * -------------------------------------------------------------------------- */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
   int             j;
   long            k;
   static long     iy=0;
   static long     iv[NTAB];
   double          temp;

   if (*idum <= 0 || !iy) {
      if (-(*idum) < 1) *idum=1;
      else *idum = -(*idum);
      for(j = NTAB + 7; j >= 0; j--) {
         k     = (*idum)/IQ;
         *idum = IA*(*idum - k*IQ) - IR*k;
         if(*idum < 0) *idum += IM;
         if(j < NTAB) iv[j] = *idum;
      }
      iy = iv[0];
   }
   k     = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   j     = iy/NDIV;
   iy    = iv[j];
   iv[j] = *idum;
   if((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}

