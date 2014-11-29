#include <iostream>
#include <cmath>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
using std::cout;
using std::endl;
int main(int args, char *argv[]) {
  double a = 0;
  int numThreads = atoi(argv[1]);
#pragma omp parallel reduction(+:a) num_threads(numThreads)
  cout << "Hello world from thread number " << omp_get_thread_num() << endl;
#pragma omp for
for(unsigned int i=0; i<100000000; i++) {
    a += pow(1.0/sqrt(i),1.0/3);
  }

  cout << a << endl;
}
