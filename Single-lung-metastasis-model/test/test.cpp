#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <random>

#include "../BioFVM.h"

using namespace BioFVM;

int omp_num_threads = 4; // set number of threads for parallel computing
// set this to # of CPU cores x 2 (for hyperthreading)
double pi= 3.1415926535897932384626433832795;

int main( int argc, char* argv[] )
{
    // openmp setup
    omp_set_num_threads(omp_num_threads);


    return 0;
}
