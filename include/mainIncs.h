/*
 * mainIncs.h
 *
 *  Created on: Aug 12, 2014
 *      Author: ushnish
 *
 	Copyright (c) 2014 Ushnish Ray
	All rights reserved.
 */

#ifndef MAININCS_H_
#define MAININCS_H_

//Basic headers
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>

//Libraries
#include <cmath>
#include <ctime>
#include <string>
#include <cstring>
#include <sstream>
#include <cassert>
#include <exception>

//Advanced Data Structures
#include <algorithm>
#include <set>
#include <vector>
#include <unordered_map>
#include <map>

//Random number generator
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//Optimization library
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

//Parallel Libraries
#include "omp.h"
#include "mpi.h"

//Helper Routines of importance
#include "constants.h"

using namespace std;
using namespace __gnu_cxx;

#endif /* MAININCS_H_ */
