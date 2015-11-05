/*
 * mainIncs.h
 *
 *  Created on: Aug 12, 2014
 *      Author: ushnish
 *
 	Copyright (c) 2014 Ushnish Ray
	All rights reserved.

	This program and the accompanying materials are made available under explicit agreement
	between Ushnish Ray and the end user. You may not redistribute the code and
	accompanying material to anyone.

	On the event that the software is used to generate data that is used implicitly or explicitly
	for research purposes, proper acknowledgment must provided  in the citations section of
	publications.

	This software is cannot be used for commercial purposes in any way whatsoever.
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

//Advanced Data Structures
#include <set>
#include <vector>
#include <unordered_map>
#include <map>

//Random number generator
#include <gsl/gsl_rng.h>
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
