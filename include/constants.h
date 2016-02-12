/*
 * constants.h
 *
 *  Created on: Aug 12, 2014
 *      Author: ushnish
 *
 	Copyright (c) 2014 Ushnish Ray
	All rights reserved.
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#define DIMENSION 1

#define ZEROTOL 1.0e-9

typedef double Real;
typedef int Integer;

//Float Precsions and Widths
#define FIELDPRECISION 9
#define FIELDWIDTH 17
#define FIELDFORMAT std::ios_base::scientific


//Error Codes
#define FILENOTFOUND -1
#define NOTALLOWED -2
#define SUCCESS 1
#define FAIL 0
#define DIMERROR -3

//MPI Status Message
#define MPISTATUSFINISH 1
#define MPISTATUSSTART 2

//MPI Message
#define MPIBRANCH 100
#define MPIBINDONE 101

//debug flags
//#define DEBUG

#endif /* CONSTANTS_H_ */
