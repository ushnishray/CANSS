/*
 * constants.h
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

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#define DIMENSION 1

#define DIVISOR 1e5
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

//debug flags
//#define DEBUG

#endif /* CONSTANTS_H_ */
