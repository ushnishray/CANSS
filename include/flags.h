/*
 * flags.h
 *
 *  Created on: Aug 28, 2015
 *      Author: ushnish
 *
 	Copyright (c) 2014 Ushnish Ray
	All rights reserved.
 */

#ifndef FLAGS_H_
#define FLAGS_H_

///////////////////////////////////////////////////////
//branching options
///////////////////////////////////////////////////////

//MINIMUM PROBABILITY TO KEEP log(1.0e-5)
#define MINBRANCHWEIGHT -11.51292546
//Branching Every (%) states
#define BRANCHPERCENT 0.05
//Constant Population Branching
#define CONSTPOPBRANCH
//#define CPB1
#define CPB2

//Turnoff Branching - on by default
//#define NOBRANCH

//Reindex - off by default
//#define REINDEX

///////////////////////////////////////////////////////
//Profiler
///////////////////////////////////////////////////////

#define PROFILER

//Debug
//#define DEBUG

#endif /* FLAGS_H_ */
