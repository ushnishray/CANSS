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
//Branching Every (%) states or specs
#define BRANCHPERCENT (-1.0)
//Constant Population Branching
#define CONSTPOPBRANCH

///////////////////////////////////////////////////////
//Debug output options
///////////////////////////////////////////////////////
//Can produce large files so watch out!
//Level 2 is Observables
//Level 3 is BRANCHING
#define DEBUG 1

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
