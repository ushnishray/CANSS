/*
 * MPIObservable.cpp
 *
 *  Created on: Aug 24, 2014
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

#ifndef MPIOBSERVABLE_H_
#define MPIOBSERVABLE_H_

namespace measures
{

class MPIObservable
{
protected:

	int processId;
	int procCount;
	bool MPIEnabled;

public:

	MPIObservable()
	{
		processId = 0;
		procCount = 0;
		MPIEnabled = false;
	}

	MPIObservable(int pid,int nprocs):processId(pid),procCount(nprocs) {MPIEnabled = true;}

	virtual int parallelSend() = 0; //To be called by slaves
	virtual int parallelReceive() = 0; //To be called by master
};

}
#endif /* MPIOBSERVABLE_H_ */
