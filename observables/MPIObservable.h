/*
 * MPIObservable.cpp
 *
 *  Created on: Aug 24, 2014
 *      Author: ushnish
 *
 	Copyright (c) 2014 Ushnish Ray
	All rights reserved.
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
