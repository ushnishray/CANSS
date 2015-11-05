/*
 * MPIBasicRunner.h
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

#ifndef MPIBASICRUNNER_H_
#define MPIBASICRUNNER_H_

namespace runners {

struct MPIBRParams
{
	int bins;
	int nsteps;

	MPIBRParams(int _bins, int _nsteps):bins(_bins),nsteps(_nsteps) {}
};

template <class T>
struct Walkers {
	vector<core::Walker<T>*>* walkerCollection;

	double maxValue;
	double minValue;

	int nextid;
	int walkerCount;

	Walkers()
	{
		maxValue = minValue = 0.0;
		nextid = walkerCount = 0;
		walkerCollection = NULL;
	}

	Walkers(vector<core::Walker<T>*>* _wc, double maxv, double minv):walkerCollection(_wc),maxValue(maxv),minValue(minv)
	{
		nextid = 0;
		walkerCount = walkerCollection->size();
	}

	core::Walker<T>* operator[](std::size_t idx) { return (*walkerCollection)[idx];}
};

template <class T>
class MPIBasicRunner {

private:

	//Mover
	Mover<T>* mover;

	//Run parameters
	MPIBRParams& runParams;

	//For slaves
	Walkers<T>& walkers;

	//For master
	int procCount;

	//For all
	vector<measures::Observable<T>*>& observablesCollection;
	vector<measures::MPIObservable*>& MPIobservablesCollection;

	//For output
	FILE* log;

public:

	//For master
	MPIBasicRunner(FILE* _log, int _pcount,
			int _bins, int _nsteps,
			vector<measures::Observable<T>*>& _oc, vector<measures::MPIObservable*>& _moc
			):
		log(_log),procCount(_pcount),mover(NULL),
		runParams(*(new MPIBRParams(_bins,_nsteps))),
		walkers(*(new Walkers<T>)),
		observablesCollection(_oc),MPIobservablesCollection(_moc)
	{ }

	//For slaves
	MPIBasicRunner(FILE* _log, int _pcount, Mover<T>* _mover,
			int _bins, int _nsteps,
			double maxv,double minv,vector<core::Walker<T>*>* _wc,
			vector<measures::Observable<T>*>& _oc, vector<measures::MPIObservable*>& _moc
			):
		log(_log),procCount(_pcount),mover(_mover),
		runParams(*(new MPIBRParams(_bins,_nsteps))),
		walkers(*(new Walkers<T>(_wc,maxv,minv))),
		observablesCollection(_oc),MPIobservablesCollection(_moc)
	{ }

	~MPIBasicRunner()
	{
		delete &runParams;
		delete &walkers;
	}

	//For slaves
	virtual void initialize();
	void run();
	void finalize();

	//For master
	void masterRun();
	void masterFinalize();

	//Brancing
	void branch();
};

} /* namespace runners */
#endif /* MPIBASICRUNNER_H_ */
