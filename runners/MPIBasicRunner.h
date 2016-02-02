/*
 * MPIBasicRunner.h
 *
 *  Created on: Aug 24, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
 */

#ifndef MPIBASICRUNNER_H_
#define MPIBASICRUNNER_H_

namespace runners {

struct MPIBRParams
{
	int eSteps;

	int bins;
	int nsteps;

	MPIBRParams(int _eSteps, int _bins, int _nsteps):eSteps(_eSteps),bins(_bins),nsteps(_nsteps) {}
};

template <class T>
struct Walkers {
	//vector<core::Walker<T>*>* walkerCollection;
	NumMap<Walker<T>>* walkerCollection;

	double maxValue;
	double minValue;

	int walkerCount;
	int maxWalkerCount;

	/////////////////////////////
	int initWalkerCount;
	long lastIndex;

	Walkers()
	{
		maxValue = minValue = 0.0;
		lastIndex = initWalkerCount = walkerCount = maxWalkerCount = 0;
		walkerCollection = NULL;
	}

	Walkers(NumMap<Walker<T>>* _wc, int _maxWalkerCount, double maxv, double minv):walkerCollection(_wc),maxValue(maxv),minValue(minv)
	{
		maxWalkerCount = _maxWalkerCount;
		lastIndex = initWalkerCount = walkerCount = walkerCollection->size();
	}

	void resetWalkers()
	{
#ifndef REINDEX
		//Not re-indexing: Faster but for very large walker populations lastIndex will be an issue
		//Need to do this if branching happened
		if(walkerCount == initWalkerCount)
			return;

		if(walkerCount > initWalkerCount)
		{
			int sizeToTrim = walkerCount - initWalkerCount;
			int i = 0;
			for(typename NumMap<Walker<T>>::iterator it = walkerCollection->begin();it!=walkerCollection->end() && i<sizeToTrim;++it)
			{
				walkerCollection->erase(it);
				i++;
			}
		}
		else
		{
			int sizeToAdd = initWalkerCount - walkerCount;
			Walker<T>* wadd = walkerCollection->begin()->second->duplicate();

			for(int i=0;i<sizeToAdd;i++)
				(*walkerCollection)[lastIndex++] = wadd->duplicate();

			delete wadd;
		}
#else
		//Reindex - slower but is suited for possibly very large number of walkers
		Walker<T>* wadd = walkerCollection->begin()->second->duplicate();
		lastIndex = 0;
		walkerCollection->clear();
		for(int i=0;i<initWalkerCount;i++)
			(*walkerCollection)[i] = wadd->duplicate();
		delete wadd;
#endif
		walkerCount = initWalkerCount;
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
			int _esteps, int _bins, int _nsteps,
			vector<measures::Observable<T>*>& _oc, vector<measures::MPIObservable*>& _moc
			):
		log(_log),procCount(_pcount),mover(NULL),
		runParams(*(new MPIBRParams(_esteps,_bins,_nsteps))),
		walkers(*(new Walkers<T>)),
		observablesCollection(_oc),MPIobservablesCollection(_moc)
	{ }

	//For slaves
	MPIBasicRunner(FILE* _log, int _pcount, Mover<T>* _mover,
			int _esteps, int _bins, int _nsteps,
			int maxwc, double maxv,double minv,NumMap<Walker<T>>* _wc,
			vector<measures::Observable<T>*>& _oc, vector<measures::MPIObservable*>& _moc
			):
		log(_log),procCount(_pcount),mover(_mover),
		runParams(*(new MPIBRParams(_esteps,_bins,_nsteps))),
		walkers(*(new Walkers<T>(_wc,maxwc,maxv,minv))),
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
