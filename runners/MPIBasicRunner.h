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

template <class T, class U>
struct Walkers {
	//vector<core::Walker<T,U>*>* walkerCollection;
	NumMap<Walker<T,U>>* walkerCollection;

	double maxValue;
	double minValue;

	int walkerCount;
	int maxWalkerCount;

	/////////////////////////////
	int initWalkerCount;
	long lastIndex;

	Walkers(int _walkerCount)
	{
		maxValue = minValue = 0.0;
		lastIndex = 0;
		initWalkerCount = walkerCount = maxWalkerCount = _walkerCount;
		walkerCollection = NULL;
	}

	Walkers(NumMap<Walker<T,U>>* _wc, int _maxWalkerCount, double maxv, double minv):walkerCollection(_wc),maxValue(maxv),minValue(minv)
	{
		maxWalkerCount = _maxWalkerCount;
		lastIndex = initWalkerCount = walkerCount = walkerCollection->size();
	}

	void displayWalkers(FILE* out)
	{
		fprintf(out,"********************************************************************\n");
		for(typename NumMap<Walker<T,U>>::iterator it = walkerCollection->begin();it!=walkerCollection->end();++it)
		{
			fprintf(out,"Walker id: %d\n",it->first);
			it->second->display();
		}
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
			for(typename NumMap<Walker<T,U>>::iterator it = walkerCollection->begin();it!=walkerCollection->end() && i<sizeToTrim;++it)
			{
				walkerCollection->erase(it);
				i++;
			}
		}
		else
		{
			int sizeToAdd = initWalkerCount - walkerCount;
			Walker<T,U>* wadd = walkerCollection->begin()->second->duplicate();

			for(int i=0;i<sizeToAdd;i++)
				(*walkerCollection)[lastIndex++] = wadd->duplicate();

			delete wadd;
		}
#else
		//Reindex - slower but is suited for possibly very large number of walkers
		Walker<T,U>* wadd = walkerCollection->begin()->second->duplicate();
		lastIndex = 0;
		walkerCollection->clear();
		for(int i=0;i<initWalkerCount;i++)
			(*walkerCollection)[i] = wadd->duplicate();
		delete wadd;
#endif
		walkerCount = initWalkerCount;
	}

	core::Walker<T,U>* operator[](std::size_t idx) { return (*walkerCollection)[idx];}
};

template <class T, class U>
class MPIBasicRunner {

private:

	//Mover
	Mover<T,U>* mover;

	//Run parameters
	MPIBRParams& runParams;

	//For slaves
	Walkers<T,U>& walkers;

	//For master
	int procCount;

	//For all
	vector<measures::Observable<T,U>*>& observablesCollection;
	vector<measures::MPIObservable*>& MPIobservablesCollection;

	//For output
	FILE* log;

	//Local statistical variables
	unsigned int branchcount;
	long nclones;
	long nelims;

	//Display parameters
	void displayBranchStat(int);

public:

	//For master
	MPIBasicRunner(FILE* _log, int _pcount,
			int _esteps, int _bins, int _nsteps,int _wcount,
			vector<measures::Observable<T,U>*>& _oc, vector<measures::MPIObservable*>& _moc
			):
		log(_log),procCount(_pcount),mover(NULL),
		runParams(*(new MPIBRParams(_esteps,_bins,_nsteps))),
		walkers(*(new Walkers<T,U>(_wcount))),
		observablesCollection(_oc),MPIobservablesCollection(_moc)
	{
		branchcount = 0;
		nclones = 0;
		nelims = 0;
	}

	//For slaves
	MPIBasicRunner(FILE* _log, int _pcount, Mover<T,U>* _mover,
			int _esteps, int _bins, int _nsteps,
			int maxwc, double maxv,double minv,NumMap<Walker<T,U>>* _wc,
			vector<measures::Observable<T,U>*>& _oc, vector<measures::MPIObservable*>& _moc
			):
		log(_log),procCount(_pcount),mover(_mover),
		runParams(*(new MPIBRParams(_esteps,_bins,_nsteps))),
		walkers(*(new Walkers<T,U>(_wc,maxwc,maxv,minv))),
		observablesCollection(_oc),MPIobservablesCollection(_moc)
	{
		branchcount = 0;
		nclones = 0;
		nelims = 0;
	}

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

	//Branching
	void branch(int);
	void masterBranch();

	//branching algorithms
	void branchLimited();
	void masterBranchLimited(float);
	void branchFull();
	void masterBranchFull();

	//Helpers
	void shortWalkerDisplay()
	{
		//Report currents
		fprintf(this->log,"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
		for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
//			fprintf(this->log,"Walker %d Current %10.6e, %10.6e, %10.6e\n",it->first,it->second->state.Q.x.logValue(),it->second->state.Q.y.logValue(),it->second->state.Q.z.logValue());
			fprintf(this->log,"Walker %d Current %10.6e\n",it->first,it->second->state.weight.logValue());
		fprintf(this->log,"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
		fflush(this->log);
	}
};

} /* namespace runners */
#endif /* MPIBASICRUNNER_H_ */
