/*
 * MPIBasicRunner<T>.cpp
 *
 *  Created on: Aug 24, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
*/

#include "dmc.h"

namespace runners {

template <class T>
void MPIBasicRunner<T>::initialize()
{

}

template <class T>
void MPIBasicRunner<T>::masterRun()
{
	for(int m=0;m<this->runParams.bins;m++)
	{
		fprintf(this->log,"Bin: %d\n",m);
		fflush(this->log);

		//Global Gather
		for(int o=0;o<this->MPIobservablesCollection.size();o++)
			this->MPIobservablesCollection[o]->parallelReceive();

		fprintf(this->log,"Parallel gather done.\n");
		fprintf(this->log,"Beginning Write of %d observables.\n",this->observablesCollection.size());
		fflush(this->log);

		//Write to file
		for(int o=0;o<this->observablesCollection.size();o++)
		{
//			fprintf(this->log,"Writing observable %d.\n",o);
//			fflush(this->log);
			this->observablesCollection[o]->writeViaIndex(m);
		}

		fprintf(this->log,"Write Done.\n");
		fflush(this->log);
	}
}

template <class T>
void MPIBasicRunner<T>::run()
{
//	double gaccept = 0.0;
//	double divisor = 1.0/this->runParams.nsteps/this->walkers.walkerCount/walkers[0]->state.particleCount;

	for(int m=0;m<this->runParams.bins;m++)
	{
		walkers.resetWalkers();

		fprintf(this->log,"Starting Bin: %d with %d walkers\n",m,walkers.walkerCount);
		fflush(this->log);
		//Initalize walkers
		for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
			mover->initialize(it->second);

		//Propagate in time
		for(int i=0;i<this->runParams.nsteps;i++)
		{
#if 0
			fprintf(this->log,"Step: %d\n",i);
			fflush(this->log);
#endif
			for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				mover->move(it->second);

			//Now do measure for each walker
			for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				it->second->measure();
#ifndef NOBRANCH
			//Compact
			branch();
#endif
		}

		fprintf(this->log,"Performing local gather\n");
		fflush(this->log);

		//Local Gather
		//Trace over all observables
		for(int o=0;o<this->observablesCollection.size();o++)
		{
			//Trace over all walkers and accumulate required observable
			for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				it->second->observablesCollection[o]->gather((void *) this->observablesCollection[o]);
		}

		fprintf(this->log,"Local gather done\n");
		fflush(this->log);

		//Global Gather
		for(int o=0;o<this->MPIobservablesCollection.size();o++)
			this->MPIobservablesCollection[o]->parallelSend();

		fprintf(this->log,"Parallel gather done\n");
		fflush(this->log);

		fprintf(this->log,"Ending Bin: %d with %d walkers\n",m,walkers.walkerCount);
		fflush(this->log);
	}
}

template <class T>
void MPIBasicRunner<T>::finalize()
{
	for(int o=0;o<this->observablesCollection.size();o++)
		this->observablesCollection[o]->clear();
}

template <class T>
void MPIBasicRunner<T>::masterFinalize()
{
	//Notify slaves to finish wait
	int tag, statusRun = MPISTATUSFINISH;
	for(int i=1;i<this->procCount;i++)
	{
		MPI_Send(&statusRun,1,MPI_INT,i,tag,MPI_COMM_WORLD);
	}
}

template <class T>
void MPIBasicRunner<T>::branch()
{
	vector<int> ridx;

	int count = walkers.walkerCount, i = 0;
	for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end() && i<walkers.walkerCount;++it)
	{
		if(it->second->state.weight.logValue() >= walkers.maxValue && count<walkers.maxWalkerCount)
		{
			it->second->state.weight.update(0.5);

			Walker<T>* wcpy = it->second->duplicate();
			(*walkers.walkerCollection)[walkers.lastIndex++] = wcpy;
			wcpy->state.weight.update(0.5);

			i++;
			count++;
		}
		else if(it->second->state.weight.logValue() <= walkers.minValue)
			ridx.push_back(it->first);
	}

	for(int t=0;t<ridx.size();t+=2)
	{
		double w1 = walkers[ridx[t]]->state.weight.logValue();
		double w2 = 0.0;
		if(t+1<ridx.size())
			w2 = walkers[ridx[t+1]]->state.weight.logValue();

		if(gsl_rng_uniform(walkers[ridx[t]]->rgenref)<w1/(w1+w2))
			walkers.walkerCollection->erase(ridx[t]);
	}

	walkers.walkerCount = walkers.walkerCollection->size();
}

///////////////////////////////////////////////////////////////////////////////
template void MPIBasicRunner<int>::initialize();
template void MPIBasicRunner<int>::masterRun();
template void MPIBasicRunner<int>::run();
template void MPIBasicRunner<int>::finalize();
template void MPIBasicRunner<int>::masterFinalize();
template void MPIBasicRunner<int>::branch();

} /* namespace runners */
