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
			
		//Do equilibration
		for(int i=0;i<this->runParams.eSteps;i++)
		{
			for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				mover->move(it->second);
		}	

		//Reset walker 
		for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
			it->second->reset();	

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
			{	
//				fprintf(this->log,"Step: %d\n",i); 
				it->second->measure();
			}
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
	paircomp paircompobj;
	vector<pair<int,double>> widx;

	int scount = 0, ccount = 0;
	for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
	{
		widx.push_back(pair<int,double>(it->first,it->second->state.weight.logValue()));

		if(it->second->state.weight.logValue() >= walkers.maxValue)
			scount++;
		else if(it->second->state.weight.logValue() <= walkers.minValue)
			ccount++;
	}

	//Merging involves 2 files so round up to an even number
	ccount = int(ccount/2+0.5)*2;

	//walkers to keep untouched
	int kcount = walkers.walkerCount - scount - ccount;
	int dpop = kcount + 2*scount + int(ccount/2+0.5) - walkers.walkerCount;

	sort(widx.begin(),widx.end(),paircompobj);

	//if dpop>0 then there are more splits than merges so we need to remove low weights
	//if dpop<0 then there are more merges than splits so we need to split large weights
	if(dpop>0)
		ccount = 2*scount;
	else
		scount = ccount/2;

	//Now merge and split. For every 2 merged we can do 1 split
	int p = 0;
	for(int i = widx.size()-ccount;i<widx.size();i+=2)
	{
		walkers[widx[p].first]->state.weight.update(0.5);
		Walker<T>* wcpy = walkers[widx[p].first]->duplicate();

		double w1 = walkers[widx[i].first]->state.weight.logValue();
		double w2 = walkers[widx[i+1].first]->state.weight.logValue();
		int idx = 0;
		if(gsl_rng_uniform(walkers[i]->rgenref)<w1/(w1+w2))
			idx = widx[i].first;
		else
			idx = widx[i+1].first;

		walkers.walkerCollection->erase(idx);
		(*walkers.walkerCollection)[idx] = wcpy;
		p++;
	}
}

///////////////////////////////////////////////////////////////////////////////
template void MPIBasicRunner<int>::initialize();
template void MPIBasicRunner<int>::masterRun();
template void MPIBasicRunner<int>::run();
template void MPIBasicRunner<int>::finalize();
template void MPIBasicRunner<int>::masterFinalize();
template void MPIBasicRunner<int>::branch();

} /* namespace runners */
