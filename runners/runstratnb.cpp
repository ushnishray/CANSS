/*
 * runstratnb.cpp
 *
 *  Created on: Mar 3, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */
#include "dmc.h"

template <class T, class U>
void MPIBasicRunner<T,U>::masterRunNB()
{
	char* rmsg = new char[this->procCount];
	char smsg = 0, msg = 0;
	for(int i = 0;i<procCount;i++)
		rmsg[i] = 0;
	int tag;
	MPI_Status stat;

	for(int m=0;m<this->runParams.bins;m++)
	{
		fprintf(this->log,"Bin: %d\n",m);
		fflush(this->log);

		for(int i=0;i<this->runParams.nsteps;i++)
		{
			//Global Gather
			for(int o=0;o<this->MPIobservablesCollection.size();o++)
				this->MPIobservablesCollection[o]->parallelReceive();

#if DEBUG>=4
			fprintf(this->log,"Parallel gather done.\n");
			fprintf(this->log,"Beginning Write of %d observables.\n",this->observablesCollection.size());
			fflush(this->log);
#endif
			MPI_Barrier(MPI_COMM_WORLD);
		}

		//Write to file
		for(int o=0;o<this->observablesCollection.size();o++)
			this->observablesCollection[o]->writeViaIndex(m);

		fprintf(this->log,"Write Done.\n");
		fflush(this->log);

		//Synchronize
		MPI_Barrier(MPI_COMM_WORLD);
	}

	delete[] rmsg;
}

template <class T, class U>
void MPIBasicRunner<T,U>::runNB()
{
	for(int m=0;m<this->runParams.bins;m++)
	{
		int nbranches = 0;
		walkers.resetWalkers();

		fprintf(this->log,"Starting Bin: %d with %d walkers\n",m,walkers.walkerCount);
		fflush(this->log);
		//Initalize walkers
		for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
			mover->initialize(it->second);

#if 0
		//Display Init State
		for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
			it->second->display();
#endif

		//Do equilibration which is about reaching stead state
		for(int i=0;i<this->runParams.eSteps;i++)
		{
			for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				mover->move(it->second);
		}

		//Reset walker times etc. equilibration is just about eliminating crazy transients.
		for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
			it->second->reset();

		//Propagate in time
		for(int i=0;i<this->runParams.nsteps;i++)
		{
			for(int b=0;b<this->runParams.branchStep;b++)
			{
				for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
					mover->move(it->second);

				//Now do measure for each walker
				for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
					it->second->measure();
			}

			//Local Gather
			//Trace over all observables
			for(int o=0;o<this->observablesCollection.size();o++)
			{
				//Trace over all walkers and accumulate required observable
				for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
					it->second->observablesCollection[o]->gather((void *) this->observablesCollection[o]);
			}

			//Global Gather
			for(int o=0;o<this->MPIobservablesCollection.size();o++)
				this->MPIobservablesCollection[o]->parallelSend();

#if DEBUG>=4
			fprintf(this->log,"Parallel gather done\n");
			fflush(this->log);
#endif

			//Reset walker state.
			for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				it->second->reset();

			//Synchronize
			MPI_Barrier(MPI_COMM_WORLD);
		}

		fprintf(this->log,"Ending Bin: %d with %d walkers\n",m,walkers.walkerCount);
		fflush(this->log);

		//Synchronize
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
template void MPIBasicRunner<int,stringstream>::masterRunNB();
template void MPIBasicRunner<int,stringstream>::runNB();

template void MPIBasicRunner<float,stringstream>::masterRunNB();
template void MPIBasicRunner<float,stringstream>::runNB();


