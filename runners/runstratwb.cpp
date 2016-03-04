/*
 * runstratwb.cpp
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
void MPIBasicRunner<T,U>::masterRunWB()
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

		//Expect branching
		do{
#if DEBUG >= 3
			fprintf(this->log,"\nIn Control Loop\n");
#endif
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Gather(&smsg,1,MPI_CHAR,rmsg,1,MPI_CHAR,0,MPI_COMM_WORLD);
			msg = rmsg[1];
			for(int i=2;i<procCount;i++)
			{
				msg &= rmsg[i];
#if DEBUG >= 3
				fprintf(this->log,"Received %d from %d\n",rmsg[i],i);
#endif
			}
#if DEBUG >= 3
			fflush(this->log);
#endif

			if(msg==MPIBRANCH)
				masterBranch();
		}while(msg != MPIBINDONE);

		//Global Gather
		for(int o=0;o<this->MPIobservablesCollection.size();o++)
			this->MPIobservablesCollection[o]->parallelReceive();

		fprintf(this->log,"Parallel gather done.\n");
		fprintf(this->log,"Beginning Write of %d observables.\n",this->observablesCollection.size());
		fflush(this->log);

		//Write to file
		for(int o=0;o<this->observablesCollection.size();o++)
		{
			//Ugly but easy
			if( BasicObs<T,U>* obj = dynamic_cast<BasicObs<T,U>*>(this->observablesCollection[o]))
				obj->freeEnergy.copy(FreeEnergy);

			this->observablesCollection[o]->writeViaIndex(m);
		}

		FreeEnergy.resetValue(); //IMPORTANT
		fprintf(this->log,"Write Done.\n");
		fflush(this->log);

		//Synchronize
		MPI_Barrier(MPI_COMM_WORLD);
	}

	delete[] rmsg;
}

template <class T, class U>
void MPIBasicRunner<T,U>::runWB()
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
#if 0
			fprintf(this->log,"Step: %d\n",i);
			fflush(this->log);
#endif

			for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				mover->move(it->second);

			//Now do measure for each walker
			for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				it->second->measure();

			//branch
			if((i+1)%runParams.branchStep == 0)
			{
				nbranches++;
				branch(i);
			}
		}

		this->displayBranchStat(nbranches);
		this->nclones = this->nelims = 0;
		fprintf(this->log,"\nPerforming local gather\n");
		fflush(this->log);

		//Local Gather
		//Trace over all observables
		for(int o=0;o<this->observablesCollection.size();o++)
		{
			//Trace over all walkers and accumulate required observable
			for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				it->second->observablesCollection[o]->gather((void *) this->observablesCollection[o]);
		}

		fprintf(this->log,"Local gather done\n");
		fflush(this->log);

		//Master needs to be notified that process has finished a bin
		int tag, smsg = MPIBINDONE, rmsg = 0;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(&smsg,1,MPI_CHAR,&rmsg,1,MPI_CHAR,0,MPI_COMM_WORLD);

		//Global Gather
		for(int o=0;o<this->MPIobservablesCollection.size();o++)
			this->MPIobservablesCollection[o]->parallelSend();

		fprintf(this->log,"Parallel gather done\n");
		fflush(this->log);

		fprintf(this->log,"Ending Bin: %d with %d walkers\n",m,walkers.walkerCount);
		fflush(this->log);

		//Synchronize
		MPI_Barrier(MPI_COMM_WORLD);
	}
}


//////////////////////////////////////////////////////////////////////////////
template void MPIBasicRunner<int,stringstream>::masterRunWB();
template void MPIBasicRunner<int,stringstream>::runWB();

template void MPIBasicRunner<float,stringstream>::masterRunWB();
template void MPIBasicRunner<float,stringstream>::runWB();
