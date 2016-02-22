/*
 * MPIBasicRunner<T,U>.cpp
 *
 *  Created on: Aug 24, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
*/

#include "dmc.h"

namespace runners {

template <class T, class U>
void MPIBasicRunner<T,U>::displayBranchStat(int nbranches)
{
	fprintf(this->log,"\n***************************************\n");
	fprintf(this->log,"Total branches: %d\n",nbranches);
	fprintf(this->log,"No. of splits: %d\n",this->nclones);
	fprintf(this->log,"No. of eliminates: %d\n",this->nelims);
}

template <class T, class U>
void MPIBasicRunner<T,U>::initialize()
{

}

template <class T, class U>
void MPIBasicRunner<T,U>::masterRun()
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

#if !defined NOBRANCH && !defined LOCALBRANCH
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
#endif

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
void MPIBasicRunner<T,U>::run()
{
//	double gaccept = 0.0;
//	double divisor = 1.0/this->runParams.nsteps/this->walkers.walkerCount/walkers[0]->state.particleCount;

//	int eqbranch = this->runParams.eSteps*0.90*BRANCHPERCENT;
//	eqbranch = (eqbranch<1) ? 1 : eqbranch;
	int databranch = runParams.nsteps*BRANCHPERCENT;
	databranch = (databranch<1) ? 1 : databranch;

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
#if 0
		//Spend 90% of time getting to weighted distribution
		for(int i=0;i<this->runParams.eSteps*0.90;i++)
		{
			for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				mover->move(it->second);

			if((i+1)%eqbranch == 0)
			{
				nbranches++;
				walkers.displayWalkers(this->log);
				branch();
				walkers.displayWalkers(this->log);
			}
		}
#endif

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
			{	
//				fprintf(this->log,"Step: %d\n",i); 
				it->second->measure();
			}

#ifndef NOBRANCH
			//branch
			if((i+1)%databranch == 0)
			{
				nbranches++;
				//shortWalkerDisplay();
				branch(i);
				//shortWalkerDisplay();
			}
#endif
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

#if !defined NOBRANCH && !defined LOCALBRANCH
		//Master needs to be notified that process has finished a bin
		int tag, smsg = MPIBINDONE, rmsg = 0;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(&smsg,1,MPI_CHAR,&rmsg,1,MPI_CHAR,0,MPI_COMM_WORLD);
#endif

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

template <class T, class U>
void MPIBasicRunner<T,U>::finalize()
{
	for(int o=0;o<this->observablesCollection.size();o++)
		this->observablesCollection[o]->clear();
}

template <class T, class U>
void MPIBasicRunner<T,U>::masterFinalize()
{
	//Notify slaves to finish wait
	int tag, statusRun = MPISTATUSFINISH;
	for(int i=1;i<this->procCount;i++)
	{
		MPI_Send(&statusRun,1,MPI_INT,i,tag,MPI_COMM_WORLD);
	}
}

#if !defined CONSTPOPBRANCH

#else
template <class T, class U>
void MPIBasicRunner<T,U>::branch(int step)
{
#if !defined LOCALBRANCH
	////////////////////////////////////////////////////////////
	//Tell Master to Start expecting branch requests
	//Keep this header no matter what. It pairs with the loop control
	//dispatch of masterRun.
	////////////////////////////////////////////////////////////
	fprintf(this->log,"\nBranching started: %d\n",branchcount++);
	char smsg = MPIBRANCH, rmsg = 0;
	int tag = 0;
	MPI_Status stat;

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(&smsg,1,MPI_CHAR,&rmsg,1,MPI_CHAR,0,MPI_COMM_WORLD);
//#if DEBUG >= 3
	fprintf(this->log,"Notified master to start expecting data.\n");
	fflush(this->log);
//#endif
	////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////
	//Only one process needs to tell master what the current sweep is
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	fprintf(this->log,"Rank %d.\n",rank);
	fflush(this->log);
	if(rank == 1)
		MPI_Send(&step,1,MPI_INT,0,tag,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////
	//Can generalize to a sweep schedule but fixed for now
	//Not a good idea to branch while Q is small wait
	//for some Q accummulation before pruning.
	////////////////////////////////////////////////////////////
	if((float) step/this->runParams.nsteps > 0.10)
		branchLimited();
#else
	//Local branching
	if((float) step/this->runParams.nsteps > 0.10)
		branchLocal(0.25); // Do not prune more than specified percent
	else if((float) step/this->runParams.nsteps > 0.60)
		branchLocal(0.50);
	else if((float) step/this->runParams.nsteps > 0.90)
		branchLocal(0.90);
#endif
}
#endif

template <class T, class U>
void MPIBasicRunner<T,U>::masterBranch()
{
#if !defined LOCALBRANCH
	////////////////////////////////////////////////////////////
	//Get sweep number
	MPI_Status stat;
	int step, tag;
	MPI_Recv(&step,1,MPI_INT,1,tag,MPI_COMM_WORLD,&stat);
	MPI_Barrier(MPI_COMM_WORLD);
	////////////////////////////////////////////////////////////
	fprintf(this->log,"\nBranching started: %d\n",branchcount++);
	fflush(this->log);

	////////////////////////////////////////////////////////////
	//Can generalize to a sweep schedule but fixed for now
	//Not a good idea to branch while Q is small wait
	//for some Q accummulation before pruning.
	////////////////////////////////////////////////////////////
	if((float) step/this->runParams.nsteps > 0.10)
		this->masterBranchLimited(1.0); // Do not prune more than specified percent

#else

#endif
}

///////////////////////////////////////////////////////////////////////////////
template void MPIBasicRunner<int,stringstream>::initialize();
template void MPIBasicRunner<int,stringstream>::masterRun();
template void MPIBasicRunner<int,stringstream>::run();
template void MPIBasicRunner<int,stringstream>::finalize();
template void MPIBasicRunner<int,stringstream>::masterFinalize();
template void MPIBasicRunner<int,stringstream>::branch(int);
template void MPIBasicRunner<int,stringstream>::masterBranch();
template void MPIBasicRunner<int,stringstream>::displayBranchStat(int);

} /* namespace runners */
