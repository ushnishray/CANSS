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
#ifdef NOBRANCH
	masterRunNB();
#else
	masterRunWB();
#endif
}

template <class T, class U>
void MPIBasicRunner<T,U>::run()
{
#ifdef NOBRANCH
	runNB();
#else
	runWB();
#endif
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
	//if((float) step/this->runParams.nsteps > 0.10)
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
	//if((float) step/this->runParams.nsteps > 0.10)
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

template void MPIBasicRunner<float,stringstream>::initialize();
template void MPIBasicRunner<float,stringstream>::masterRun();
template void MPIBasicRunner<float,stringstream>::run();
template void MPIBasicRunner<float,stringstream>::finalize();
template void MPIBasicRunner<float,stringstream>::masterFinalize();
template void MPIBasicRunner<float,stringstream>::branch(int);
template void MPIBasicRunner<float,stringstream>::masterBranch();
template void MPIBasicRunner<float,stringstream>::displayBranchStat(int);
} /* namespace runners */
