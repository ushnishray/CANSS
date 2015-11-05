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
		fprintf(this->log,"Bin: %d\n",m);
		fflush(this->log);

		//Initalize walkers
		for(int t=0;t<this->walkers.walkerCount;t++)
			mover->initialize(this->walkers[t]);

		//Propagate in time
		for(int i=0;i<this->runParams.nsteps;i++)
		{
#if 0
			fprintf(this->log,"Step: %d\n",i);
			fflush(this->log);
#endif
			for(int t=0;t<this->walkers.walkerCount;t++)
				mover->move(this->walkers[t]);

			//Now do measure for each walker
			for(int t=0;t<this->walkers.walkerCount;t++)
				this->walkers[t]->measure();

			//Compact
			branch();
		}

		fprintf(this->log,"Performing local gather\n");
		fflush(this->log);

		//Local Gather
		//Trace over all observables
		for(int o=0;o<this->observablesCollection.size();o++)
		{
			//Trace over all walkers and accumulate required observable
			for(int t=0;t<this->walkers.walkerCount;t++)
			{
				//fprintf(log,"%d %d\n",o,t);
				//fflush(log);
				this->walkers[t]->observablesCollection[o]->gather((void *) this->observablesCollection[o]);
			}
		}

		fprintf(this->log,"Local gather done\n");
		fflush(this->log);

		//Global Gather
		for(int o=0;o<this->MPIobservablesCollection.size();o++)
			this->MPIobservablesCollection[o]->parallelSend();

		fprintf(this->log,"Parallel gather done\n");
		fflush(this->log);

		//Acceptance Ratios
//		double accept = mover->getAccept();
//		fprintf(this->log,"------------------------------------------------------\n");
//		fprintf(this->log,"Acceptance Ratio for Current Bin: %10.6e\n",accept*divisor);
//		fprintf(this->log,"------------------------------------------------------\n");
//		fflush(this->log);
//		gaccept += accept;
	}

//	fprintf(this->log,"------------------------------------------------------\n");
//	fprintf(this->log,"Global Acceptance Ratio: %10.6e\n",gaccept*divisor/this->runParams.bins);
//	fprintf(this->log,"------------------------------------------------------\n");
//	fflush(this->log);

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

}

///////////////////////////////////////////////////////////////////////////////
template void MPIBasicRunner<int>::initialize();
template void MPIBasicRunner<int>::masterRun();
template void MPIBasicRunner<int>::run();
template void MPIBasicRunner<int>::finalize();
template void MPIBasicRunner<int>::masterFinalize();
template void MPIBasicRunner<int>::branch();

} /* namespace runners */
