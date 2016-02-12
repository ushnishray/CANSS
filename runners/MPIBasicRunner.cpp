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
void MPIBasicRunner<T>::displayBranchStat(int nbranches)
{
	fprintf(this->log,"\n***************************************\n");
	fprintf(this->log,"Total branches: %d\n",nbranches);
	fprintf(this->log,"No. of splits: %d\n",this->nclones);
	fprintf(this->log,"No. of eliminates: %d\n",this->nelims);
}

template <class T>
void MPIBasicRunner<T>::initialize()
{

}

template <class T>
void MPIBasicRunner<T>::masterRun()
{
	int msg,tag;
	MPI_Status stat;

	for(int m=0;m<this->runParams.bins;m++)
	{
		fprintf(this->log,"Bin: %d\n",m);
		fflush(this->log);

		//Expect branching
		do{
			//Bin done
			for(int procId=1;procId<this->procCount;procId++)
				MPI_Recv(&msg,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);

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
//			fprintf(this->log,"Writing observable %d.\n",o);
//			fflush(this->log);
			this->observablesCollection[o]->writeViaIndex(m);
		}

		fprintf(this->log,"Write Done.\n");
		fflush(this->log);

		//Synchronize
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

template <class T>
void MPIBasicRunner<T>::run()
{
//	double gaccept = 0.0;
//	double divisor = 1.0/this->runParams.nsteps/this->walkers.walkerCount/walkers[0]->state.particleCount;

	int eqbranch = this->runParams.eSteps*0.90*BRANCHPERCENT;
	eqbranch = (eqbranch<1) ? 1 : eqbranch;
	int databranch = runParams.nsteps*0.90*BRANCHPERCENT;
	databranch = (databranch<1) ? 1 : databranch;

	for(int m=0;m<this->runParams.bins;m++)
	{
		int nbranches = 0;
		walkers.resetWalkers();

		fprintf(this->log,"Starting Bin: %d with %d walkers\n",m,walkers.walkerCount);
		fflush(this->log);
		//Initalize walkers
		for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
			mover->initialize(it->second);
			
		//Do equilibration
#if defined CONSTPOPBRANCH
		//Spend 10% of time getting to steady state
		for(int i=0;i<this->runParams.eSteps*0.10;i++)
		{
			for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				mover->move(it->second);
		}

		//Spend 90% of time getting to weighted distribution
		for(int i=0;i<this->runParams.eSteps*0.90;i++)
		{
			for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				mover->move(it->second);

			if((i+1)%eqbranch == 0)
			{
				nbranches++;
				branch();
			}
		}
#else
		for(int i=0;i<this->runParams.eSteps;i++)
		{
			for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				mover->move(it->second);
		}
#endif

		//Reset walker times
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
			//branch
			if((i+1)%databranch == 0)
			{
				nbranches++;
				branch();
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
			for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
				it->second->observablesCollection[o]->gather((void *) this->observablesCollection[o]);
		}

		fprintf(this->log,"Local gather done\n");
		fflush(this->log);

		//Master needs to be notified that process has finished a bin
		int tag, msg = MPIBINDONE;
		MPI_Send(&msg,1,MPI_INT,0,tag,MPI_COMM_WORLD);

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

#if !defined CONSTPOPBRANCH
template <class T>
void MPIBasicRunner<T>::branch()
{
#if defined CONSTPOPBRANCH
	paircomp paircompobj;
	vector<pair<int,double>> widx;

	Weight localNorm;
	int scount = 0, ccount = 0;

	for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
	{
		widx.push_back(pair<int,double>(it->first,it->second->state.weight.logValue()));
		localNorm.add(it->second->state.weight);
#if defined CPB1
		if(it->second->state.weight.logValue() >= walkers.maxValue)
			scount++;
		else if(it->second->state.weight.logValue() <= walkers.minValue)
			ccount++;
#endif
	}

	//Merging involves 2 files so round up to an even number
	ccount = int(ccount/2+0.5)*2;

	//walkers to keep untouched
	int kcount = walkers.walkerCount - scount - ccount;

#if defined CPB1
	int dpop = kcount + 2*scount + ccount - walkers.walkerCount;
#elif defined CPB2
	//int dpop = kcount + 2*scount - walkers.walkerCount;
#endif

	sort(widx.begin(),widx.end(),paircompobj);

	//if dpop>0 then there are more splits than merges so we need to remove low weights
	//if dpop<0 then there are more merges than splits so we need to split large weights
#if defined CPB1
	if(dpop>0)
		ccount = 2*scount;
	else
		scount = ccount/2;
#elif defined CPB2
	//Eliminate 25% of the low weighted walkers
	//Provided they have a small weight
	scount = ccount = walkers.walkerCount*0.25;
	Weight discardedwt;
	int tsi = 0;
	for(int i = widx.size()-1;i>=widx.size()-scount;i--,tsi++)
	{
		discardedwt.add(walkers[widx[i].first]->state.weight);
		double cprob = discardedwt.logValue() - localNorm.logValue();
		if(cprob>=MINBRANCHWEIGHT)
			break;
	}
	scount = ccount = tsi;

#endif
	this->nclones += scount;
	this->nelims += ccount;

	//Now merge and split. For every 2 merged we can do 1 split
	int p = 0;
#if defined CPB1
	for(int i = widx.size()-ccount;i<widx.size();i+=2)
#elif defined CPB2
	for(int i = widx.size()-ccount;i<widx.size();i++)
#endif
	{
#if 0
		//This SHOULD WORK but unfortunately c++ garbage collection is hideously bad
		//Ends in memory over-allocation :(
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
#else
#if defined CPB1
		walkers[widx[p].first]->state.weight.update(0.5);

		double w1 = walkers[widx[i].first]->state.weight.logValue();
		double w2 = walkers[widx[i+1].first]->state.weight.logValue();
		int idx = 0;
		if(gsl_rng_uniform(walkers[i]->rgenref)<w1/(w1+w2))
			idx = widx[i].first;
		else
			idx = widx[i+1].first;
		(*walkers.walkerCollection)[idx]->copy(*walkers[widx[p].first]);
#elif defined CPB2
		(*walkers.walkerCollection)[widx[i].first]->copy(*walkers[widx[p].first]);
#endif
#endif
		p++;
	}

	widx.clear();
	//walkers.displayWalkers(this->log);
#endif
}
#else
template <class T>
void MPIBasicRunner<T>::branch()
{
	////////////////////////////////////////////////////////////
	//Tell Master to Start expecting branch requests
	////////////////////////////////////////////////////////////
	int msg = MPIBRANCH, tag;
	MPI_Send(&msg,1,MPI_INT,0,tag,MPI_COMM_WORLD);
	////////////////////////////////////////////////////////////

	int* idx = new int[walkers.walkerCount];
	int* ni = new int[walkers.walkerCount];
	bool* zvals = new bool[walkers.walkerCount];

	int Nw = walkers.walkerCount*this->procCount; //Total no. of walkers

	//Send to master
	Weight lw;
	for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
		lw.add(it->second->state.weight);
	lw.mpiSend(0);
	lw.resetValue();
	lw.mpiBcast(0);

	//Now compute new population of walkers and figure out all the walkers that will be written over
	int i = 0;
	int ccount = 0, zc=0;
	for(typename NumMap<Walker<T>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
	{
		idx[i] = it->first;
		ni[i] = (it->second->state.weight/lw).value()*Nw;

		//Keep track of zero values
		if(ni[i]==0)
			zvals[zc++] = i;

		ccount += ni[i++];
	}

	//Tell master how many walkers are in excess or under
	int rem = ccount - walkers.walkerCount;
	MPI_Send(&rem,1,MPI_INT,0,tag,MPI_COMM_WORLD);

	if(rem==0) //Rare but happens and lucky if it does!
	{
		int p = 0;
		for(int i=0;i<walkers.walkerCount && p<zc;i++)
		{
			while(ni[i]>1)
			{
				(*walkers.walkerCollection)[zvals[p++]]->copy(*walkers[idx[i]]);
				ni[i]--;
			}
		}
		return;
	}

	////////////////////////////////////////////////////////////
	//Now start communicating with other processes
	////////////////////////////////////////////////////////////

	//Receive the update table
	MPI_Status stat;
	int tt;
	MPI_Recv(&tt,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
	int* updateTable = new int[2*tt];
	MPI_Recv(updateTable,2*tt,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);

	if(rem<0)
	{
		//Accept walkers from other nodes
		//This is going to be slow!


	}
	else
	{
		//Send walkers to other nodes
		//This is also going to be slow!

	}

	delete[] updateTable;
}
#endif

template <class T>
void MPIBasicRunner<T>::masterBranch()
{
	int tag;
	MPI_Status stat;
	Weight Z;

	//Accummulate Weights and Transmit
	for(int i = 1;i<this->procCount;i++)
		Z.mpiReceive(i);
	Z.mpiBcast(0);

	//Accummulate senders and receivers
	vector<int> sendProcs, sendProcCount, recvProcs, recvProcCount;
	for(int i = 1;i<this->procCount;i++)
	{
		int rem;
		MPI_Recv(&rem,1,MPI_INT,i,tag,MPI_COMM_WORLD,&stat);
		if(rem>0)
		{
			sendProcs.push_back(i);
			sendProcCount.push_back(rem);
		}
		else if(rem<0)
		{
			recvProcs.push_back(i);
			recvProcCount.push_back(-rem);
		}
	}

	//Assemble messages for senders and receivers
	vector<int>* mproc = new vector<int>[sendProcs.size()+recvProcs.size()];
	vector<int>* mcount = new vector<int>[sendProcs.size()+recvProcs.size()];
	int sendp = 0,recvp = 0;

	while(sendp<sendProcs.size() && recvp<recvProcs.size())
	{
		if(recvProcCount[recvp]<sendProcCount[sendp])
		{
			//Assign receivers which processors to expect data from and data count
			mproc[recvp+sendProcs.size()].push_back(sendProcs[sendp]);
			mcount[recvp+sendProcs.size()].push_back(recvProcCount[recvp]);

			//Assign senders which processors to send data to and data count
			mproc[sendp].push_back(recvProcs[recvp]);
			mcount[sendp].push_back(recvProcCount[recvp]);

			sendProcCount[sendp] -= recvProcCount[recvp];
			recvp++;
		}
		else
		{
			//Assign receivers which processors to expect data from and data count
			mproc[recvp+sendProcs.size()].push_back(sendProcs[sendp]);
			mcount[recvp+sendProcs.size()].push_back(sendProcCount[sendp]);

			//Assign senders which processors to send data to and data count
			mproc[sendp].push_back(recvProcs[recvp]);
			mcount[sendp].push_back(sendProcCount[sendp]);

			recvProcCount[recvp] -= sendProcCount[sendp];
			//sendProcCount[sendp] = 0;

			sendp++;
			if(recvProcCount[recvp]==0)
				recvp++;
		}
	}

	//Now inform the processes first senders and then receivers
	int tsp = sendProcs.size();
	for(int i=0;i<tsp;i++)
	{
		int tt = mproc[i].size();
		int* updateTable = new int[tt*2];
		memcpy(updateTable,mproc[i].data(),tt*sizeof(int));
		memcpy(updateTable+tt,mcount[i].data(),tt*sizeof(int));
		MPI_Send(&tt,1,MPI_INT,sendProcs[i],tag,MPI_COMM_WORLD);
		MPI_Send(updateTable,2*tt,MPI_INT,sendProcs[i],tag,MPI_COMM_WORLD);
		delete[] updateTable;
	}

	for(int i=0;i<recvProcs.size();i++)
	{
		int tt = mproc[tsp+i].size();
		int* updateTable = new int[tt*2];
		memcpy(updateTable,mproc[tsp+i].data(),tt*sizeof(int));
		memcpy(updateTable+tt,mcount[tsp+i].data(),tt*sizeof(int));
		MPI_Send(&tt,1,MPI_INT,sendProcs[i],tag,MPI_COMM_WORLD);
		MPI_Send(updateTable,2*tt,MPI_INT,recvProcs[i],tag,MPI_COMM_WORLD);
		delete[] updateTable;
	}

	//Clean up
	delete[] mproc;
	delete[] mcount;
}

///////////////////////////////////////////////////////////////////////////////
template void MPIBasicRunner<int>::initialize();
template void MPIBasicRunner<int>::masterRun();
template void MPIBasicRunner<int>::run();
template void MPIBasicRunner<int>::finalize();
template void MPIBasicRunner<int>::masterFinalize();
template void MPIBasicRunner<int>::branch();
template void MPIBasicRunner<int>::masterBranch();
template void MPIBasicRunner<int>::displayBranchStat(int);

} /* namespace runners */
