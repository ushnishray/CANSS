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

		//Expect branching
		do{
			fprintf(this->log,"\nIn Control Loop\n");

			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Gather(&smsg,1,MPI_CHAR,rmsg,1,MPI_CHAR,0,MPI_COMM_WORLD);
			msg = rmsg[1];
			for(int i=2;i<procCount;i++)
			{
				msg &= rmsg[i];
				fprintf(this->log,"Received %d from %d\n",rmsg[i],i);
			}
			fflush(this->log);

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
				branch();
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
template <class T, class U>
void MPIBasicRunner<T,U>::branch()
{
#if defined CONSTPOPBRANCH
	paircomp paircompobj;
	vector<pair<int,double>> widx;

	Weight localNorm;
	int scount = 0, ccount = 0;

	for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
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
		Walker<T,U>* wcpy = walkers[widx[p].first]->duplicate();

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
template <class T, class U>
void MPIBasicRunner<T,U>::branch()
{

	fprintf(this->log,"\nBranching started: %d\n",branchcount++);
	////////////////////////////////////////////////////////////
	//Tell Master to Start expecting branch requests
	////////////////////////////////////////////////////////////
	char smsg = MPIBRANCH, rmsg = 0;
	int tag = 0, rem;
	MPI_Status stat;

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(&smsg,1,MPI_CHAR,&rmsg,1,MPI_CHAR,0,MPI_COMM_WORLD);
#if DEBUG >= 3
	fprintf(this->log,"Notified master to start expecting data.\n");
	fflush(this->log);
#endif
	////////////////////////////////////////////////////////////

	int* idx = new int[walkers.walkerCount];
	int* ni = new int[walkers.walkerCount];
	int* zvals = new int[walkers.walkerCount];

	int Nw = walkers.walkerCount*(this->procCount-1); //Total no. of walkers recall master has no walkers

	//Send to master
	Weight lw;
	for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
		lw.add(it->second->state.weight);
	MPI_Recv(&rem,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
#if DEBUG >= 3
	fprintf(this->log,"Received Z information send request.\n");
#endif
	lw.mpiSend(0);
#if DEBUG >= 3
	fprintf(this->log,"Sent Z information. log(Value) was: %10.6e\n",lw.logValue());
#endif
	lw.resetValue();
	lw.mpiBcast(0);
#if DEBUG >= 3
	fprintf(this->log,"BCast from 0. log(Value) is %10.6e\n",lw.logValue());
	fflush(this->log);
#endif

	//Now compute new population of walkers and figure out all the walkers that will be written over
	int i = 0;
	int ccount = 0, zc=0;
//	fprintf(this->log,"++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
	{
		idx[i] = it->first;
		double probi = (it->second->state.weight/lw).value();
		ni[i] = int(probi*Nw+0.5);

//		fprintf(this->log,"Walker %d, new population %d [prob %10.6e]\n",idx[i],ni[i],probi);

		//Keep track of zero values
		if(ni[i]==0)
			zvals[zc++] = idx[i];

		ccount += ni[i++];
	}
//	fprintf(this->log,"++++++++++++++++++++++++++++++++++++++++++++++++++\n");
//	fprintf(this->log,"New count of walkers: %d\n",ccount);
//	for(int i = 0;i<zc;i++)
//		fprintf(this->log,"Zero value %d at %d\n",i,zvals[i]);
//	walkers.displayWalkers(this->log);

	//Tell master how many walkers are in excess or under
	MPI_Recv(&rem,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
	fprintf(this->log,"Received delta{pop} sending request.\n");
	rem = ccount - walkers.walkerCount;
	fprintf(this->log,"Excess/Deficient walkers: %d\n",rem);
	MPI_Send(&rem,1,MPI_INT,0,tag,MPI_COMM_WORLD);
	fprintf(this->log,"Sent delta{pop}.\n");
	fflush(this->log);

	if(rem==0) //Rare but happens and lucky if it does!
	{
		int p = 0;
		for(int i=0;i<walkers.walkerCount && p<zc;i++)
		{
			if(ni[i]>1) //Have to do this so that weight update doesn't end up getting stuck
			{
				double uw = 1.0/ni[i];
				walkers[idx[i]]->state.weight.multUpdate(uw);

				while(ni[i]>1)
				{
					(*walkers.walkerCollection)[zvals[p++]]->copy(*walkers[idx[i]]);
					ni[i]--;
				}
			}
		}
	}
	else
	{
		////////////////////////////////////////////////////////////
		//Now start communicating with other processes
		////////////////////////////////////////////////////////////

		//Receive the update table
		int tt;
		MPI_Recv(&tt,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
		int* updateTable = new int[2*tt];
		MPI_Recv(updateTable,2*tt,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
		fprintf(this->log,"Received update table of size %d from 0\n",tt);
		fflush(this->log);

		if(rem<0)
		{
			//Accept walkers from other nodes
			//This is going to be slow!
#if DEBUG >= 3
			fprintf(this->log,"Will accept from:\n");
			for(int i=0;i<tt;i++)
				fprintf(this->log,"%d %d walkers\n",updateTable[i],updateTable[i+tt]);
			fflush(this->log);
#endif
			//Now receive from other procs
			int p = 0;
			Serializer<stringstream> ser;
			for(int i=0;i<tt;i++)
			{
				int tsize = 0;
				MPI_Send(&tsize,1,MPI_INT,updateTable[i],tag,MPI_COMM_WORLD); //Notify and then
				MPI_Recv(&tsize,1,MPI_INT,updateTable[i],tag,MPI_COMM_WORLD,&stat); //Receive
				char* data = new char[tsize];
				MPI_Recv(data,tsize,MPI_CHAR,updateTable[i],tag,MPI_COMM_WORLD,&stat);

#if DEBUG >= 3
				fprintf(this->log,"Data received from %d of size %d:\n",updateTable[i],tsize);
				fflush(this->log);
#endif

				ser.write(data,tsize);
				while(ser && p<zc)
				{
					int copies;
					ser>>copies;
					if(!ser)
						break;

#if DEBUG >= 3
					fprintf(this->log,"Copies: %d\n",copies);
					//fprintf(this->log,"Copying into %d from proc\n",zvals[p]);
					fflush(this->log);
#endif
					Walker<T,U>* copyw = walkers[zvals[p++]];
					copyw->unserialize(ser);
					//fflush(this->log);
					//copyw->display();
					for(int i = 1;i<copies;i++)
					{
						//fprintf(this->log,"Copying into %d from proc\n",zvals[p]);
						walkers[zvals[p++]]->copy(*copyw);
						//fflush(this->log);
					}
					//fflush(this->log);

					if(!ser)
						break;
				}
				ser.str("");
				ser.clear();
				delete[] data;

			}

			//Redistribute locally as needed
			for(int i=0;i<walkers.walkerCount && p<zc;i++)
			{
				if(ni[i]>1)
				{
					double uw = 1.0/ni[i];
					walkers[idx[i]]->state.weight.multUpdate(uw);
					while(ni[i]>1)
					{
						//fprintf(this->log,"Copying into %d from %d\n",zvals[p],idx[i]);
						(*walkers.walkerCollection)[zvals[p++]]->copy(*walkers[idx[i]]);
						ni[i]--;
					}
				}
			}
			//fflush(this->log);

		}
		else
		{
			//Send walkers to other nodes
			//This is also going to be slow!
#if DEBUG >= 3
			fprintf(this->log,"Will send to:\n");
			for(int i=0;i<tt;i++)
				fprintf(this->log,"%d %d walkers\n",updateTable[i],updateTable[i+tt]);
			fflush(this->log);
#endif

			//Also update weight
			int* nisend = new int[walkers.walkerCount];
			for(int i = 0;i<walkers.walkerCount;i++)
			{
				nisend[i] = 0;
				if(ni[i]>1)
				{
					double uw = 1.0/ni[i];
					walkers[idx[i]]->state.weight.multUpdate(uw);
				}
			}


			int count = 0;
			while(count<rem)
			{
				for(int i = 0;i<walkers.walkerCount && count<rem;i++)
				{
					if(ni[i]>1)
					{
						ni[i]--;
						nisend[i]++;
						count++;
					}
				}
			}

			/*
			fprintf(this->log,"++++++++ Distribution Send +++++++++++++++++++++++++++++++++\n");
			for(int i = 0;i<walkers.walkerCount;i++)
				fprintf(this->log,"Walker %d, local population %d send population %d\n",idx[i],ni[i],nisend[i]);
			fprintf(this->log,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
			fflush(this->log);
			*/

			//Accummulate and send
			int lcount = 0, startat = 0;
			for(int i=0;i<tt;i++)
			{
				Serializer<stringstream> ser;
				lcount = 0;
				for(int j=startat;j<walkers.walkerCount;j++)
				{
					//Skip 0's
					if(nisend[j]==0)
						continue;

					if(lcount+nisend[j]>updateTable[i+tt])
					{
						int tosend = updateTable[i+tt] - lcount;
						nisend[j] -= tosend;

						ser<<tosend;
						walkers[idx[j]]->serialize(ser);
#if DEBUG >= 3
						fprintf(this->log,"Sending to %d, %d copies of walker %d\n",updateTable[i],tosend,idx[j]);
#endif
						startat = j;
						break;
					}
					else if(lcount+nisend[j] == updateTable[i+tt])
					{
						ser<<nisend[j];
						walkers[idx[j]]->serialize(ser);
#if DEBUG >= 3
						fprintf(this->log,"Sending to %d, %d copies of walker %d\n",updateTable[i],nisend[j],idx[j]);
#endif
						startat = j+1;
						break;
					}

					lcount += nisend[j];
					//walkers[idx[j]]->display();
					ser<<nisend[j];
					walkers[idx[j]]->serialize(ser);
#if DEBUG >= 3
					fprintf(this->log,"Sending to %d, %d copies of walker %d\n",updateTable[i],nisend[j],idx[j]);
#endif
				}

				//Transmit Data
				ser.seekg(0,ser.end);
				unsigned int tsize = ser.tellg();
				ser.seekg(0,ser.beg);
				char* data = new char[tsize];
				ser.read(data,tsize); //Read from stream

				int temp;
				MPI_Recv(&temp,1,MPI_INT,updateTable[i],tag,MPI_COMM_WORLD,&stat); //Receive notification
				MPI_Send(&tsize,1,MPI_UNSIGNED,updateTable[i],tag,MPI_COMM_WORLD); //Send
				MPI_Send(data,tsize,MPI_CHAR,updateTable[i],tag,MPI_COMM_WORLD); //Send
#if DEBUG >= 3
				fprintf(this->log,"Data sent to %d of size %d:\n",updateTable[i],tsize);
				fflush(this->log);
#endif
				delete[] data;
			}

			//Now redistribute locally
			int p = 0;
			for(int i=0;i<walkers.walkerCount && p<zc;i++)
			{
				while(ni[i]>1)
				{
					//fprintf(this->log,"Copying %d into %d\n",idx[i],zvals[p]);
					(*walkers.walkerCollection)[zvals[p++]]->copy(*walkers[idx[i]]);
					ni[i]--;
				}
			}
			//fflush(this->log);
			delete[] nisend;
		}

		delete[] updateTable;
	}
#if DEBUG >= 3
	fprintf(this->log,"Notifying master that branching is done.\n");
	fflush(this->log);
#endif
	//walkers.displayWalkers(this->log);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(&smsg,1,MPI_CHAR,&rmsg,1,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
#if DEBUG >= 3
	fprintf(this->log,"Notified master that branching is done.\n");
	fprintf(this->log,"Branching done\n");
	fprintf(this->log,"========================================================\n");
	fflush(this->log);
#endif

#if 0
	//Send to master
	Weight ulw;
	for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
		ulw.add(it->second->state.weight);
	MPI_Recv(&rem,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
#if DEBUG >= 3
	fprintf(this->log,"Received Z information send request.\n");
#endif
	ulw.mpiSend(0);
#if DEBUG >= 3
	fprintf(this->log,"Sent Z information. log(Value) was: %10.6e\n",ulw.logValue());
#endif
	ulw.resetValue();
	ulw.mpiBcast(0);
#if DEBUG >= 3
	fprintf(this->log,"BCast from 0. log(Value) is %10.6e\n",lw.logValue());
	fflush(this->log);
#endif
#endif
	delete[] idx;
	delete[] ni;
	delete[] zvals;
}
#endif

template <class T, class U>
void MPIBasicRunner<T,U>::masterBranch()
{
	fprintf(this->log,"\nBranching started: %d\n",branchcount++);
	int tag = 0, rem;
	MPI_Status stat;
	Weight Z;
	//Accummulate Weights and Transmit
	for(int i = 1;i<this->procCount;i++)
	{
		MPI_Send(&rem,1,MPI_INT,i,tag,MPI_COMM_WORLD);
		Z.mpiReceive(i);
#if DEBUG >= 3
		fprintf(this->log,"Received Z information from %d\n",i);
#endif
	}
	Z.mpiBcast(0);
#if DEBUG >= 3
	fprintf(this->log,"BCast Z information.\n");
	fflush(this->log);
#endif
	int totalsend = 0, totalrecv = 0;
	//Accummulate senders and receivers
	vector<int> sendProcs, sendProcCount, recvProcs, recvProcCount;
	for(int i = 1;i<this->procCount;i++)
	{
		rem = MPISTATUSSTART;
		MPI_Send(&rem,1,MPI_INT,i,tag,MPI_COMM_WORLD);
		MPI_Recv(&rem,1,MPI_INT,i,tag,MPI_COMM_WORLD,&stat);
		if(rem>0)
		{
			sendProcs.push_back(i);
			sendProcCount.push_back(rem);
			totalsend += rem;
		}
		else if(rem<0)
		{
			recvProcs.push_back(i);
			recvProcCount.push_back(-rem);
			totalrecv += (-rem);
		}
	}
#if DEBUG >= 3
	fprintf(this->log,"Received delta{pop} information.\nTotal sends: %d, total recvs: %d\n",totalsend,totalrecv);
	fflush(this->log);
#endif

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
#if DEBUG >= 3
	fprintf(this->log,"Calculated delta{pop}.\n");
	fflush(this->log);
	fprintf(this->log,"+++++++++++++++++++++++++++++++++++++++++++++++++++\n");

	for(int i=0;i<sendProcs.size();i++)
		for(int j=0;j<mcount[i].size();j++)
			fprintf(this->log,"%d will send to %d, walkers = %d\n",sendProcs[i],mproc[i][j],mcount[i][j]);

	int shift = sendProcs.size();
	for(int i=0;i<recvProcs.size();i++)
		for(int j=0;j<mcount[i+shift].size();j++)
			fprintf(this->log,"%d will get from %d, walkers = %d\n",recvProcs[i],mproc[shift+i][j],mcount[shift+i][j]);
	fprintf(this->log,"+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif

	//Walker count assert
	int swalkers = 0, rwalkers = 0;
	for(int i=0;i<sendProcs.size();i++)
		for(int j=0;j<mcount[i].size();j++)
			swalkers += mcount[i][j];

	for(int i=sendProcs.size();i<sendProcs.size()+recvProcs.size();i++)
		for(int j=0;j<mcount[i].size();j++)
			rwalkers += mcount[i][j];

#if DEBUG >= 3
	fprintf(this->log,"Send/Receive Walkers: %d/%d\n",swalkers,rwalkers);
#endif
	assert(swalkers == rwalkers);

	//Now inform the processes first senders and then receivers
	int tsp = sendProcs.size();
	for(int i=0;i<tsp;i++)
	{
		int tt = mproc[i].size();
		int* updateTable = new int[tt*2];
		memcpy(updateTable,mproc[i].data(),tt*sizeof(int));
		memcpy(updateTable+tt,mcount[i].data(),tt*sizeof(int));

#if DEBUG >= 3
		fprintf(this->log,"Sending update table to: %d\n",sendProcs[i]);
#endif
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
#if DEBUG >= 3
		fprintf(this->log,"Sending update table to: %d\n",recvProcs[i]);
#endif
		MPI_Send(&tt,1,MPI_INT,recvProcs[i],tag,MPI_COMM_WORLD);
		MPI_Send(updateTable,2*tt,MPI_INT,recvProcs[i],tag,MPI_COMM_WORLD);
		delete[] updateTable;
	}

	//Clean up
	delete[] mproc;
	delete[] mcount;
	fflush(this->log);
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	char* rmsg = new char[this->procCount];
	char smsg = 0;
	for(int i = 0;i<procCount;i++)
		rmsg[i] = 0;
#if DEBUG >= 3
	fprintf(this->log,"Expecting to received info from processes to notify completion of branching.\n");
	fflush(this->log);
#endif
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(&smsg,1,MPI_CHAR,rmsg,1,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

#if DEBUG >= 3
	fprintf(this->log,"Received info from processes that branching is done.\n");
	fflush(this->log);
#endif
	delete[] rmsg;
	fprintf(this->log,"Branching done\n");
	fprintf(this->log,"========================================================\n");
	fflush(this->log);


#if 0
	//Accummulate Weights and Transmit
	Weight ZC;
	for(int i = 1;i<this->procCount;i++)
	{
		MPI_Send(&rem,1,MPI_INT,i,tag,MPI_COMM_WORLD);
		ZC.mpiReceive(i);
#if DEBUG >= 3
		fprintf(this->log,"Received Z information from %d\n",i);
#endif
	}
	ZC.mpiBcast(0);
#if DEBUG >= 3
	fprintf(this->log,"BCast Z information.\n");
	fflush(this->log);
#endif
#endif
}

///////////////////////////////////////////////////////////////////////////////
template void MPIBasicRunner<int,stringstream>::initialize();
template void MPIBasicRunner<int,stringstream>::masterRun();
template void MPIBasicRunner<int,stringstream>::run();
template void MPIBasicRunner<int,stringstream>::finalize();
template void MPIBasicRunner<int,stringstream>::masterFinalize();
template void MPIBasicRunner<int,stringstream>::branch();
template void MPIBasicRunner<int,stringstream>::masterBranch();
template void MPIBasicRunner<int,stringstream>::displayBranchStat(int);

} /* namespace runners */
