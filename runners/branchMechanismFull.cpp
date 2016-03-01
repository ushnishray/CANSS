/*
 * branchMechanismFull.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */

#include "dmc.h"

template <class T, class U>
void MPIBasicRunner<T,U>::branchFull()
{
	int tag = 0, temp;
	MPI_Status stat;

	int* idx = new int[walkers.walkerCount];
	int* ni = new int[walkers.walkerCount];
	int* zvals = new int[walkers.walkerCount];

	int Nw = walkers.walkerCount*(this->procCount-1); //Total no. of walkers recall master has no walkers

	//Send to master
	Weight lw;
	for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
		lw.add(it->second->state.weight);
	MPI_Recv(&temp,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
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
	fflush(this->log);
#endif

	//Now compute new population of walkers and figure out all the walkers that will be written over
	int i = 0;
	int ccount = 0, zc=0;
	double tprob = 0.0;
#if DEBUG>=4
	fprintf(this->log,"++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif
	for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
	{
		idx[i] = it->first;
		double probi = (it->second->state.weight/lw).value(); tprob += probi;
		it->second->state.weight.resetValue(); //Reset for next round

		ni[i] = probi*Nw + gsl_rng_uniform(walkers[0]->rgenref);
#if DEBUG>=4
		fprintf(this->log,"Walker %d, new population %d [prob %10.6e]\n",idx[i],ni[i],probi);
#endif
		//Keep track of zero values
		if(ni[i]==0)
			zvals[zc++] = idx[i];

		ccount += ni[i++];
	}
#if DEBUG>=4
	fprintf(this->log,"++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(this->log,"New count of walkers: %d\n",ccount);
	fflush(this->log);
#elif DEBUG>=3
	fprintf(this->log,"Net prob of walkers: %10.6e\n",tprob);
	fflush(this->log);
#endif

//	for(int i = 0;i<zc;i++)
//		fprintf(this->log,"Zero value %d at %d\n",i,zvals[i]);
//	walkers.displayWalkers(this->log);

	//Tell master how many walkers are in excess or under
	MPI_Recv(&temp,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
#if DEBUG >= 3
	fprintf(this->log,"Received delta{pop} sending request.\n");
#endif
	int rem = ccount - walkers.walkerCount;
#if DEBUG >= 3
	fprintf(this->log,"Excess/Deficient walkers: %d\n",rem);
#endif
	MPI_Send(&rem,1,MPI_INT,0,tag,MPI_COMM_WORLD);
#if DEBUG >= 3
	fprintf(this->log,"Sent delta{pop}.\n");
	fflush(this->log);
#endif

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

			//Redistribute locally as needed
			int p = 0;
			for(int i=0;i<walkers.walkerCount && p<zc;i++)
			{
				while(ni[i]>1)
				{
					//fprintf(this->log,"Copying into %d from %d\n",zvals[p],idx[i]);
					(*walkers.walkerCollection)[zvals[p++]]->copy(*walkers[idx[i]]);
					ni[i]--;
				}
			}

			//Now receive from other procs
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


#if DEBUG >= 3
				fprintf(this->log,"Randomly copying remaining %d walkers\n",zc-p);
				fflush(this->log);
#endif
			//Now make up for fractional occupation by randomly increasing walker populations
			for(;p<zc;p++)
			{
				int ridx = gsl_rng_uniform_int(walkers[0]->rgenref,walkers.walkerCount);
				(*walkers.walkerCollection)[zvals[p]]->copy(*walkers[idx[ridx]]);
				ni[ridx]++;
			}

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

			int* nisend = new int[walkers.walkerCount];
			for(int i = 0;i<walkers.walkerCount;i++)
				nisend[i] = 0;

			//Update new populations due to sending
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
	char smsg,rmsg;
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

template <class T, class U>
void MPIBasicRunner<T,U>::masterBranchFull()
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

	//Update master's estimation of C.G.F.
	this->FreeEnergy.add(Z.logValue());

#if DEBUG >= 3
	fprintf(this->log,"BCast Z information.\n");
	fflush(this->log);
#endif

#if DEBUG >= 1
	//This debug is particularly useful to see how the total weight of walkers
	//changes over the course of the sampling. We always want this to increase
	fprintf(this->log,"BCast from 0. log(Value) is %10.6e\n",Z.logValue());
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////
template void MPIBasicRunner<int,stringstream>::branchFull();
template void MPIBasicRunner<int,stringstream>::masterBranchFull();
template void MPIBasicRunner<float,stringstream>::branchFull();
template void MPIBasicRunner<float,stringstream>::masterBranchFull();
