/*
 * CloneMultiplicity.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
 */

#include "dmc.h"

namespace measures {

template <class T,class U>
void CloneMultiplicity<T,U>::measure() {

}

template <class T, class U>
void CloneMultiplicity<T,U>::writeViaIndex(int idx) {
#if DEBUG >= 2
	fprintf(this->log,"CloneMultiplicity Write\n");
	fflush(this->log);
#endif
	stringstream s;
	s<<idx;
	string fname = this->baseFileName + "_" + s.str();

	ofstream wif(fname,std::ofstream::out);
	int lsize = idc.size();
	for(int i = 0;i<lsize;i++)
		wif<<i<<" "<<idc[i].size()<<endl;
	wif.close();
	clear();
}

template <class T, class U>
void CloneMultiplicity<T,U>::clear()
{
	idc.clear();
}

template <class T, class U>
void CloneMultiplicity<T,U>::branchGather(void* p)
{
#ifndef NOBRANCH
	CloneMultiplicity<T,U>* obj = (CloneMultiplicity<T,U>*)p;
	int lsize = this->state.idhistory.size();

	if(obj->idc.size()<lsize)
		obj->idc.resize(lsize);

	for(int i=0;i<lsize;i++)
		obj->idc[i].insert(this->state.idhistory[i]);
#endif
}

template <class T, class U>
void CloneMultiplicity<T,U>::gather(void* p)
{
	//Store into global accumulator
	CloneMultiplicity<T,U>* obj = (CloneMultiplicity<T,U>*)p;
}


template <class T, class U>
void CloneMultiplicity<T,U>::copy(void* p)
{
	CloneMultiplicity<T,U>* obj = (CloneMultiplicity<T,U>*)p;
}

template <class T, class U>
Observable<T,U>* CloneMultiplicity<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	CloneMultiplicity<T,U>* newo = new CloneMultiplicity<T,U>(this->processId,this->procCount,this->totalWalkers,ws,
			this->baseFileName,this->log);
	return newo;
}

template <class T, class U>
void CloneMultiplicity<T,U>::display()
{
}

template <class T, class U>
int CloneMultiplicity<T,U>::parallelSend()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	int tag, recv;
	MPI_Status stat;

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);

	//Wait to be notified by master
	MPI_Recv(&recv,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
#if DEBUG >= 2
	fprintf(this->log,"CloneMultiplicity Received notification from master\n");
	fflush(this->log);
#endif

#ifndef NOBRANCH
	int lsize = this->idc.size();
	MPI_Send(&lsize,1,MPI_INT,0,tag,MPI_COMM_WORLD);
	for(int i = 0;i<lsize;i++)
	{
		int ssize = idc[i].size();
		MPI_Send(&ssize,1,MPI_INT,0,tag,MPI_COMM_WORLD);
		int* data = new int[ssize];
		int j = 0;
		for(set<int>::iterator it = idc[i].begin();it!=idc[i].end();++it)
			data[j++] = *it;
		MPI_Send(data,ssize,MPI_INT,0,tag,MPI_COMM_WORLD);
		delete[] data;
	}
#endif

#if DEBUG >= 2
	fprintf(this->log,"CloneMultiplicity Transfer Complete\n");
	fflush(this->log);
#endif

	clear();
	return SUCCESS;
}

template <class T, class U>
int CloneMultiplicity<T,U>::parallelReceive()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	clear();

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);
	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
#if DEBUG >= 2
		fprintf(this->log,"CloneMultiplicity Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif

#ifndef NOBRANCH
		int lsize = 0;
		MPI_Recv(&lsize,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		if(idc.size()<lsize)
			idc.resize(lsize);

		for(int i = 0;i<lsize;i++)
		{
			int ssize;
			MPI_Recv(&ssize,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
			int* ldata = new int[ssize];
			MPI_Recv(ldata,ssize,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
			for(int j=0;j<ssize;j++)
				idc[i].insert(ldata[j]);
			delete[] ldata;
		}
#endif

#if DEBUG >= 2
		fprintf(this->log,"CloneMultiplicity Finished receiving from process:%d\n",procId);
		fflush(this->log);
#endif
	}

#ifndef NOBRANCH
	ofstream wif(this->baseFileName + "P",std::ofstream::app);
	int lsize = idc.size();
	for(int i = 0;i<lsize;i++)
	{
		wif<<i<<" "<<idc[i].size()<<endl;
	}
	wif<<"---------------------------------\n";
	wif.close();
#endif

}

template <class T, class U>
void CloneMultiplicity<T,U>::serialize(Serializer<U>& obj)
{

}

template <class T, class U>
void CloneMultiplicity<T,U>::unserialize(Serializer<U>& obj)
{

}
///////////////////////////////////////////////////////////////////////////

template void CloneMultiplicity<int,stringstream>::measure();
template void CloneMultiplicity<int,stringstream>::writeViaIndex(int idx);
template void CloneMultiplicity<int,stringstream>::clear();
template void CloneMultiplicity<int,stringstream>::gather(void* p);
template Observable<int,stringstream>* CloneMultiplicity<int,stringstream>::duplicate(core::WalkerState<int,stringstream>&);
template void CloneMultiplicity<int,stringstream>::copy(void* p);
template int CloneMultiplicity<int,stringstream>::parallelSend();
template int CloneMultiplicity<int,stringstream>::parallelReceive();
template void CloneMultiplicity<int,stringstream>::serialize(Serializer<stringstream>&);
template void CloneMultiplicity<int,stringstream>::unserialize(Serializer<stringstream>&);

template void CloneMultiplicity<float,stringstream>::measure();
template void CloneMultiplicity<float,stringstream>::writeViaIndex(int idx);
template void CloneMultiplicity<float,stringstream>::clear();
template void CloneMultiplicity<float,stringstream>::gather(void* p);
template Observable<float,stringstream>* CloneMultiplicity<float,stringstream>::duplicate(core::WalkerState<float,stringstream>&);
template void CloneMultiplicity<float,stringstream>::copy(void* p);
template int CloneMultiplicity<float,stringstream>::parallelSend();
template int CloneMultiplicity<float,stringstream>::parallelReceive();
template void CloneMultiplicity<float,stringstream>::serialize(Serializer<stringstream>&);
template void CloneMultiplicity<float,stringstream>::unserialize(Serializer<stringstream>&);

} /* namespace measures */


