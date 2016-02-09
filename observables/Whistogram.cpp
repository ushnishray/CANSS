/*
 * Whistogram.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
 */

#include "dmc.h"

namespace measures {

template <class T>
void Whistogram<T>::measure() {
	ltime = this->state.ltime;
	localWeight = this->state.weight;
//	fprintf(this->log,"%10.6e\n",localWeight.logValue()); 
}

template <class T>
void Whistogram<T>::writeViaIndex(int idx) {

	fprintf(this->log,"Whistogram Write\n");
	fflush(this->log);

	stringstream s;
	s<<idx;
	string fname = this->baseFileName + "_" + s.str();
	
	ofstream wif(fname);
	wif.precision(FIELDPRECISION);
	wif.width(FIELDWIDTH);
	wif.setf(FIELDFORMAT);
	wif.fill(' ');

	for(int i=0;i<Wcollection.size();i++)
		wif<<Wcollection[i]<<endl;
	wif.close();

	//Reset for next bin
	Wcollection.clear();
	clear();
}

template <class T>
void Whistogram<T>::clear()
{
	ltime = 0;
	localWeight.resetValue();
}

template <class T>
void Whistogram<T>::gather(void* p)
{
	Whistogram<T>* obj = (Whistogram<T>*)p;
	obj->ltime = ltime;
	double it = 1.0/(ltime*dt);
	obj->Wcollection.push_back(localWeight.logValue());

	localWeight.resetValue();	
	ltime = 0;
}

template <class T>
Observable<T>* Whistogram<T>::duplicate(core::WalkerState<T>& ws)
{
	Whistogram<T>* newo = new Whistogram<T>(this->processId,this->procCount,ws,
			this->baseFileName,this->log,this->dt);
	newo->ltime = this->ltime;
	newo->localWeight.copy(this->localWeight);
	newo->Wcollection = this->Wcollection;
	return newo;
}

template <class T>
void Whistogram<T>::copy(void* p)
{
	Whistogram<T>* obj = (Whistogram<T>*)p;
	this->ltime = obj->ltime;
	this->localWeight.copy(obj->localWeight);
	this->Wcollection = obj->Wcollection;
}

template <class T>
int Whistogram<T>::parallelSend()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	int tag, recv;
	MPI_Status stat;

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);

	//Wait to be notified by master
	MPI_Recv(&recv,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
//	fprintf(this->log,"Whistogram Received notification from master\n");
//	fflush(this->log);

	//First send size	
	int lsize = this->Wcollection.size();
	MPI_Send(&lsize,1,MPI_INT,0,tag,MPI_COMM_WORLD);

	//Then wait for notification to send data
	MPI_Recv(&recv,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
	MPI_Send(&this->Wcollection.front(),this->Wcollection.size(),MPI_DOUBLE,0,tag,MPI_COMM_WORLD);	

	ltime = 0;
	localWeight.resetValue();
	Wcollection.clear();
	return SUCCESS;
}

template <class T>
int Whistogram<T>::parallelReceive()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);

	long tsize = 0;
	int* psizes = new int[procCount-1];
	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
//		fprintf(this->log,"Whistogram Sending notification to process:%d\n",procId);
//		fflush(this->log);
		int lsize = 0;		
		MPI_Recv(&lsize,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		tsize += lsize;
		psizes[procId-1] = lsize;
	}

	this->Wcollection.resize(tsize); //need to resize first
	double* data = this->Wcollection.data();

	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
//		fprintf(this->log,"Whistogram Sending notification to process:%d\n",procId);
//		fflush(this->log);
		MPI_Recv(data,psizes[procId-1],MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		data += psizes[procId-1];
	}
	
	delete[] psizes;
	return SUCCESS;
}

///////////////////////////////////////////////////////////////////////////

template void Whistogram<int>::measure();
template void Whistogram<int>::writeViaIndex(int idx);
template void Whistogram<int>::clear();
template void Whistogram<int>::gather(void* p);
template Observable<int>* Whistogram<int>::duplicate(core::WalkerState<int>&);
template void Whistogram<int>::copy(void* p);
template int Whistogram<int>::parallelSend();
template int Whistogram<int>::parallelReceive();

} /* namespace measures */


