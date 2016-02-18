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

template <class T,class U>
void Whistogram<T,U>::display()
{
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"W-Histogram Observable\n");
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"ltime: %d\n",ltime);
	fprintf(this->log,"[log] local Weight: %16.10e\n",localWeight.logValue());
	for(int i = 0;i<Wcollection.size();i++)
		fprintf(this->log,"%16.10e\n",Wcollection[i]);
	fprintf(this->log,"==============================================\n");
}

template <class T,class U>
void Whistogram<T,U>::measure() {
	ltime = this->state.ltime;
	localWeight = this->state.weight;
//	fprintf(this->log,"%10.6e\n",localWeight.logValue()); 
}

template <class T,class U>
void Whistogram<T,U>::writeViaIndex(int idx) {
#if DEBUG >= 2
	fprintf(this->log,"Whistogram Write\n");
	fflush(this->log);
#endif
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

template <class T,class U>
void Whistogram<T,U>::clear()
{
	ltime = 0;
	localWeight.resetValue();
}

template <class T,class U>
void Whistogram<T,U>::gather(void* p)
{
	Whistogram<T,U>* obj = (Whistogram<T,U>*)p;
	obj->ltime = ltime;
	double it = 1.0/(ltime*dt);
	obj->Wcollection.push_back(localWeight.logValue());

	localWeight.resetValue();	
	ltime = 0;
}

template <class T,class U>
Observable<T,U>* Whistogram<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	Whistogram<T,U>* newo = new Whistogram<T,U>(this->processId,this->procCount,ws,
			this->baseFileName,this->log,this->dt);
	newo->ltime = this->ltime;
	newo->localWeight.copy(this->localWeight);
	newo->Wcollection = this->Wcollection;
	return newo;
}

template <class T,class U>
void Whistogram<T,U>::copy(void* p)
{
	Whistogram<T,U>* obj = (Whistogram<T,U>*)p;
	this->ltime = obj->ltime;
	this->localWeight.copy(obj->localWeight);
	this->Wcollection = obj->Wcollection;
}

template <class T,class U>
int Whistogram<T,U>::parallelSend()
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
	fprintf(this->log,"Whistogram Received notification from master\n");
	fflush(this->log);
#endif
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

template <class T,class U>
int Whistogram<T,U>::parallelReceive()
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
#if DEBUG >= 2
		fprintf(this->log,"Whistogram Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif
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
#if DEBUG >= 2
		fprintf(this->log,"Whistogram Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif
		MPI_Recv(data,psizes[procId-1],MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		data += psizes[procId-1];
	}
	
	delete[] psizes;
	return SUCCESS;
}


template <class T,class U>
void Whistogram<T,U>::serialize(Serializer<U>& obj)
{
	obj<<dt<<localWeight<<ltime<<Wcollection;
}

template <class T,class U>
void Whistogram<T,U>::unserialize(Serializer<U>& obj)
{
	obj>>dt>>localWeight>>ltime>>Wcollection;
}
///////////////////////////////////////////////////////////////////////////

template void Whistogram<int,stringstream>::measure();
template void Whistogram<int,stringstream>::writeViaIndex(int idx);
template void Whistogram<int,stringstream>::clear();
template void Whistogram<int,stringstream>::gather(void* p);
template Observable<int,stringstream>* Whistogram<int,stringstream>::duplicate(core::WalkerState<int,stringstream>&);
template void Whistogram<int,stringstream>::copy(void* p);
template int Whistogram<int,stringstream>::parallelSend();
template int Whistogram<int,stringstream>::parallelReceive();
template void Whistogram<int,stringstream>::serialize(Serializer<stringstream>&);
template void Whistogram<int,stringstream>::unserialize(Serializer<stringstream>&);

} /* namespace measures */


