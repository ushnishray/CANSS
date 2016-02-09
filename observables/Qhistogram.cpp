/*
 * Qhistogram.cpp
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
void Qhistogram<T>::measure() {
	Q.x += this->state.dQ.x;
	Q.y += this->state.dQ.y;
	Q.z += this->state.dQ.z;

	ltime = this->state.ltime;
}

template <class T>
void Qhistogram<T>::writeViaIndex(int idx) {

	fprintf(this->log,"Qhistogram Write\n");
	fflush(this->log);

	stringstream s;
	s<<idx;
	string fname = this->baseFileName + "_" + s.str();
	
	ofstream wif(fname,std::ofstream::app);
	wif.precision(FIELDPRECISION);
	wif.width(FIELDWIDTH);
	wif.setf(FIELDFORMAT);
	wif.fill(' ');

	for(int i=0;i<Qcollection.size();i++)
		wif<<Qcollection[i].x<<" "<<Qcollection[i].y<<" "<<Qcollection[i].z<<endl;
	wif.close();

	//Reset for next bin
	Qcollection.clear();
	clear();
}

template <class T>
void Qhistogram<T>::clear()
{
	Q.x = 0.0;
	Q.y = 0.0;
	Q.z = 0.0;
	ltime = 0;
}

template <class T>
void Qhistogram<T>::gather(void* p)
{
	Qhistogram<T>* obj = (Qhistogram<T>*)p;
	obj->ltime = ltime;
	double it = 1.0/(ltime*dt);
	Q.x *= it;
	Q.y *= it;
	Q.z *= it;
	obj->Qcollection.push_back(Q);
	
	Q.x = Q.y = Q.z = 0;
	ltime = 0;
}

template <class T>
Observable<T>* Qhistogram<T>::duplicate(core::WalkerState<T>& ws)
{
	Qhistogram<T>* newo = new Qhistogram<T>(this->processId,this->procCount,ws,
			this->baseFileName,this->log,this->dt);
	newo->ltime = this->ltime;
	newo->Q.x = this->Q.x;
	newo->Q.y = this->Q.y;
	newo->Q.z = this->Q.z;
	newo->Qcollection = this->Qcollection;
	return newo;
}

template <class T>
void Qhistogram<T>::copy(void* p)
{
	Qhistogram<T>* obj = (Qhistogram<T>*) p;
	this->Qcollection.clear();

	this->ltime = obj->ltime;
	this->Q.x = obj->Q.x;
	this->Q.y = obj->Q.y;
	this->Q.z = obj->Q.z;
	this->Qcollection = obj->Qcollection;
}

template <class T>
int Qhistogram<T>::parallelSend()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	int tag, recv;
	MPI_Status stat;

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);

	//Wait to be notified by master
	MPI_Recv(&recv,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
//	fprintf(this->log,"Qhistogram Received notification from master\n");
//	fflush(this->log);

	//First send size	
	int lsize = this->Qcollection.size();
	MPI_Send(&lsize,1,MPI_INT,0,tag,MPI_COMM_WORLD);

	//Then wait for notification to send data
	MPI_Recv(&recv,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
	MPI_Send(&this->Qcollection.front(),this->Qcollection.size()*3,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);	

	ltime = 0;
	Q.x = Q.y = Q.z = 0;
	Qcollection.clear();
	return SUCCESS;
}

template <class T>
int Qhistogram<T>::parallelReceive()
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
//		fprintf(this->log,"Qhistogram Sending notification to process:%d\n",procId);
//		fflush(this->log);
		int lsize = 0;		
		MPI_Recv(&lsize,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		tsize += lsize;
		psizes[procId-1] = lsize;
	}

	this->Qcollection.resize(tsize); //need to resize first
	vect<double>* data = this->Qcollection.data();

	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
//		fprintf(this->log,"Qhistogram Sending notification to process:%d\n",procId);
//		fflush(this->log);
		MPI_Recv(data,psizes[procId-1]*3,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		data += psizes[procId-1];
	}
	
	delete[] psizes;
	return SUCCESS;
}

///////////////////////////////////////////////////////////////////////////

template void Qhistogram<int>::measure();
template void Qhistogram<int>::writeViaIndex(int idx);
template void Qhistogram<int>::clear();
template void Qhistogram<int>::gather(void* p);
template Observable<int>* Qhistogram<int>::duplicate(core::WalkerState<int>&);
template void Qhistogram<int>::copy(void* p);
template int Qhistogram<int>::parallelSend();
template int Qhistogram<int>::parallelReceive();

} /* namespace measures */


