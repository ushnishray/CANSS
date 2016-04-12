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

#ifndef NOBRANCH
	fname = this->baseFileName + "E_" + s.str();
	wif.open(fname);
	wif.precision(FIELDPRECISION);
	wif.width(FIELDWIDTH);
	wif.setf(FIELDFORMAT);
	wif.fill(' ');
	for(int i=0;i<Wacollection.size();i++)
		wif<<Wacollection[i]<<endl;
	wif.close();
#endif

	//Reset for next bin
	clear();
}

template <class T,class U>
void Whistogram<T,U>::clear()
{
	ltime = 0;
	localWeight.resetValue();
	this->Wcollection.clear();
#ifndef NOBRANCH
	this->Wacollection.clear();
#endif
}

template <class T,class U>
void Whistogram<T,U>::gather(void* p)
{
	Whistogram<T,U>* obj = (Whistogram<T,U>*)p;
	obj->ltime = ltime;
//	double it = 1.0/(ltime*dt);
	obj->Wcollection.push_back(localWeight.value());

#ifdef NOBRANCH
	localWeight.resetValue();	
	ltime = 0;
#endif
}

template <class T,class U>
void Whistogram<T,U>::branchGather(void* p)
{
	Whistogram<T,U>* obj = (Whistogram<T,U>*)p;
	obj->ltime = ltime;
	obj->Wacollection.push_back(localWeight.value());

	localWeight.resetValue();
	ltime = 0;
}

template <class T,class U>
Observable<T,U>* Whistogram<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	Whistogram<T,U>* newo = new Whistogram<T,U>(this->processId,this->procCount,this->totalWalkers,ws,
			this->baseFileName,this->log,this->dt);
	newo->ltime = this->ltime;
	newo->localWeight.copy(this->localWeight);
	newo->Wcollection = this->Wcollection;
#ifndef NOBRANCH
	newo->Wacollection = this->Wacollection;
#endif
	return newo;
}

template <class T,class U>
void Whistogram<T,U>::copy(void* p)
{
	Whistogram<T,U>* obj = (Whistogram<T,U>*)p;
	this->ltime = obj->ltime;
	this->localWeight.copy(obj->localWeight);
	this->Wcollection.clear();
	this->Wcollection = obj->Wcollection;
#ifndef NOBRANCH
	this->Wacollection.clear();
	this->Wacollection = obj->Wacollection;
#endif
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

#ifndef NOBRANCH
	MPI_Send(&this->Wacollection.front(),this->Wacollection.size(),MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
#endif

	ltime = 0;
	localWeight.resetValue();
	Wcollection.clear();
#ifndef NOBRANCH
	Wacollection.clear();
#endif
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

	vector<double> lWcollection, lWacollection;
	lWcollection.resize(tsize); //need to resize first
	double* data = lWcollection.data();
#ifndef NOBRANCH
	lWacollection.resize(tsize); //need to resize first
	double* data1 = lWacollection.data();
#endif

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

#ifndef NOBRANCH
		MPI_Recv(data1,psizes[procId-1],MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		data1 += psizes[procId-1];
#endif
	}
	
	delete[] psizes;

	//Compute variance per observation
	double x = 0.0, x2 = 0.0;
#ifndef NOBRANCH
	double y = 0.0, y2 = 0.0;
#endif
	for(int i = 0;i<tsize;i++)
	{
		x += lWcollection[i];
		x2 += lWcollection[i]*lWcollection[i];
#ifndef NOBRANCH
		y += lWacollection[i];
		y2 += lWacollection[i]*lWacollection[i];
#endif
	}
	x /= tsize; x2 /= tsize; x2 = (x2 - x*x)/tsize;
#ifndef NOBRANCH
	y /= tsize; y2 /= tsize; y2 = (y2 - y*y)/tsize;
#endif

#ifdef NOBRANCH
	ofstream wof(this->baseFileName + "_gather",std::ofstream::app);
	wof<<x<<"\t"<<sqrt(x2)<<endl;
	wof.close();
#else
	ofstream wof(this->baseFileName + "_gather",std::ofstream::app);
	wof<<x<<"\t"<<sqrt(x2)<<"\t"<<y<<"\t"<<sqrt(y2)<<endl;
	wof.close();
#endif

	//Local gathering done
	//Now put into global collector
	int osize = Wcollection.size();
	this->Wcollection.resize(osize+tsize);
	data = Wcollection.data() + osize;
	memcpy(data,lWcollection.data(),sizeof(double)*tsize);
#ifndef NOBRANCH
	osize = Wacollection.size();
	this->Wacollection.resize(osize+tsize);
	data1 = Wacollection.data() + osize;
	memcpy(data1,lWacollection.data(),sizeof(double)*tsize);
#endif

	return SUCCESS;
}


template <class T,class U>
void Whistogram<T,U>::serialize(Serializer<U>& obj)
{
	obj<<dt<<localWeight<<ltime<<Wcollection;
#ifndef NOBRANCH
	obj<<Wacollection;
#endif
}

template <class T,class U>
void Whistogram<T,U>::unserialize(Serializer<U>& obj)
{
	Wcollection.clear();
	obj>>dt>>localWeight>>ltime>>Wcollection;
#ifndef NOBRANCH
	Wacollection.clear();
	obj>>Wacollection;
#endif
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

template void Whistogram<float,stringstream>::measure();
template void Whistogram<float,stringstream>::writeViaIndex(int idx);
template void Whistogram<float,stringstream>::clear();
template void Whistogram<float,stringstream>::gather(void* p);
template Observable<float,stringstream>* Whistogram<float,stringstream>::duplicate(core::WalkerState<float,stringstream>&);
template void Whistogram<float,stringstream>::copy(void* p);
template int Whistogram<float,stringstream>::parallelSend();
template int Whistogram<float,stringstream>::parallelReceive();
template void Whistogram<float,stringstream>::serialize(Serializer<stringstream>&);
template void Whistogram<float,stringstream>::unserialize(Serializer<stringstream>&);
} /* namespace measures */


