/*
 * AutoCorr.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */

#ifndef AUTOCORR_CPP_
#define AUTOCORR_CPP_

#include "dmc.h"

namespace measures {

template <class T,class U>
void AutoCorr<T,U>::display()
{
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"Auto-Correlation Observable\n");
	fprintf(this->log,"==============================================\n");
	for(int i = 0;i<Wcollection.size();i++)
		fprintf(this->log,"%16.10e\n",Wcollection[i]);
	fprintf(this->log,"==============================================\n");
}

template <class T,class U>
void AutoCorr<T,U>::measure() {
	this->Wcollection.push_back(this->state.dweight);
}

template <class T,class U>
void AutoCorr<T,U>::writeViaIndex(int idx) {
#if DEBUG >= 2
	fprintf(this->log,"AutoCorr Write\n");
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

	double div = 1.0/this->Zcount;
	for(int i=0;i<lsize;i++)
	{
		autocorr[i] *= div;
		autocorr2[i] *= div;
		autocorr2[i] = sqrt((autocorr2[i] - autocorr[i]*autocorr[i])/(this->Zcount-1));

		wif<<i<<"\t"<<autocorr[i]<<"\t"<<autocorr2[i]<<endl;
	}
	wif.close();

	//Reset for next bin
	clear();
}

template <class T,class U>
void AutoCorr<T,U>::clear()
{
	Zcount = 0;
	if(allocated)
	{
		allocated = false;
		lsize = 0;
		delete[] autocorr;
		delete[] autocorr2;
	}

	this->Wcollection.clear();
}

template <class T,class U>
void AutoCorr<T,U>::gather(void* p)
{
	AutoCorr<T,U>* obj = (AutoCorr<T,U>*)p;

	int l = this->Wcollection.size();

	if(!obj->allocated)
	{
		obj->lsize = l;
		obj->autocorr = new double[l];
		obj->autocorr2 = new double[l];
		obj->allocated = true;

		for(int i = 0;i<l;i++)
			obj->autocorr[i] = obj->autocorr2[i] = 0.0;
	}

	double avg = 0.0;
	for(int i = 0;i<l;i++)
		avg += this->Wcollection[i];
	avg /= l;

	for(int i = 0;i<l;i++)
	{
		//Compute local a/c
		double lval = 0.0;
		for(int j = i;j<l;j++)
			lval += this->Wcollection[j-i]*this->Wcollection[j];
		lval /= (l-i);

		double t = (lval-avg*avg);
		obj->autocorr[i] += t;
		obj->autocorr2[i] += t*t;
	}
	obj->Zcount += 1;

	this->Wcollection.clear();
}

template <class T,class U>
Observable<T,U>* AutoCorr<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	AutoCorr<T,U>* newo = new AutoCorr<T,U>(this->processId,this->procCount,ws,
			this->baseFileName,this->log);
	newo->Wcollection = this->Wcollection;
	return newo;
}

template <class T,class U>
void AutoCorr<T,U>::copy(void* p)
{
	AutoCorr<T,U>* obj = (AutoCorr<T,U>*)p;
	this->Wcollection.clear();
	this->Wcollection = obj->Wcollection;
}

template <class T,class U>
int AutoCorr<T,U>::parallelSend()
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
	fprintf(this->log,"AutoCorr Received notification from master\n");
	fflush(this->log);
#endif
	//First send size
	MPI_Send(&this->Zcount,1,MPI_INT,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->lsize,1,MPI_INT,0,tag,MPI_COMM_WORLD);
	MPI_Send(this->autocorr,lsize,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(this->autocorr2,lsize,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);

	clear();
	return SUCCESS;
}

template <class T,class U>
int AutoCorr<T,U>::parallelReceive()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	int tag, recv;
	MPI_Status stat;

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);

	//Get from 1 for size
	MPI_Send(&recv,1,MPI_INT,1,tag,MPI_COMM_WORLD);

	int lzc;
	MPI_Recv(&lzc,1,MPI_INT,1,tag,MPI_COMM_WORLD,&stat);
	this->Zcount += lzc;

	MPI_Recv(&lsize,1,MPI_INT,1,tag,MPI_COMM_WORLD,&stat);
	double* temp = new double[lsize];
	if(!allocated)
	{
		allocated = true;
		this->autocorr = new double[lsize];
		this->autocorr2 = new double[lsize];
		MPI_Recv(this->autocorr,lsize,MPI_DOUBLE,1,tag,MPI_COMM_WORLD,&stat);
		MPI_Recv(this->autocorr2,lsize,MPI_DOUBLE,1,tag,MPI_COMM_WORLD,&stat);
	}
	else
	{
		MPI_Recv(temp,lsize,MPI_DOUBLE,1,tag,MPI_COMM_WORLD,&stat);
		for(int i = 0;i<lsize;i++)
			this->autocorr[i] += temp[i];

		MPI_Recv(temp,lsize,MPI_DOUBLE,1,tag,MPI_COMM_WORLD,&stat);
		for(int i = 0;i<lsize;i++)
			this->autocorr2[i] += temp[i];
	}

	for(int procId=2;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
#if DEBUG >= 2
		fprintf(this->log,"AutoCorr Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif
		MPI_Recv(&lzc,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		this->Zcount += lzc;

		MPI_Recv(&lsize,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		MPI_Recv(temp,lsize,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		for(int i = 0;i<lsize;i++)
			this->autocorr[i] += temp[i];

		MPI_Recv(temp,lsize,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		for(int i = 0;i<lsize;i++)
			this->autocorr2[i] += temp[i];

	}

	delete[] temp;
	return SUCCESS;
}


template <class T,class U>
void AutoCorr<T,U>::serialize(Serializer<U>& obj)
{
	obj<<Wcollection;
}

template <class T,class U>
void AutoCorr<T,U>::unserialize(Serializer<U>& obj)
{
	Wcollection.clear();
	obj>>Wcollection;
}
///////////////////////////////////////////////////////////////////////////

template void AutoCorr<int,stringstream>::measure();
template void AutoCorr<int,stringstream>::writeViaIndex(int idx);
template void AutoCorr<int,stringstream>::clear();
template void AutoCorr<int,stringstream>::gather(void* p);
template Observable<int,stringstream>* AutoCorr<int,stringstream>::duplicate(core::WalkerState<int,stringstream>&);
template void AutoCorr<int,stringstream>::copy(void* p);
template int AutoCorr<int,stringstream>::parallelSend();
template int AutoCorr<int,stringstream>::parallelReceive();
template void AutoCorr<int,stringstream>::serialize(Serializer<stringstream>&);
template void AutoCorr<int,stringstream>::unserialize(Serializer<stringstream>&);

template void AutoCorr<float,stringstream>::measure();
template void AutoCorr<float,stringstream>::writeViaIndex(int idx);
template void AutoCorr<float,stringstream>::clear();
template void AutoCorr<float,stringstream>::gather(void* p);
template Observable<float,stringstream>* AutoCorr<float,stringstream>::duplicate(core::WalkerState<float,stringstream>&);
template void AutoCorr<float,stringstream>::copy(void* p);
template int AutoCorr<float,stringstream>::parallelSend();
template int AutoCorr<float,stringstream>::parallelReceive();
template void AutoCorr<float,stringstream>::serialize(Serializer<stringstream>&);
template void AutoCorr<float,stringstream>::unserialize(Serializer<stringstream>&);
} /* namespace measures */


#endif /* AUTOCORR_CPP_ */
