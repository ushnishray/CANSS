/*
 * AutoCorrFull.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */

#ifndef AutoCorrFull_CPP_
#define AutoCorrFull_CPP_

#include "dmc.h"

namespace measures {

template <class T,class U>
void AutoCorrFull<T,U>::display()
{
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"Auto-Correlation Observable\n");
	fprintf(this->log,"==============================================\n");
	for(int i = 0;i<Wcollection.size();i++)
		fprintf(this->log,"%16.10e\n",Wcollection[i]);
	fprintf(this->log,"==============================================\n");
}

template <class T,class U>
void AutoCorrFull<T,U>::measure() {
	this->Wcollection.push_back(this->state.dweight);
}

template <class T,class U>
void AutoCorrFull<T,U>::writeViaIndex(int idx) {
#if DEBUG >= 1
	fprintf(this->log,"AutoCorrFull Write\n");
	fprintf(this->log,"AutoCorrFull Size: %d %d\n",ac.size(),this->lsize);
	fflush(this->log);
#endif


	//Compute autocorrelation function
	double* acl = new double[this->lsize];
	double* acl2 = new double[this->lsize];
	for(int i = 0;i<this->lsize;i++)
		acl[i] = acl2[i] = 0.0;

	for(int p=0;p<this->ac.size();p++)
	{
		for(int i = 0;i<this->lsize;i++)
		{
			//Compute local a/c
			double lavg1 = 0.0, lavg2 = 0.0;
			double lval = 0.0;
			for(int j = i;j<this->lsize;j++)
			{
				lval += this->ac[p][j-i]*this->ac[p][j];
				lavg1 += this->ac[p][j-i];
				lavg2 += this->ac[p][j];
			}
			lavg1 /= (lsize-i);
			lavg2 /= (lsize-i);
			lval /= (lsize-i);

			double t = (lval-lavg1*lavg2);
			//double t = lval;
			acl[i] += t;
			acl2[i] += t*t;
		}
	}

	double norm = 1.0/ac.size();
	for(int i = 0;i<this->lsize;i++)
	{
		acl[i] *= norm;
		acl2[i] *= norm;
		acl2[i] = sqrt((acl2[i]-acl[i]*acl[i])/(ac.size()-1));
	}

	stringstream s;
	s<<idx;
	string fname = this->baseFileName + "_" + s.str();

	ofstream wif(fname);
	wif.precision(FIELDPRECISION);
	wif.width(FIELDWIDTH);
	wif.setf(FIELDFORMAT);
	wif.fill(' ');

	for(int i=0;i<lsize;i++)
		wif<<i<<"\t"<<acl[i]<<"\t"<<acl2[i]<<endl;
	wif.close();

	//Reset for next bin
	clear();
}

template <class T,class U>
void AutoCorrFull<T,U>::clear()
{
	if(lsize>0)
	{
		for(int i = 0;i<ac.size();i++)
			delete[] ac[i];
		lsize = 0;
	}

	ac.clear();
	Wcollection.clear();
}

template <class T,class U>
void AutoCorrFull<T,U>::gather(void* p)
{
	AutoCorrFull<T,U>* obj = (AutoCorrFull<T,U>*)p;

	int l = this->Wcollection.size();
	double* a = new double[l];
	memcpy(a,this->Wcollection.data(),sizeof(double)*l);

	obj->ac.push_back(a);
	obj->lsize = l;

	this->Wcollection.clear();
}

template <class T,class U>
Observable<T,U>* AutoCorrFull<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	AutoCorrFull<T,U>* newo = new AutoCorrFull<T,U>(this->processId,this->procCount,ws,
			this->baseFileName,this->log);
	newo->Wcollection = this->Wcollection;
	return newo;
}

template <class T,class U>
void AutoCorrFull<T,U>::copy(void* p)
{
	AutoCorrFull<T,U>* obj = (AutoCorrFull<T,U>*)p;
	this->Wcollection.clear();
	this->Wcollection = obj->Wcollection;
}

template <class T,class U>
int AutoCorrFull<T,U>::parallelSend()
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
	fprintf(this->log,"AutoCorrFull Received notification from master\n");
	fflush(this->log);
#endif

	Serializer<U> ser;
	ser<<(int) this->ac.size()<<this->lsize;
	for(int i = 0;i<this->ac.size();i++)
		for(int j=0;j<this->lsize;j++)
			ser<<this->ac[i][j];

	ser.seekg(0,ser.end);
	unsigned int tsize = ser.tellg();
	ser.seekg(0,ser.beg);
	char* data = new char[tsize];
	ser.read(data,tsize); //Read from stream

	//First send size
	MPI_Send(&tsize,1,MPI_UNSIGNED,0,tag,MPI_COMM_WORLD); //Send
	MPI_Send(data,tsize,MPI_CHAR,0,tag,MPI_COMM_WORLD); //Send

	clear();
	return SUCCESS;
}

template <class T,class U>
int AutoCorrFull<T,U>::parallelReceive()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	int tag, recv;
	MPI_Status stat;

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);

	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
#if DEBUG >= 2
		fprintf(this->log,"AutoCorrFull Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif

		unsigned int tsize;
		MPI_Recv(&tsize,1,MPI_UNSIGNED,procId,tag,MPI_COMM_WORLD,&stat); //Receive
		char* data = new char[tsize];
		MPI_Recv(data,tsize,MPI_CHAR,procId,tag,MPI_COMM_WORLD,&stat);

		Serializer<stringstream> ser;
		ser.write(data,tsize);
		int vsize,hsize;
		ser>>vsize>>hsize;
		if(procId == 1)
			this->lsize += hsize;

		if(ac.size()<vsize*procId)
		{
			//This means first gather
			for(int i = 0;i<vsize;i++)
			{
				double* a = new double[hsize];
				for(int j = 0;j<hsize;j++)
					ser>>a[j];
				ac.push_back(a);
			}
		}
		else
		{
			//This means data exists so add by extending
			for(int i = 0;i<vsize;i++)
			{
				double* a = new double[lsize];
				memcpy(a,ac[(procId-1)*vsize + i],sizeof(double)*(lsize-hsize));

				for(int j = 0;j<hsize;j++)
					ser>>a[lsize-hsize+j];

				delete[] ac[(procId-1)*vsize + i];
				ac[(procId-1)*vsize + i] = a;
			}
		}

		delete[] data;
	}

	return SUCCESS;
}


template <class T,class U>
void AutoCorrFull<T,U>::serialize(Serializer<U>& obj)
{
	obj<<Wcollection;
}

template <class T,class U>
void AutoCorrFull<T,U>::unserialize(Serializer<U>& obj)
{
	Wcollection.clear();
	obj>>Wcollection;
}
///////////////////////////////////////////////////////////////////////////

template void AutoCorrFull<int,stringstream>::measure();
template void AutoCorrFull<int,stringstream>::writeViaIndex(int idx);
template void AutoCorrFull<int,stringstream>::clear();
template void AutoCorrFull<int,stringstream>::gather(void* p);
template Observable<int,stringstream>* AutoCorrFull<int,stringstream>::duplicate(core::WalkerState<int,stringstream>&);
template void AutoCorrFull<int,stringstream>::copy(void* p);
template int AutoCorrFull<int,stringstream>::parallelSend();
template int AutoCorrFull<int,stringstream>::parallelReceive();
template void AutoCorrFull<int,stringstream>::serialize(Serializer<stringstream>&);
template void AutoCorrFull<int,stringstream>::unserialize(Serializer<stringstream>&);

template void AutoCorrFull<float,stringstream>::measure();
template void AutoCorrFull<float,stringstream>::writeViaIndex(int idx);
template void AutoCorrFull<float,stringstream>::clear();
template void AutoCorrFull<float,stringstream>::gather(void* p);
template Observable<float,stringstream>* AutoCorrFull<float,stringstream>::duplicate(core::WalkerState<float,stringstream>&);
template void AutoCorrFull<float,stringstream>::copy(void* p);
template int AutoCorrFull<float,stringstream>::parallelSend();
template int AutoCorrFull<float,stringstream>::parallelReceive();
template void AutoCorrFull<float,stringstream>::serialize(Serializer<stringstream>&);
template void AutoCorrFull<float,stringstream>::unserialize(Serializer<stringstream>&);
} /* namespace measures */


#endif /* AutoCorrFull_CPP_ */
