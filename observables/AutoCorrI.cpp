/*
 * AutoCorrI.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */

#ifndef AutoCorrI_CPP_
#define AutoCorrI_CPP_

#include "dmc.h"

namespace measures {

template <class T,class U>
void AutoCorrI<T,U>::display()
{
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"Auto-Correlation Observable\n");
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"==============================================\n");
}

template <class T,class U>
void AutoCorrI<T,U>::measure() {

}

template <class T,class U>
void AutoCorrI<T,U>::writeViaIndex(int idx) {
#if DEBUG >= 1
	fprintf(this->log,"AutoCorrI Write\n");
	fprintf(this->log,"AutoCorrI Size: %d %d\n",ac.size(),this->hsize);
	fflush(this->log);
#endif


	//Compute autocorrelation function
	double* acl = new double[this->hsize];
	double* acl2 = new double[this->hsize];
	for(int i = 0;i<this->hsize;i++)
		acl[i] = acl2[i] = 0.0;

	for(int p=0;p<this->ac.size();p++)
	{
		for(int i = 0;i<this->hsize;i++)
		{
			//Compute local a/c
			double lavg1 = 0.0, lavg2 = 0.0;
			double lval = 0.0;
			for(int j = i;j<this->hsize;j++)
			{
				double v1 = (*this->ac[p])[j-i].logValue();
				double v2 = (*this->ac[p])[j].logValue();

				lval += v1*v2;
				lavg1 += v1;
				lavg2 += v2;
			}

			lavg1 /= (hsize-i);
			lavg2 /= (hsize-i);
			lval /= (hsize-i);

			double t = (lval-lavg1*lavg2);
			//double t = lval;
			acl[i] += t;
			acl2[i] += t*t;
		}
	}

	double norm = 1.0/ac.size();
	for(int i = 0;i<this->hsize;i++)
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

	for(int i=0;i<hsize;i++)
		wif<<i<<"\t"<<acl[i]<<"\t"<<acl2[i]<<endl;
	wif.close();

	//Reset for next bin
	delete[] acl;
	delete[] acl2;
	clear();
}

template <class T,class U>
void AutoCorrI<T,U>::clear()
{
	lg.clear();

	hsize = 0;
	for(int i = 0;i<ac.size();i++)
	{
		ac[i]->clear();
		//delete[] ac[i];
	}
	ac.clear();
}

template <class T,class U>
void AutoCorrI<T,U>::gather(void* p)
{
	AutoCorrI<T,U>* obj = (AutoCorrI<T,U>*)p;
	obj->lg.push_back(this->state.weight);
}

template <class T,class U>
Observable<T,U>* AutoCorrI<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	AutoCorrI<T,U>* newo = new AutoCorrI<T,U>(this->processId,this->procCount,this->totalWalkers,ws,
			this->baseFileName,this->log);
	newo->lg = this->lg;
	return newo;
}

template <class T,class U>
void AutoCorrI<T,U>::copy(void* p)
{
	AutoCorrI<T,U>* obj = (AutoCorrI<T,U>*)p;

	this->lg.clear();
	this->lg = obj->lg;
}

template <class T,class U>
int AutoCorrI<T,U>::parallelSend()
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
	fprintf(this->log,"AutoCorrI Received notification from master\n");
	fflush(this->log);
#endif

	Serializer<U> ser;
	ser<<(int) this->lg.size();
	for(int i = 0;i<this->lg.size();i++)
		ser<<this->lg[i];

	ser.seekg(0,ser.end);
	unsigned int tsize = ser.tellg();
	ser.seekg(0,ser.beg);
	char* data = new char[tsize];
	ser.read(data,tsize); //Read from stream

	//First send size
	MPI_Send(&tsize,1,MPI_UNSIGNED,0,tag,MPI_COMM_WORLD); //Send
	MPI_Send(data,tsize,MPI_CHAR,0,tag,MPI_COMM_WORLD); //Send

	lg.clear();
	return SUCCESS;
}

template <class T,class U>
int AutoCorrI<T,U>::parallelReceive()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	int tag, recv;
	MPI_Status stat;

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);

	this->hsize++;
	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
#if DEBUG >= 2
		fprintf(this->log,"AutoCorrI Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif

		unsigned int tsize;
		MPI_Recv(&tsize,1,MPI_UNSIGNED,procId,tag,MPI_COMM_WORLD,&stat); //Receive
		char* data = new char[tsize];
		MPI_Recv(data,tsize,MPI_CHAR,procId,tag,MPI_COMM_WORLD,&stat);

		Serializer<stringstream> ser;
		ser.write(data,tsize);
		int vsize;
		ser>>vsize;

		if(ac.size()<vsize*procId)
		{
			//This means first gather
			for(int i = 0;i<vsize;i++)
			{
				Weight lw;
				ser>>lw;

				vector<Weight>* ww = new vector<Weight>;
				ww->push_back(lw);
				ac.push_back(ww);
			}
		}
		else
		{
			//This means data exists so just add
			for(int i = 0;i<vsize;i++)
			{
				Weight lw;
				ser>>lw;
				ac[(procId-1)*vsize + i]->push_back(lw);
			}
		}

		delete[] data;
	}

	return SUCCESS;
}


template <class T,class U>
void AutoCorrI<T,U>::serialize(Serializer<U>& obj)
{
	obj<<lg;
}

template <class T,class U>
void AutoCorrI<T,U>::unserialize(Serializer<U>& obj)
{
	lg.clear();
	obj>>lg;
}
///////////////////////////////////////////////////////////////////////////

template void AutoCorrI<int,stringstream>::measure();
template void AutoCorrI<int,stringstream>::writeViaIndex(int idx);
template void AutoCorrI<int,stringstream>::clear();
template void AutoCorrI<int,stringstream>::gather(void* p);
template Observable<int,stringstream>* AutoCorrI<int,stringstream>::duplicate(core::WalkerState<int,stringstream>&);
template void AutoCorrI<int,stringstream>::copy(void* p);
template int AutoCorrI<int,stringstream>::parallelSend();
template int AutoCorrI<int,stringstream>::parallelReceive();
template void AutoCorrI<int,stringstream>::serialize(Serializer<stringstream>&);
template void AutoCorrI<int,stringstream>::unserialize(Serializer<stringstream>&);

template void AutoCorrI<float,stringstream>::measure();
template void AutoCorrI<float,stringstream>::writeViaIndex(int idx);
template void AutoCorrI<float,stringstream>::clear();
template void AutoCorrI<float,stringstream>::gather(void* p);
template Observable<float,stringstream>* AutoCorrI<float,stringstream>::duplicate(core::WalkerState<float,stringstream>&);
template void AutoCorrI<float,stringstream>::copy(void* p);
template int AutoCorrI<float,stringstream>::parallelSend();
template int AutoCorrI<float,stringstream>::parallelReceive();
template void AutoCorrI<float,stringstream>::serialize(Serializer<stringstream>&);
template void AutoCorrI<float,stringstream>::unserialize(Serializer<stringstream>&);
} /* namespace measures */


#endif /* AutoCorrI_CPP_ */
