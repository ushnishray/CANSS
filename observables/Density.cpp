/*
 * Density.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
 */

#include "dmc.h"

namespace measures {

template <class T, class U>
void Density<T,U>::display()
{
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"Density Observable\n");
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"Zcount: %d\n",Zcount);
	for(typename vectToValue<T>::iterator itr=rho.begin();itr!=rho.end();++itr)
		fprintf(this->log,"%9.6e %9.6e %9.6e -> %24.16e\n",(float) itr->first.z,(float) itr->first.y,(float) itr->first.x,itr->second);
	fprintf(this->log,"==============================================\n");
}

template <class T, class U>
void Density<T,U>::measure() {
	for(typename PtclMap<T>::iterator it = this->state.Rcurr->begin(); it!=this->state.Rcurr->end();++it)
		rho[it->first] += 1.0;
	Zcount++;
}

template <class T, class U>
void Density<T,U>::writeViaIndex(int idx) {
#if DEBUG >= 2
	fprintf(this->log,"Density Write\n");
	fflush(this->log);
#endif
	stringstream s;
	s<<idx;
	string fname = this->baseFileName + "_" + s.str();

	double izc = 1.0/Zcount;
	ofstream wif(fname);
	wif.precision(FIELDPRECISION);
	wif.width(FIELDWIDTH);
	wif.setf(FIELDFORMAT);

	for(typename vectToValue<T>::iterator itr=rho.begin();itr!=rho.end();++itr)
	{
		wif<<itr->first.z<<" "<<itr->first.y<<" "<<itr->first.x<<" "<<itr->second*izc<<endl;
	}
	wif.close();
	rho.clear();

	Zcount = 0;
}

template <class T, class U>
void Density<T,U>::clear()
{
	rho.clear();
	Zcount = 0;
}

template <class T, class U>
void Density<T,U>::gather(void* p)
{
	Density<T,U>* obj = (Density<T,U>*)p;
	for(typename vectToValue<T>::iterator itr=rho.begin();itr!=rho.end();++itr)
		obj->rho[itr->first] = obj->rho[itr->first] + itr->second;
	obj->Zcount+=Zcount;

	Zcount = 0;
	rho.clear();
}

template <class T, class U>
Observable<T,U>* Density<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	Density<T,U>* newo = new Density<T,U>(this->processId,this->procCount,ws,this->baseFileName,this->log);
	newo->Zcount = this->Zcount;

	for(typename PtclMap<T>::iterator it = this->state.Rcurr->begin(); it!=this->state.Rcurr->end();++it)
		newo->rho[it->first] = it->second;

	return newo;
}

template <class T, class U>
void Density<T,U>::copy(void* p)
{
	Density<T,U>* w = (Density<T,U>*)p;
	this->Zcount = w->Zcount;
	this->rho.clear();
	for(typename PtclMap<T>::iterator it = w->state.Rcurr->begin(); it!=w->state.Rcurr->end();++it)
		this->rho[it->first] = it->second;
}

template <class T, class U>
int Density<T,U>::parallelSend()
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
	fprintf(this->log,"Density Received notification from master\n");
	fflush(this->log);
#endif
	//send Zcount
	MPI_Send(&this->Zcount,1,MPI_INT,0,tag,MPI_COMM_WORLD);

	//Send rho
	int rows = rho.size();
	MPI_Send(&rows,1,MPI_INT,0,tag,MPI_COMM_WORLD);
	double* r_is = new double[3*rows];
	double* values = new double[rows];

	int i = 0;
	for(typename vectToValue<T>::iterator it=rho.begin();it!=rho.end();++it)
	{
		r_is[3*i] = (double) it->first.z;
		r_is[3*i+1] = (double) it->first.y;
		r_is[3*i+2] = (double) it->first.x;
		values[i++] = it->second;
	}
	MPI_Send(r_is,rows*3,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(values,rows,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
#if DEBUG >= 2
	fprintf(this->log,"Density Transfer Complete\n");
	fflush(this->log);
#endif
	rho.clear();
	Zcount = 0;

	delete[] r_is;
	delete[] values;

	return SUCCESS;
}

template <class T, class U>
int Density<T,U>::parallelReceive()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);

	this->Zcount = 0;
	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
#if DEBUG >= 2
		fprintf(this->log,"Density Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif
		int zcountrec;
		MPI_Recv(&zcountrec,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		this->Zcount += zcountrec;

		//Receive rho
		int rows=0;
		MPI_Recv(&rows,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		double* r_is = new double[3*rows];
		double* values = new double[rows];
		MPI_Recv(r_is,rows*3,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		MPI_Recv(values,rows,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);

		for(int i=0;i<rows;i++)
		{
			vect<T> r_i(r_is[3*i],r_is[3*i+1],r_is[3*i+2]);
			this->rho[r_i] = this->rho[r_i] + values[i];
		}
		delete[] r_is;
		delete[] values;
#if DEBUG >= 2
		fprintf(this->log,"Density Finished receiving from process:%d\n",procId);
		fflush(this->log);
#endif
	}
}


template <class T, class U>
void Density<T,U>::serialize(Serializer<U>& obj)
{
	obj<<rho<<Zcount;
}

template <class T, class U>
void Density<T,U>::unserialize(Serializer<U>& obj)
{
	rho.clear();
	obj>>rho>>Zcount;
}
///////////////////////////////////////////////////////////////////////////

template void Density<int,stringstream>::measure();
template void Density<int,stringstream>::writeViaIndex(int idx);
template void Density<int,stringstream>::clear();
template void Density<int,stringstream>::gather(void* p);
template Observable<int,stringstream>* Density<int,stringstream>::duplicate(core::WalkerState<int,stringstream>&);
template void Density<int,stringstream>::copy(void* p);
template int Density<int,stringstream>::parallelSend();
template int Density<int,stringstream>::parallelReceive();
template void Density<int,stringstream>::serialize(Serializer<stringstream>&);
template void Density<int,stringstream>::unserialize(Serializer<stringstream>&);

template void Density<float,stringstream>::measure();
template void Density<float,stringstream>::writeViaIndex(int idx);
template void Density<float,stringstream>::clear();
template void Density<float,stringstream>::gather(void* p);
template Observable<float,stringstream>* Density<float,stringstream>::duplicate(core::WalkerState<float,stringstream>&);
template void Density<float,stringstream>::copy(void* p);
template int Density<float,stringstream>::parallelSend();
template int Density<float,stringstream>::parallelReceive();
template void Density<float,stringstream>::serialize(Serializer<stringstream>&);
template void Density<float,stringstream>::unserialize(Serializer<stringstream>&);
} /* namespace measures */


