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

template <class T,class U>
void Qhistogram<T,U>::display()
{
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"Q-Histogram Observable\n");
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"ltime: %d\n",ltime);
	fprintf(this->log,"Q: %16.10e %16.10e %16.10e\n",Q.x,Q.y,Q.z);
	for(vector<vect<double>>::iterator itr=Qcollection.begin();itr!=Qcollection.end();++itr)
		fprintf(this->log,"%16.10e %16.10e %16.10e\n",itr->x,itr->y,itr->z);
	fprintf(this->log,"==============================================\n");
}

template <class T,class U>
void Qhistogram<T,U>::measure() {
	Q.x += this->state.dQ.x;
	Q.y += this->state.dQ.y;
	Q.z += this->state.dQ.z;

	ltime = this->state.ltime;
}

template <class T,class U>
void Qhistogram<T,U>::writeViaIndex(int idx) {
#if DEBUG >= 2
	fprintf(this->log,"Qhistogram Write\n");
	fflush(this->log);
#endif
	stringstream s;
	s<<idx;
	string fname = this->baseFileName + "_" + s.str();
	
	ofstream wif(fname,std::ofstream::app);
	wif.precision(FIELDPRECISION);
	wif.width(FIELDWIDTH);
	wif.setf(FIELDFORMAT);
	wif.fill(' ');

	for(int i=0;i<Qcollection.size();i++)
	{
#ifdef NOBRANCH
		wif<<Qcollection[i].x<<" "<<Qcollection[i].y<<" "<<Qcollection[i].z<<endl;
#else
		wif<<Qcollection[i].x<<" "<<Qcollection[i].y<<" "<<Qcollection[i].z<<"\t#\t"<<
				Qacollection[i].x<<" "<<Qacollection[i].y<<" "<<Qacollection[i].z<<endl;
#endif
	}
	wif.close();

	//Reset for next bin
	clear();
}

template <class T,class U>
void Qhistogram<T,U>::clear()
{
	Q.x = 0.0;
	Q.y = 0.0;
	Q.z = 0.0;
	ltime = 0;
	this->Qcollection.clear();
#ifndef NOBRANCH
	this->Qacollection.clear();
#endif
}

template <class T,class U>
void Qhistogram<T,U>::gather(void* p)
{
	Qhistogram<T,U>* obj = (Qhistogram<T,U>*)p;
	obj->ltime = ltime;
	double it = 1.0/(ltime*dt);
	Q.x *= it;
	Q.y *= it;
	Q.z *= it;
	obj->Qcollection.push_back(Q);
	
#ifdef NOBRANCH
	Q.x = Q.y = Q.z = 0;
	ltime = 0;
#endif
}

template <class T,class U>
void Qhistogram<T,U>::branchGather(void* p)
{
#ifndef NOBRANCH
	Qhistogram<T,U>* obj = (Qhistogram<T,U>*)p;
	Qacollection.push_back(Q);

	int lsize = Qacollection.size();
	int osize = obj->Qacollection.size();
	obj->Qacollection.resize(osize+lsize);
	memcpy(obj->Qacollection.data()+osize,Qacollection.data(),sizeof(double)*lsize*3);

	Q.x = Q.y = Q.z = 0;
	ltime = 0;
#endif
}

template <class T,class U>
Observable<T,U>* Qhistogram<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	Qhistogram<T,U>* newo = new Qhistogram<T,U>(this->processId,this->procCount,this->totalWalkers,ws,
			this->baseFileName,this->log,this->dt);
	newo->ltime = this->ltime;
	newo->Q.x = this->Q.x;
	newo->Q.y = this->Q.y;
	newo->Q.z = this->Q.z;
	newo->Qcollection = this->Qcollection;
#ifndef NOBRANCH
	newo->Qacollection = this->Qacollection;
#endif
	return newo;
}

template <class T,class U>
void Qhistogram<T,U>::copy(void* p)
{
	Qhistogram<T,U>* obj = (Qhistogram<T,U>*) p;

	this->ltime = obj->ltime;
	this->Q.x = obj->Q.x;
	this->Q.y = obj->Q.y;
	this->Q.z = obj->Q.z;

	this->Qcollection.clear();
	this->Qcollection = obj->Qcollection;
#ifndef NOBRANCH
	this->Qacollection.clear();
	this->Qacollection = obj->Qacollection;
#endif
}

template <class T,class U>
int Qhistogram<T,U>::parallelSend()
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
	fprintf(this->log,"Qhistogram Received notification from master\n");
	fflush(this->log);
#endif
	//First send size	
	int lsize = this->Qcollection.size();
	MPI_Send(&lsize,1,MPI_INT,0,tag,MPI_COMM_WORLD);
#ifndef NOBRANCH
	lsize = Qacollection.size();
	MPI_Send(&lsize,1,MPI_INT,0,tag,MPI_COMM_WORLD);
#endif

	//Then wait for notification to send data
	MPI_Recv(&recv,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
	MPI_Send(&this->Qcollection.front(),this->Qcollection.size()*3,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);	
#ifndef NOBRANCH
	MPI_Send(&this->Qacollection.front(),this->Qacollection.size()*3,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
#endif

	ltime = 0;
	Q.x = Q.y = Q.z = 0;
	Qcollection.clear();
#ifndef NOBRANCH
	Qacollection.clear();
#endif
	return SUCCESS;
}

template <class T,class U>
int Qhistogram<T,U>::parallelReceive()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);

	long tsize = 0;
	int* psizes = new int[procCount-1];
#ifndef NOBRANCH
	int* pasizes = new int[procCount-1];
	long tasize = 0;
#endif

	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
#if DEBUG >= 2
		fprintf(this->log,"Qhistogram Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif
		int lsize = 0;		
		MPI_Recv(&lsize,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		tsize += lsize;
		psizes[procId-1] = lsize;

#ifndef NOBRANCH
		MPI_Recv(&lsize,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		tasize += lsize;
		pasizes[procId-1] = lsize;
#endif
	}

	vector<vect<double>> lQcollection;
	lQcollection.resize(tsize); //need to resize first
	vect<double>* data = lQcollection.data();
#ifndef NOBRANCH
	Qacollection.clear();
	Qacollection.resize(tasize); //need to resize first
	vect<double>* data1 = Qacollection.data();
#endif

	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
#if DEBUG >= 2
		fprintf(this->log,"Qhistogram Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif
		MPI_Recv(data,psizes[procId-1]*3,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		data += psizes[procId-1];
#ifndef NOBRANCH
		MPI_Recv(data1,pasizes[procId-1]*3,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		data1 += pasizes[procId-1];
#endif
	}
	
#ifndef NOBRANCH
	delete[] psizes;
	delete[] pasizes;
#endif

	//Local gathering done
	//Now put into global collector
	int osize = Qcollection.size();
	this->Qcollection.resize(osize+tsize);
	data = Qcollection.data() + osize;
	memcpy(data,lQcollection.data(),sizeof(double)*tsize*3);
	lQcollection.clear();

	return SUCCESS;
}

template <class T,class U>
void Qhistogram<T,U>::serialize(Serializer<U>& obj)
{
	obj<<dt<<ltime<<Q;
#ifndef NOBRANCH
	obj<<Qacollection;
#endif
}

template <class T,class U>
void Qhistogram<T,U>::unserialize(Serializer<U>& obj)
{
	Qcollection.clear();
	obj>>dt>>ltime>>Q;
#ifndef NOBRANCH
	Qacollection.clear();
	obj>>Qacollection;
#endif
}
///////////////////////////////////////////////////////////////////////////

template void Qhistogram<int,stringstream>::measure();
template void Qhistogram<int,stringstream>::writeViaIndex(int idx);
template void Qhistogram<int,stringstream>::clear();
template void Qhistogram<int,stringstream>::gather(void* p);
template Observable<int,stringstream>* Qhistogram<int,stringstream>::duplicate(core::WalkerState<int,stringstream>&);
template void Qhistogram<int,stringstream>::copy(void* p);
template int Qhistogram<int,stringstream>::parallelSend();
template int Qhistogram<int,stringstream>::parallelReceive();
template void Qhistogram<int,stringstream>::serialize(Serializer<stringstream>&);
template void Qhistogram<int,stringstream>::unserialize(Serializer<stringstream>&);

template void Qhistogram<float,stringstream>::measure();
template void Qhistogram<float,stringstream>::writeViaIndex(int idx);
template void Qhistogram<float,stringstream>::clear();
template void Qhistogram<float,stringstream>::gather(void* p);
template Observable<float,stringstream>* Qhistogram<float,stringstream>::duplicate(core::WalkerState<float,stringstream>&);
template void Qhistogram<float,stringstream>::copy(void* p);
template int Qhistogram<float,stringstream>::parallelSend();
template int Qhistogram<float,stringstream>::parallelReceive();
template void Qhistogram<float,stringstream>::serialize(Serializer<stringstream>&);
template void Qhistogram<float,stringstream>::unserialize(Serializer<stringstream>&);
} /* namespace measures */


