/*
 * Phistogram.cpp
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
void Phistogram<T,U>::display()
{
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"P-Histogram Observable\n");
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"ltime: %d\n",ltime);
	fprintf(this->log,"[log] local Weight: %16.10e\n",localWeight.logValue());
	for(int i = 0;i<Pcollection.size();i++)
		fprintf(this->log,"%16.10e\n",Pcollection[i]);
	fprintf(this->log,"==============================================\n");
}

template <class T,class U>
void Phistogram<T,U>::measure() {
	ltime = this->state.ltime;
	localWeight = this->state.weight;
//	fprintf(this->log,"%10.6e\n",localWeight.logValue()); 
}

template <class T,class U>
void Phistogram<T,U>::writeViaIndex(int idx) {
#if DEBUG >= 2
	fprintf(this->log,"Phistogram Write\n");
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

	for(int i=0;i<Pcollection.size();i++)
		wif<<Pcollection[i]<<endl;
	wif.close();

	//Reset for next bin
	Pcollection.clear();
	clear();
}

template <class T,class U>
void Phistogram<T,U>::clear()
{
	ltime = 0;
	localWeight.resetValue();
	this->Pcollection.clear();
}

template <class T,class U>
void Phistogram<T,U>::gather(void* p)
{
	Phistogram<T,U>* obj = (Phistogram<T,U>*)p;
	obj->ltime = ltime;
	double it = 1.0/(ltime*dt);
	obj->Pcollection.push_back(localWeight.value());

	localWeight.resetValue();	
	ltime = 0;
}

template <class T,class U>
Observable<T,U>* Phistogram<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	Phistogram<T,U>* newo = new Phistogram<T,U>(this->processId,this->procCount,ws,
			this->baseFileName,this->log,this->dt);
	newo->ltime = this->ltime;
	newo->localWeight.copy(this->localWeight);
	newo->Pcollection = this->Pcollection;
	return newo;
}

template <class T,class U>
void Phistogram<T,U>::copy(void* p)
{
	Phistogram<T,U>* obj = (Phistogram<T,U>*)p;
	this->ltime = obj->ltime;
	this->localWeight.copy(obj->localWeight);
	this->Pcollection.clear();
	this->Pcollection = obj->Pcollection;
}

template <class T,class U>
int Phistogram<T,U>::parallelSend()
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
	fprintf(this->log,"Phistogram Received notification from master\n");
	fflush(this->log);
#endif
	//First send size	
	int lsize = this->Pcollection.size();
	MPI_Send(&lsize,1,MPI_INT,0,tag,MPI_COMM_WORLD);

	//Then wait for notification to send data
	MPI_Recv(&recv,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
	MPI_Send(&this->Pcollection.front(),this->Pcollection.size(),MPI_DOUBLE,0,tag,MPI_COMM_WORLD);

	ltime = 0;
	localWeight.resetValue();
	Pcollection.clear();
	return SUCCESS;
}

template <class T,class U>
int Phistogram<T,U>::parallelReceive()
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
		fprintf(this->log,"Phistogram Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif
		int lsize = 0;		
		MPI_Recv(&lsize,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		tsize += lsize;
		psizes[procId-1] = lsize;
	}

	vector<double> lPcollection;
	lPcollection.resize(tsize); //need to resize first
	double* data = lPcollection.data();

	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
#if DEBUG >= 2
		fprintf(this->log,"Phistogram Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif
		MPI_Recv(data,psizes[procId-1],MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		data += psizes[procId-1];
	}
	
	delete[] psizes;

	//Local gathering done
	//Convert to probabilities
	Weight totalF;
	for(int i = 0;i<tsize;i++)
		totalF.addUpdate(lPcollection[i]);

	ofstream wof(this->baseFileName + "_gather",std::ofstream::app);
	wof<<totalF.logValue()/(tsize*this->dt)<<endl;
	wof.close();

	for(int i = 0;i<tsize;i++)
	{
		//Weight lv(lPcollection[i]);
		//fprintf(this->log,"%10.6e ",lPcollection[i]);
		lPcollection[i] = exp(log(lPcollection[i]) - totalF.logValue());//(lv/totalF).value();
		//fprintf(this->log,"%10.6e %10.6e\n",lPcollection[i],totalF.logValue());
	}
	//fprintf(this->log,"XXXXXXXXXXXXXXXXXXX\n");

	//Now put into global collector
	int osize = Pcollection.size();
	this->Pcollection.resize(osize+tsize);
	data = Pcollection.data() + osize;
	memcpy(data,lPcollection.data(),sizeof(double)*tsize);

	lPcollection.clear();
	return SUCCESS;
}


template <class T,class U>
void Phistogram<T,U>::serialize(Serializer<U>& obj)
{
	obj<<dt<<localWeight<<ltime<<Pcollection;
}

template <class T,class U>
void Phistogram<T,U>::unserialize(Serializer<U>& obj)
{
	Pcollection.clear();
	obj>>dt>>localWeight>>ltime>>Pcollection;
}
///////////////////////////////////////////////////////////////////////////

template void Phistogram<int,stringstream>::measure();
template void Phistogram<int,stringstream>::writeViaIndex(int idx);
template void Phistogram<int,stringstream>::clear();
template void Phistogram<int,stringstream>::gather(void* p);
template Observable<int,stringstream>* Phistogram<int,stringstream>::duplicate(core::WalkerState<int,stringstream>&);
template void Phistogram<int,stringstream>::copy(void* p);
template int Phistogram<int,stringstream>::parallelSend();
template int Phistogram<int,stringstream>::parallelReceive();
template void Phistogram<int,stringstream>::serialize(Serializer<stringstream>&);
template void Phistogram<int,stringstream>::unserialize(Serializer<stringstream>&);

template void Phistogram<float,stringstream>::measure();
template void Phistogram<float,stringstream>::writeViaIndex(int idx);
template void Phistogram<float,stringstream>::clear();
template void Phistogram<float,stringstream>::gather(void* p);
template Observable<float,stringstream>* Phistogram<float,stringstream>::duplicate(core::WalkerState<float,stringstream>&);
template void Phistogram<float,stringstream>::copy(void* p);
template int Phistogram<float,stringstream>::parallelSend();
template int Phistogram<float,stringstream>::parallelReceive();
template void Phistogram<float,stringstream>::serialize(Serializer<stringstream>&);
template void Phistogram<float,stringstream>::unserialize(Serializer<stringstream>&);
} /* namespace measures */


