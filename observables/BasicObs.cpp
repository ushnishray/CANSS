/*
 * BasicObs.cpp
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
void BasicObs<T>::measure() {
	Q.x += this->state.dQ.x;
	Q.y += this->state.dQ.y;
	Q.z += this->state.dQ.z;

	ltime = this->state.ltime;
	freeEnergy = this->state.weight;
}

template <class T>
void BasicObs<T>::writeViaIndex(int idx) {

	fprintf(this->log,"BasicObs Write\n");
	fflush(this->log);

	double t = this->ltime*dt;
	double it = 1.0/t;
	double iZ = 1.0/Zcount;

	ofstream wif(this->baseFileName,std::ofstream::app);
	wif.precision(FIELDPRECISION);
	wif.width(FIELDWIDTH);
	wif.setf(FIELDFORMAT);
	wif.fill(' ');

	Q.x *= iZ;
	Q.y *= iZ;
	Q.z *= iZ;
	Q2.x *= iZ;
	Q2.y *= iZ;
	Q2.z *= iZ;
	Q2.x -= (Q.x*Q.x);
	Q2.y -= (Q.y*Q.y);
	Q2.z -= (Q.z*Q.z);

	wif<<t<<" "<<(freeEnergy.logValue()*it)<<" "<<setfill(' ')<<Q.x<<" "<<setfill(' ')<<Q.y<<" "<<setfill(' ')
			<<Q.z<<" "<<Q2.x*t<<" "<<Q2.y*t<<" "<<Q2.z*t<<endl;
	wif.close();

	clear();
}

template <class T>
void BasicObs<T>::clear()
{
	Zcount = 0;
	Q2.x = Q.x = 0.0;
	Q2.y = Q.y = 0.0;
	Q2.z = Q.z = 0.0;
	ltime = 0;
	freeEnergy.resetValue();
}

template <class T>
void BasicObs<T>::gather(void* p)
{
	BasicObs<T>* obj = (BasicObs<T>*)p;
	obj->ltime = ltime;
	obj->Zcount += 1;

	double it = 1.0/(ltime*dt);
	obj->Q.x += Q.x*it;
	obj->Q.y += Q.y*it;
	obj->Q.z += Q.z*it;
	obj->Q2.x += Q.x*Q.x*it*it;
	obj->Q2.y += Q.y*Q.y*it*it;
	obj->Q2.z += Q.z*Q.z*it*it;
	
//	obj->freeEnergy.display(this->log);
	obj->freeEnergy.add(freeEnergy);
//	freeEnergy.display(this->log);
//	fprintf(this->log,"FE Check %10.6e %10.6e\n",obj->freeEnergy.logValue(),freeEnergy.logValue());

	Zcount = 0;
	Q.x = Q.y = Q.z = 0;
	Q2.x = Q2.y = Q2.z = 0;
	ltime = 0;
	freeEnergy.resetValue();

	fflush(this->log);
}

template <class T>
Observable<T>* BasicObs<T>::duplicate(core::WalkerState<T>& ws)
{
	BasicObs<T>* newo = new BasicObs<T>(this->processId,this->procCount,ws,
			this->baseFileName,this->log,this->dt);
	newo->ltime = this->ltime;
	newo->Zcount = this->Zcount;
	newo->Q.x = this->Q.x;
	newo->Q.y = this->Q.y;
	newo->Q.z = this->Q.z;
	newo->Q2.x = this->Q2.x;
	newo->Q2.y = this->Q2.y;
	newo->Q2.z = this->Q2.z;
	newo->freeEnergy = this->freeEnergy;
	return newo;
}


template <class T>
int BasicObs<T>::parallelSend()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	int tag, recv;
	MPI_Status stat;

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);

	//Wait to be notified by master
	MPI_Recv(&recv,1,MPI_INT,0,tag,MPI_COMM_WORLD,&stat);
	fprintf(this->log,"BasicObs Received notification from master\n");
	fflush(this->log);

	//send time
	MPI_Send(&this->ltime,1,MPI_LONG,0,tag,MPI_COMM_WORLD);

	//send Zcount
	MPI_Send(&this->Zcount,1,MPI_INT,0,tag,MPI_COMM_WORLD);

	//Q
	MPI_Send(&this->Q.x,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->Q.y,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->Q.z,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);

	//Q2
	MPI_Send(&this->Q2.x,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->Q2.y,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->Q2.z,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);

	//FreeEnergy
	this->freeEnergy.mpiSend(0);

	fprintf(this->log,"BasicObs Transfer Complete\n");
	fflush(this->log);

	Zcount = 0;
	ltime = 0;
	Q.x = Q.y = Q.z = 0;
	Q2.x = Q2.y = Q2.z = 0;
	freeEnergy.resetValue();

	return SUCCESS;
}

template <class T>
int BasicObs<T>::parallelReceive()
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
		fprintf(this->log,"BasicObs Sending notification to process:%d\n",procId);
		fflush(this->log);

		MPI_Recv(&ltime,1,MPI_LONG,procId,tag,MPI_COMM_WORLD,&stat);

		int zcountrec;
		MPI_Recv(&zcountrec,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		this->Zcount += zcountrec;

		//Receive Q
		double temp=0;
		MPI_Recv(&temp,1,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		this->Q.x += temp;
		MPI_Recv(&temp,1,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		this->Q.y += temp;
		MPI_Recv(&temp,1,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		this->Q.z += temp;

		//Receive Q2
		MPI_Recv(&temp,1,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		this->Q2.x += temp;
		MPI_Recv(&temp,1,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		this->Q2.y += temp;
		MPI_Recv(&temp,1,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		this->Q2.z += temp;

		//FreeEnergy
		this->freeEnergy.mpiReceive(procId);

		fprintf(this->log,"BasicObs Finished receiving from process:%d\n",procId);
		fflush(this->log);
	}
}

///////////////////////////////////////////////////////////////////////////

template void BasicObs<int>::measure();
template void BasicObs<int>::writeViaIndex(int idx);
template void BasicObs<int>::clear();
template void BasicObs<int>::gather(void* p);
template Observable<int>* BasicObs<int>::duplicate(core::WalkerState<int>&);
template int BasicObs<int>::parallelSend();
template int BasicObs<int>::parallelReceive();

} /* namespace measures */


