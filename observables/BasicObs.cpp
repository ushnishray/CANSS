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

template <class T,class U>
void BasicObs<T,U>::measure() {
	Q.x += this->state.dQ.x;
	Q.y += this->state.dQ.y;
	Q.z += this->state.dQ.z;

	ltime = this->state.ltime;

#if defined NOBRANCH
	this->freeEnergy.copy(this->state.weight);
#endif
}

template <class T, class U>
void BasicObs<T,U>::writeViaIndex(int idx) {
#if DEBUG >= 2
	fprintf(this->log,"BasicObs Write\n");
	fflush(this->log);
#endif
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

	Q2.x = t*(Q2.x-Q.x*Q.x);
	Q2.y = t*(Q2.y-Q.x*Q.y);
	Q2.z = t*(Q2.z-Q.x*Q.z);


	wif<<t<<" "<<(this->freeEnergy.value()*it);

	wif<<" "<<setfill(' ')<<Q.x<<" "<<setfill(' ')<<Q.y<<" "<<setfill(' ')
			<<Q.z<<" "<<Q2.x<<" "<<Q2.y<<" "<<Q2.z<<endl;
	wif.close();

	clear();
}

template <class T, class U>
void BasicObs<T,U>::clear()
{
	Zcount = 0;
	Q.x = 0.0; Q2.x = 0.0;
	Q.y = 0.0; Q2.y = 0.0;
	Q.z = 0.0; Q2.z = 0.0;	ltime = 0;
	freeEnergy.resetValue();

#ifdef NOBRANCH
	Qx.resetValue();
	Qy.resetValue();
	Qz.resetValue();

	Qx2.resetValue();
	Qy2.resetValue();
	Qz2.resetValue();

	lfE.resetValue();
#endif
}

template <class T, class U>
void BasicObs<T,U>::gather(void* p)
{
	//Store into global accumulator
	double it = 1.0/(ltime*dt);
	BasicObs<T,U>* obj = (BasicObs<T,U>*)p;
	obj->ltime = ltime;
	obj->Zcount += 1;

	vect<T> lQ;
	lQ.x = Q.x*it;
	lQ.y = Q.y*it;
	lQ.z = Q.z*it;

	obj->Q.x += lQ.x;
	obj->Q.y += lQ.y;
	obj->Q.z += lQ.z;

	obj->Q2.x += lQ.x*lQ.x;
	obj->Q2.y += lQ.y*lQ.y;
	obj->Q2.z += lQ.z*lQ.z;

#ifdef NOBRANCH
	Weight temp = freeEnergy; temp.multUpdate(lQ.x);
	obj->Qx.add(temp);
	temp = freeEnergy; temp.multUpdate(lQ.y);
	obj->Qy.add(temp);
	temp = freeEnergy; temp.multUpdate(lQ.z);
	obj->Qz.add(temp);

	temp = freeEnergy; temp.multUpdate(lQ.x*lQ.x);
	obj->Qx2.add(temp);
	temp = freeEnergy; temp.multUpdate(lQ.y*lQ.y);
	obj->Qy2.add(temp);
	temp = freeEnergy; temp.multUpdate(lQ.z*lQ.z);
	obj->Qz2.add(temp);

	obj->freeEnergy.add(freeEnergy);
	freeEnergy.resetValue();
#endif

	//Reset walker observables
	Zcount = 0;
	Q.x = Q.y = Q.z = 0;
	ltime = 0;
}

template <class T, class U>
void BasicObs<T,U>::copy(void* p)
{
	BasicObs<T,U>* obj = (BasicObs<T,U>*)p;
	this->ltime = obj->ltime;
	this->Zcount = obj->Zcount;
	this->Q.x = obj->Q.x;
	this->Q.y = obj->Q.y;
	this->Q.z = obj->Q.z;

	this->Q2.x = obj->Q2.x;
	this->Q2.y = obj->Q2.y;
	this->Q2.z = obj->Q2.z;

	this->freeEnergy.copy(obj->freeEnergy);

#ifdef NOBRANCH
	this->lfE = obj->lfE;
	this->Qx = obj->Qx;
	this->Qy = obj->Qy;
	this->Qz = obj->Qz;
	this->Qx2 = obj->Qx;
	this->Qy2 = obj->Qy;
	this->Qz2 = obj->Qz;
#endif
}

template <class T, class U>
Observable<T,U>* BasicObs<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	BasicObs<T,U>* newo = new BasicObs<T,U>(this->processId,this->procCount,ws,
			this->baseFileName,this->log,this->dt);
	newo->ltime = this->ltime;
	newo->Zcount = this->Zcount;
	newo->Q.x = this->Q.x;
	newo->Q.y = this->Q.y;
	newo->Q.z = this->Q.z;

	newo->Q2.x = this->Q2.x;
	newo->Q2.y = this->Q2.y;
	newo->Q2.z = this->Q2.z;

	newo->freeEnergy.copy(this->freeEnergy);

#ifdef NOBRANCH
	newo->lfE = this->lfE;
	newo->Qx = this->Qx;
	newo->Qy = this->Qy;
	newo->Qz = this->Qz;
	newo->Qx2 = this->Qx2;
	newo->Qy2 = this->Qy2;
	newo->Qz2 = this->Qz2;
#endif

	return newo;
}

template <class T, class U>
void BasicObs<T,U>::display()
{
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"Basic Observable\n");
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"dt: %f\n",dt);
	fprintf(this->log,"Zcount: %d\n",Zcount);
	fprintf(this->log,"Bias Free Q.x: %9.6e\tQ.y: %9.6e\tQ.z: %9.6e\n",Q.x,Q.y,Q.z);
	fprintf(this->log,"With Bias Q2.x: %9.6e\tQ2.y: %9.6e\tQ2.z: %9.6e\n",Q2.x,Q2.y,Q2.z);
	fprintf(this->log,"Free Energy: %9.6e\n",freeEnergy.logValue());
	fprintf(this->log,"==============================================\n");
	fflush(this->log);
}

template <class T, class U>
int BasicObs<T,U>::parallelSend()
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
	fprintf(this->log,"BasicObs Received notification from master\n");
	fflush(this->log);
#endif
	//send time
	MPI_Send(&this->ltime,1,MPI_UNSIGNED,0,tag,MPI_COMM_WORLD);

	//send Zcount
	MPI_Send(&this->Zcount,1,MPI_INT,0,tag,MPI_COMM_WORLD);

#ifndef NOBRANCH
	//Q
	MPI_Send(&this->Q.x,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->Q.y,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->Q.z,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);

	//Q2
	MPI_Send(&this->Q2.x,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->Q2.y,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->Q2.z,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
#endif

#if DEBUG >= 2
	fprintf(this->log,"BasicObs Transfer Complete\n");
	fflush(this->log);
#endif
	Zcount = 0;
	ltime = 0;
	Q.x = Q.y = Q.z = 0;
	Q2.x = Q2.y = Q2.z = 0.0;

#ifdef NOBRANCH
	//Free Energy
	freeEnergy.mpiSend(0);
	freeEnergy.resetValue();

	Qx.mpiSend(0);
	Qx.resetValue();
	Qy.mpiSend(0);
	Qy.resetValue();
	Qz.mpiSend(0);
	Qz.resetValue();

	Qx2.mpiSend(0);
	Qx2.resetValue();
	Qy2.mpiSend(0);
	Qy2.resetValue();
	Qz2.mpiSend(0);
	Qz2.resetValue();
#endif

	return SUCCESS;
}

template <class T, class U>
int BasicObs<T,U>::parallelReceive()
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
		fprintf(this->log,"BasicObs Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif
		MPI_Recv(&ltime,1,MPI_UNSIGNED,procId,tag,MPI_COMM_WORLD,&stat);

		int zcountrec;
		MPI_Recv(&zcountrec,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		this->Zcount += zcountrec;

#ifndef NOBRANCH
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
#endif

#ifdef NOBRANCH
		this->lfE.mpiReceive(procId);

		this->Qx.mpiReceive(procId);
		this->Qy.mpiReceive(procId);
		this->Qz.mpiReceive(procId);

		this->Qx2.mpiReceive(procId);
		this->Qy2.mpiReceive(procId);
		this->Qz2.mpiReceive(procId);
#endif
		//fprintf(this->log,"MPI Check %10.6Le %10.6Le %10.6Le\n",freeEnergy.value(),Qx.value(),Q2x.value());
#if DEBUG >= 2
		fprintf(this->log,"BasicObs Finished receiving from process:%d\n",procId);
		fflush(this->log);
#endif
	}

#ifdef NOBRANCH
	//Compute observables
	this->Q.x += (Qx/lfE).value();
	this->Q.y += (Qy/lfE).value();
	this->Q.z += (Qz/lfE).value();

	this->Q2.x += (Qx2/lfE).value();
	this->Q2.y += (Qy2/lfE).value();
	this->Q2.z += (Qz2/lfE).value();

	freeEnergy.addUpdate(lfE.logValue());

	//Reset for next collection
	lfE.resetValue();
	Qx.resetValue();
	Qy.resetValue();
	Qz.resetValue();
	Qx2.resetValue();
	Qy2.resetValue();
	Qz2.resetValue();
#endif
}

template <class T, class U>
void BasicObs<T,U>::serialize(Serializer<U>& obj)
{
	obj<<dt<<Zcount<<ltime<<Q<<freeEnergy<<Q2;
#ifdef NOBRANCH
	obj<<Qx<<Qy<<Qz<<Qx2<<Qy2<<Qz2;
#endif
}

template <class T, class U>
void BasicObs<T,U>::unserialize(Serializer<U>& obj)
{
	obj>>dt>>Zcount>>ltime>>Q>>freeEnergy>>Q2;
#ifdef NOBRANCH
	obj>>Qx>>Qy>>Qz>>Qx2>>Qy2>>Qz2;
#endif
}
///////////////////////////////////////////////////////////////////////////

template void BasicObs<int,stringstream>::measure();
template void BasicObs<int,stringstream>::writeViaIndex(int idx);
template void BasicObs<int,stringstream>::clear();
template void BasicObs<int,stringstream>::gather(void* p);
template Observable<int,stringstream>* BasicObs<int,stringstream>::duplicate(core::WalkerState<int,stringstream>&);
template void BasicObs<int,stringstream>::copy(void* p);
template int BasicObs<int,stringstream>::parallelSend();
template int BasicObs<int,stringstream>::parallelReceive();
template void BasicObs<int,stringstream>::serialize(Serializer<stringstream>&);
template void BasicObs<int,stringstream>::unserialize(Serializer<stringstream>&);

template void BasicObs<float,stringstream>::measure();
template void BasicObs<float,stringstream>::writeViaIndex(int idx);
template void BasicObs<float,stringstream>::clear();
template void BasicObs<float,stringstream>::gather(void* p);
template Observable<float,stringstream>* BasicObs<float,stringstream>::duplicate(core::WalkerState<float,stringstream>&);
template void BasicObs<float,stringstream>::copy(void* p);
template int BasicObs<float,stringstream>::parallelSend();
template int BasicObs<float,stringstream>::parallelReceive();
template void BasicObs<float,stringstream>::serialize(Serializer<stringstream>&);
template void BasicObs<float,stringstream>::unserialize(Serializer<stringstream>&);
} /* namespace measures */


