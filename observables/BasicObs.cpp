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
//	fprintf(this->log,"FE Check %10.6e %10.6e\n",freeEnergy.logValue(),Q.x);

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

	Qx = Qx/freeEnergy;
	Q2x = Q2x/freeEnergy;
	Qy = Qy/freeEnergy;
	Q2y = Q2y/freeEnergy;
	Qz = Qz/freeEnergy;
	Q2z = Q2z/freeEnergy;

	Q.x *= iZ;
	Q.y *= iZ;
	Q.z *= iZ;

	vect<double> Q2;
	Q2.x = t*(Q2x.value() - Qx.value()*Qx.value());
	Q2.y = t*(Q2y.value() - Qy.value()*Qy.value());
	Q2.z = t*(Q2z.value() - Qz.value()*Qz.value());

	wif<<t<<" "<<(freeEnergy.logValue()*it)<<" "<<setfill(' ')<<Qx.value()<<" "<<setfill(' ')<<Qy.value()<<" "<<setfill(' ')
			<<Qz.value()<<" "<<Q2.x<<" "<<Q2.y<<" "<<Q2.z<<" "<<Q.x<<" "<<Q.y<<" "<<Q.z<<endl;
	wif.close();

	clear();
}

template <class T>
void BasicObs<T>::clear()
{
	Zcount = 0;
	Q.x = 0.0;
	Q.y = 0.0;
	Q.z = 0.0;
	ltime = 0;
	freeEnergy.resetValue();
	Qx.resetValue();
	Qy.resetValue();
	Qz.resetValue();
	Q2x.resetValue();
	Q2y.resetValue();
	Q2z.resetValue();
}

template <class T>
void BasicObs<T>::gather(void* p)
{
	double it = 1.0/(ltime*dt);

	BasicObs<T>* obj = (BasicObs<T>*)p;
	obj->ltime = ltime;
	obj->Zcount += 1;

//	fprintf(this->log,"%d\n", obj->Zcount);
	Qx = freeEnergy;
//	Qx.display(this->log);
	Qx.update(Q.x*it);
//	Qx.display(this->log);
	obj->Qx.add(Qx);
//	obj->Qx.display(this->log);
	Q2x = freeEnergy;
	Q2x.update(Q.x*Q.x*it*it);
	obj->Q2x.add(Q2x);

	Qy = freeEnergy;
	Qy.update(Q.y*it);
	obj->Qy.add(Qy);
	Q2y = freeEnergy;
	Q2y.update(Q.y*Q.y*it*it);
	obj->Q2y.add(Q2y);

	Qz = freeEnergy;
	Qz.update(Q.z*it);
	obj->Qz.add(Qz);
	Q2z = freeEnergy;
	Q2z.update(Q.z*Q.z*it*it);
	obj->Q2z.add(Q2z);

	obj->Q.x += Q.x*it;
	obj->Q.y += Q.y*it;
	obj->Q.z += Q.z*it;

//	fprintf(this->log,"%d\n", obj->Zcount);
//	obj->freeEnergy.display(this->log);
	obj->freeEnergy.add(freeEnergy);
//	freeEnergy.display(this->log);
//	fprintf(this->log,"FE Check %10.6e %10.6e\n",freeEnergy.logValue(),Q.x);

	Zcount = 0;
	Q.x = Q.y = Q.z = 0;
	ltime = 0;
	freeEnergy.resetValue();
	Qx.resetValue();
	Qy.resetValue();
	Qz.resetValue();

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

	newo->freeEnergy = this->freeEnergy;
	newo->Qx = this->Qx;
	newo->Qy = this->Qy;
	newo->Qz = this->Qz;
	newo->Q2x = this->Q2x;
	newo->Q2y = this->Q2y;
	newo->Q2z = this->Q2z;

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

	//FreeEnergy
	this->freeEnergy.mpiSend(0);

	//Large
	this->Qx.mpiSend(0);
	this->Qy.mpiSend(0);
	this->Qz.mpiSend(0);
	this->Q2x.mpiSend(0);
	this->Q2y.mpiSend(0);
	this->Q2z.mpiSend(0);

	fprintf(this->log,"BasicObs Transfer Complete\n");
	fflush(this->log);

	Zcount = 0;
	ltime = 0;
	Q.x = Q.y = Q.z = 0;
	freeEnergy.resetValue();
	Qx.resetValue();
	Q2x.resetValue();

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

		//FreeEnergy
		this->freeEnergy.mpiReceive(procId);

		this->Qx.mpiReceive(procId);
		this->Qy.mpiReceive(procId);
		this->Qz.mpiReceive(procId);
		this->Q2x.mpiReceive(procId);
		this->Q2y.mpiReceive(procId);
		this->Q2z.mpiReceive(procId);

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


