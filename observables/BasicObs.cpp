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

#if defined NOBRANCH
	wif<<t<<" "<<(this->freeEnergy.logValue()*it);
#else
	wif<<t<<" "<<(this->freeEnergy.value()*it);
#endif

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
}

template <class T, class U>
void BasicObs<T,U>::gather(void* p)
{
	double it = 1.0/(ltime*dt);

	BasicObs<T,U>* obj = (BasicObs<T,U>*)p;
	obj->ltime = ltime;
	obj->Zcount += 1;

	obj->Q.x += Q.x*it;
	obj->Q.y += Q.y*it;
	obj->Q.z += Q.z*it;

	obj->Q2.x += Q.x*Q.x*it*it;
	obj->Q2.y += Q.y*Q.y*it*it;
	obj->Q2.z += Q.z*Q.z*it*it;

#if defined NOBRANCH
	obj->freeEnergy.add(freeEnergy);
//	fprintf(this->log,"FE: %10.6e %10.6e\n",freeEnergy.logValue(),obj->freeEnergy.logValue());
	freeEnergy.resetValue();
#endif

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

	//Q
	MPI_Send(&this->Q.x,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->Q.y,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->Q.z,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);

	//Q2
	MPI_Send(&this->Q2.x,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->Q2.y,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->Q2.z,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);


#if DEBUG >= 2
	fprintf(this->log,"BasicObs Transfer Complete\n");
	fflush(this->log);
#endif
	Zcount = 0;
	ltime = 0;
	Q.x = Q.y = Q.z = 0;
	Q2.x = Q2.y = Q2.z = 0.0;

#if defined NOBRANCH
	freeEnergy.mpiSend(0);
	freeEnergy.resetValue();
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

#if defined NOBRANCH
		this->freeEnergy.mpiReceive(procId);
		//fprintf(this->log,"FE: %10.6e\n",freeEnergy.logValue());
		//fflush(this->log);
#endif

		//fprintf(this->log,"MPI Check %10.6Le %10.6Le %10.6Le\n",freeEnergy.value(),Qx.value(),Q2x.value());
#if DEBUG >= 2
		fprintf(this->log,"BasicObs Finished receiving from process:%d\n",procId);
		fflush(this->log);
#endif
	}
}

template <class T, class U>
void BasicObs<T,U>::serialize(Serializer<U>& obj)
{
	obj<<dt<<Zcount<<ltime<<Q<<freeEnergy<<Q2;
}

template <class T, class U>
void BasicObs<T,U>::unserialize(Serializer<U>& obj)
{
	obj>>dt>>Zcount>>ltime>>Q>>freeEnergy>>Q2;
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
} /* namespace measures */


