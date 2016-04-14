/*
 * CloneMultiplicity.cpp
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
void CloneMultiplicity<T,U>::measure() {

}

template <class T, class U>
void CloneMultiplicity<T,U>::writeViaIndex(int idx) {
#if DEBUG >= 2
	fprintf(this->log,"CloneMultiplicity Write\n");
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
	double Qxe = sqrt((Q2.x - Q.x*Q.x)*iZ);
	double Qye = sqrt((Q2.y - Q.y*Q.y)*iZ);
	double Qze = sqrt((Q2.z - Q.z*Q.z)*iZ);

	V.x *= iZ;
	V.y *= iZ;
	V.z *= iZ;
	V2.x *= iZ;
	V2.y *= iZ;
	V2.z *= iZ;
	double Vxe = sqrt((V2.x - V.x*V.x)*iZ);
	double Vye = sqrt((V2.y - V.y*V.y)*iZ);
	double Vze = sqrt((V2.z - V.z*V.z)*iZ);

	double afE = fe*iZ;
	double afE2 = sqrt((this->fe2*iZ - afE*afE)*iZ);

	wif<<t*Zcount<<" "<<afE<<" "<<afE2;

	wif<<" "<<setfill(' ')<<Q.x<<" "<<Qxe<<" "<<setfill(' ')<<Q.y<<" "<<Qye<<" "<<setfill(' ')<<Q.z<<" "<<Qze
			<<" "<<V.x<<" "<<Vxe<<" "<<V.y<<" "<<Vye<<" "<<V.z<<" "<<Vze;

#ifndef NOBRANCH
	Qa.x *= iZ;
	Qa.y *= iZ;
	Qa.z *= iZ;
	Qa2.x *= iZ;
	Qa2.y *= iZ;
	Qa2.z *= iZ;
	Qxe = sqrt((Qa2.x - Qa.x*Qa.x)*iZ);
	Qye = sqrt((Qa2.y - Qa.y*Qa.y)*iZ);
	Qze = sqrt((Qa2.z - Qa.z*Qa.z)*iZ);

	Va.x *= iZ;
	Va.y *= iZ;
	Va.z *= iZ;
	Va2.x *= iZ;
	Va2.y *= iZ;
	Va2.z *= iZ;
	Vxe = sqrt((Va2.x - Va.x*Va.x)*iZ);
	Vye = sqrt((Va2.y - Va.y*Va.y)*iZ);
	Vze = sqrt((Va2.z - Va.z*Va.z)*iZ);

	afE = fea*iZ;
	afE2 = sqrt((this->fea2*iZ - afE*afE)*iZ);

	wif<<" # "<<afE<<" "<<afE2;
	wif<<" "<<setfill(' ')<<Qa.x<<" "<<Qxe<<" "<<setfill(' ')<<Qa.y<<" "<<Qye<<" "<<setfill(' ')<<Qa.z<<" "<<Qze
			<<" "<<Va.x<<" "<<Vxe<<" "<<Va.y<<" "<<Vye<<" "<<Va.z<<" "<<Vze<<endl;
#else
	wif<<endl;
#endif
	wif.close();

	clear();
}

template <class T, class U>
void CloneMultiplicity<T,U>::clear()
{

}

template <class T, class U>
void CloneMultiplicity<T,U>::branchGather(void* p)
{
#ifndef NOBRANCH
	CloneMultiplicity<T,U>* obj = (CloneMultiplicity<T,U>*)p;

#endif
}

template <class T, class U>
void CloneMultiplicity<T,U>::gather(void* p)
{
	//Store into global accumulator
	CloneMultiplicity<T,U>* obj = (CloneMultiplicity<T,U>*)p;

	//Global gather
	obj->ltime = this->state.ltime;

}


template <class T, class U>
void CloneMultiplicity<T,U>::copy(void* p)
{
	CloneMultiplicity<T,U>* obj = (CloneMultiplicity<T,U>*)p;
	this->ltime = obj->ltime;
	this->Q.x = obj->Q.x;
	this->Q.y = obj->Q.y;
	this->Q.z = obj->Q.z;
}

template <class T, class U>
Observable<T,U>* CloneMultiplicity<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	CloneMultiplicity<T,U>* newo = new CloneMultiplicity<T,U>(this->processId,this->procCount,this->totalWalkers,ws,
			this->baseFileName,this->log,this->dt);
	newo->ltime = this->ltime;
	newo->Q.x = this->Q.x;
	newo->Q.y = this->Q.y;
	newo->Q.z = this->Q.z;
	return newo;
}

template <class T, class U>
void CloneMultiplicity<T,U>::display()
{
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"Basic Observable\n");
	fprintf(this->log,"==============================================\n");
	fprintf(this->log,"dt: %f\n",dt);
	fprintf(this->log,"==============================================\n");
	fflush(this->log);
}

template <class T, class U>
int CloneMultiplicity<T,U>::parallelSend()
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
	fprintf(this->log,"CloneMultiplicity Received notification from master\n");
	fflush(this->log);
#endif

	MPI_Send(&this->ltime,1,MPI_UNSIGNED,0,tag,MPI_COMM_WORLD);
	ltime = 0;

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

#ifndef NOBRANCH
	freeEnergya.mpiSend(0);
	freeEnergya.resetValue();

	Qax.mpiSend(0);
	Qax.resetValue();
	Qay.mpiSend(0);
	Qay.resetValue();
	Qaz.mpiSend(0);
	Qaz.resetValue();

	Qax2.mpiSend(0);
	Qax2.resetValue();
	Qay2.mpiSend(0);
	Qay2.resetValue();
	Qaz2.mpiSend(0);
	Qaz2.resetValue();
#endif

#if DEBUG >= 2
	fprintf(this->log,"CloneMultiplicity Transfer Complete\n");
	fflush(this->log);
#endif

	return SUCCESS;
}

template <class T, class U>
int CloneMultiplicity<T,U>::parallelReceive()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);
	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
#if DEBUG >= 2
		fprintf(this->log,"CloneMultiplicity Sending notification to process:%d\n",procId);
		fflush(this->log);
#endif

		MPI_Recv(&this->ltime,1,MPI_UNSIGNED,procId,tag,MPI_COMM_WORLD,&stat);

		this->freeEnergy.mpiReceive(procId);

		this->Qx.mpiReceive(procId);
		this->Qy.mpiReceive(procId);
		this->Qz.mpiReceive(procId);

		this->Qx2.mpiReceive(procId);
		this->Qy2.mpiReceive(procId);
		this->Qz2.mpiReceive(procId);

#ifndef NOBRANCH
		this->freeEnergya.mpiReceive(procId);

		this->Qax.mpiReceive(procId);
		this->Qay.mpiReceive(procId);
		this->Qaz.mpiReceive(procId);

		this->Qax2.mpiReceive(procId);
		this->Qay2.mpiReceive(procId);
		this->Qaz2.mpiReceive(procId);
#endif

#if DEBUG >= 2
		fprintf(this->log,"CloneMultiplicity Finished receiving from process:%d\n",procId);
		fflush(this->log);
#endif
	}

	//Compute observables
	this->Zcount++;
	double it = 1.0/(this->ltime*this->dt);
	double lqx = (Qx/freeEnergy).value();
	double lqy = (Qy/freeEnergy).value();
	double lqz = (Qz/freeEnergy).value();
	this->Q.x += it*lqx;
	this->Q.y += it*lqy;
	this->Q.z += it*lqz;
	this->Q2.x += it*it*lqx*lqx;
	this->Q2.y += it*it*lqy*lqy;
	this->Q2.z += it*it*lqz*lqz;

	double lvx = it*((Qx2/freeEnergy).value() - lqx*lqx);
	double lvy = it*((Qy2/freeEnergy).value() - lqy*lqy);
	double lvz = it*((Qz2/freeEnergy).value() - lqz*lqz);
	this->V.x += lvx;
	this->V.y += lvy;
	this->V.z += lvz;
	this->V2.x += lvx*lvx;
	this->V2.y += lvy*lvy;
	this->V2.z += lvz*lvz;

	double offset = log(this->totalWalkers);
	double temp = it*(freeEnergy.logValue()-offset);
	this->fe += temp;
	this->fe2 += temp*temp;

#ifndef NOBRANCH
	double lqax = (Qax/freeEnergya).value();
	double lqay = (Qay/freeEnergya).value();
	double lqaz = (Qaz/freeEnergya).value();
	this->Qa.x += it*lqax;
	this->Qa.y += it*lqay;
	this->Qa.z += it*lqaz;
	this->Qa2.x += it*it*lqax*lqax;
	this->Qa2.y += it*it*lqay*lqay;
	this->Qa2.z += it*it*lqaz*lqaz;

	double lvax = it*((Qax2/freeEnergya).value() - lqax*lqax);
	double lvay = it*((Qay2/freeEnergya).value() - lqay*lqay);
	double lvaz = it*((Qaz2/freeEnergya).value() - lqaz*lqaz);
	this->Va.x += lvax;
	this->Va.y += lvay;
	this->Va.z += lvaz;
	this->Va2.x += lvax*lvax;
	this->Va2.y += lvay*lvay;
	this->Va2.z += lvaz*lvaz;

	double temp1 = it*(freeEnergya.logValue()-offset);
	this->fea += temp1;
	this->fea2 += temp1*temp1;
#endif

#ifndef NOBRANCH
	lqx*=it;
	lqax*=it;
	ofstream wif(this->baseFileName + "P",std::ofstream::app);
	wif<<it<<" "<<temp<<" "<<temp1<<" "<<lqx<<" "<<lqax<<endl;
	wif.close();
#else
	lqx*=it;
	ofstream wif(this->baseFileName + "P",std::ofstream::app);
	wif<<it<<" "<<temp<<" "<<lqx<<endl;
	wif.close();
#endif

	//Reset for next collection
	freeEnergy.resetValue();
	Qx.resetValue();
	Qy.resetValue();
	Qz.resetValue();
	Qx2.resetValue();
	Qy2.resetValue();
	Qz2.resetValue();
#ifndef NOBRANCH
	freeEnergya.resetValue();
	Qax.resetValue();
	Qay.resetValue();
	Qaz.resetValue();
	Qax2.resetValue();
	Qay2.resetValue();
	Qaz2.resetValue();
#endif

}

template <class T, class U>
void CloneMultiplicity<T,U>::serialize(Serializer<U>& obj)
{
	obj<<dt<<Q;
}

template <class T, class U>
void CloneMultiplicity<T,U>::unserialize(Serializer<U>& obj)
{
	obj>>dt>>Q;
}
///////////////////////////////////////////////////////////////////////////

template void CloneMultiplicity<int,stringstream>::measure();
template void CloneMultiplicity<int,stringstream>::writeViaIndex(int idx);
template void CloneMultiplicity<int,stringstream>::clear();
template void CloneMultiplicity<int,stringstream>::gather(void* p);
template Observable<int,stringstream>* CloneMultiplicity<int,stringstream>::duplicate(core::WalkerState<int,stringstream>&);
template void CloneMultiplicity<int,stringstream>::copy(void* p);
template int CloneMultiplicity<int,stringstream>::parallelSend();
template int CloneMultiplicity<int,stringstream>::parallelReceive();
template void CloneMultiplicity<int,stringstream>::serialize(Serializer<stringstream>&);
template void CloneMultiplicity<int,stringstream>::unserialize(Serializer<stringstream>&);

template void CloneMultiplicity<float,stringstream>::measure();
template void CloneMultiplicity<float,stringstream>::writeViaIndex(int idx);
template void CloneMultiplicity<float,stringstream>::clear();
template void CloneMultiplicity<float,stringstream>::gather(void* p);
template Observable<float,stringstream>* CloneMultiplicity<float,stringstream>::duplicate(core::WalkerState<float,stringstream>&);
template void CloneMultiplicity<float,stringstream>::copy(void* p);
template int CloneMultiplicity<float,stringstream>::parallelSend();
template int CloneMultiplicity<float,stringstream>::parallelReceive();
template void CloneMultiplicity<float,stringstream>::serialize(Serializer<stringstream>&);
template void CloneMultiplicity<float,stringstream>::unserialize(Serializer<stringstream>&);
} /* namespace measures */


