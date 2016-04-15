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
	//Average Trajectories
	double offset = log(this->totalWalkers);
	int lsize = cavgQ.size();
	double avgqx = 0.0, avgqx2 = 0.0;
	double avgqxe = 0.0, avgqx2e = 0.0;
	for(int i = 0;i<lsize;i++)
	{
		double lqx = cavgQ[i].x/totalWalkers;
		double lqx2 = (cavgQ2[i].x/totalWalkers - lqx*lqx)*it;
		avgqx2 += lqx2;
		avgqx2e += lqx2*lqx2;

		lqx *= it;
		avgqx += lqx;
		avgqxe += lqx*lqx;
	}

	double norm = 1.0/lsize;
	avgqx *= norm;
	avgqxe = sqrt((avgqxe*norm - avgqx*avgqx)/(lsize-1));

	avgqx2 *= norm;
	avgqx2e = sqrt((avgqx2e*norm - avgqx2*avgqx2)/(lsize-1));

	wif<<" # "<<t*lsize<<" "<<avgqx<<" "<<avgqxe<<" "<<avgqx2<<" "<<avgqx2e<<endl;
#else
	wif<<endl;
#endif
	wif.close();

	clear();
}

template <class T, class U>
void BasicObs<T,U>::clear()
{
	Zcount = 0;
	Q.x = 0.0; Q2.x = 0.0;
	Q.y = 0.0; Q2.y = 0.0;
	Q.z = 0.0; Q2.z = 0.0;
	ltime = 0;

	freeEnergy.resetValue();
	Qx.resetValue();
	Qy.resetValue();
	Qz.resetValue();
	Qx2.resetValue();
	Qy2.resetValue();
	Qz2.resetValue();
	fe = fe2 = 0.0;

#ifndef NOBRANCH
	cavgQ.clear();
	cavgQ2.clear();
#endif
}

template <class T, class U>
void BasicObs<T,U>::branchGather(void* p)
{
#ifndef NOBRANCH
	BasicObs<T,U>* obj = (BasicObs<T,U>*)p;
	Q.x = Q.y = Q.z = 0;

	//Averages
	if(obj->cavgQ.size()<cavgQ.size())
	{
		obj->cavgQ.resize(cavgQ.size());
		obj->cavgQ2.resize(cavgQ.size());
	}

	for(int i = 0;i<cavgQ.size();i++)
	{
		obj->cavgQ[i].x += cavgQ[i].x;
		obj->cavgQ2[i].x += cavgQ[i].x*cavgQ[i].x;
	}
#endif
}

template <class T, class U>
void BasicObs<T,U>::gather(void* p)
{
	//Store into global accumulator
	BasicObs<T,U>* obj = (BasicObs<T,U>*)p;

	//Global gather
	obj->ltime = this->state.ltime;

	Weight temp = this->state.weight; temp.multUpdate(Q.x);
	obj->Qx.add(temp);
	temp = this->state.weight; temp.multUpdate(Q.y);
	obj->Qy.add(temp);
	temp = this->state.weight; temp.multUpdate(Q.z);
	obj->Qz.add(temp);

	temp = this->state.weight; temp.multUpdate(Q.x*Q.x);
	obj->Qx2.add(temp);
	temp = this->state.weight; temp.multUpdate(Q.y*Q.y);
	obj->Qy2.add(temp);
	temp = this->state.weight; temp.multUpdate(Q.z*Q.z);
	obj->Qz2.add(temp);

	obj->freeEnergy.add(this->state.weight);

#ifdef NOBRANCH
	//Reset walker observables for next gather event
	Q.x = Q.y = Q.z = 0;
#else
	cavgQ.push_back(Q);
#endif
}


template <class T, class U>
void BasicObs<T,U>::copy(void* p)
{
	BasicObs<T,U>* obj = (BasicObs<T,U>*)p;
	this->ltime = obj->ltime;
	this->Q.x = obj->Q.x;
	this->Q.y = obj->Q.y;
	this->Q.z = obj->Q.z;
#ifndef NOBRANCH
	this->cavgQ = obj->cavgQ;
#endif
}

template <class T, class U>
Observable<T,U>* BasicObs<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	BasicObs<T,U>* newo = new BasicObs<T,U>(this->processId,this->procCount,this->totalWalkers,ws,
			this->baseFileName,this->log,this->dt);
	newo->ltime = this->ltime;
	newo->Q.x = this->Q.x;
	newo->Q.y = this->Q.y;
	newo->Q.z = this->Q.z;
#ifndef NOBRANCH
	newo->cavgQ = this->cavgQ;
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
	fprintf(this->log,"Q.x: %9.6e\tQ.y: %9.6e\tQ.z: %9.6e\n",Q.x,Q.y,Q.z);
	fprintf(this->log,"Q2.x: %9.6e\tQ2.y: %9.6e\tQ2.z: %9.6e\n",Q2.x,Q2.y,Q2.z);

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
	int lsize = this->cavgQ.size();
	MPI_Send(&lsize,1,MPI_INT,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->cavgQ.front(),lsize*3,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(&this->cavgQ2.front(),lsize*3,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);

	cavgQ.clear();
	cavgQ2.clear();
#endif


#if DEBUG >= 2
	fprintf(this->log,"BasicObs Transfer Complete\n");
	fflush(this->log);
#endif

	return SUCCESS;
}

template <class T, class U>
int BasicObs<T,U>::parallelReceive()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

	cavgQ.clear();
	cavgQ2.clear();

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);
	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
#if DEBUG >= 2
		fprintf(this->log,"BasicObs Sending notification to process:%d\n",procId);
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
		int lsize;
		MPI_Recv(&lsize,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		vector<vect<double>> tempcQ;
		vector<vect<double>> tempcQ2;
		if(this->cavgQ.size()<lsize)
		{
			cavgQ.resize(lsize);
			cavgQ2.resize(lsize);
			tempcQ.resize(lsize);
			tempcQ2.resize(lsize);
		}
		MPI_Recv(tempcQ.data(),lsize*3,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		MPI_Recv(tempcQ2.data(),lsize*3,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		for(int i = 0;i<lsize;i++)
		{
			cavgQ[i].x += tempcQ[i].x;
			cavgQ2[i].x += tempcQ2[i].x;
		}
#endif

#if DEBUG >= 2
		fprintf(this->log,"BasicObs Finished receiving from process:%d\n",procId);
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
	//Average Trajectory
	int lsize = cavgQ.size();
	double ait = it/lsize;
	double avgqx = 0.0, avgqx2 = 0.0;

#if 1
	ofstream dif(this->baseFileName + "D",std::ofstream::app);
#endif
	for(int i = 0;i<lsize;i++)
	{
		double lqx = cavgQ[i].x/this->totalWalkers;
		avgqx += lqx;

		double lqx2 = (cavgQ2[i].x/totalWalkers - lqx*lqx);
		avgqx2 += lqx2;

#if 1
		dif<<i<<" "<<lqx*it<<" "<<lqx2*it<<endl;
#endif
	}
#if 1
	dif<<"--------------------------\n";
	dif.close();
#endif

	avgqx2 *= ait;
	avgqx *= ait;

#endif

	lqx*=it;
	ofstream wif(this->baseFileName + "P",std::ofstream::app);
	wif<<it<<" "<<temp<<" "<<lqx<<endl;
	wif.close();
#ifndef NOBRANCH
	ofstream aif(this->baseFileName + "A",std::ofstream::app);
	aif<<(this->ltime*this->dt)*lsize<<" "<<avgqx<<" "<<avgqx2<<endl;
	aif.close();
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
#endif

}

template <class T, class U>
void BasicObs<T,U>::serialize(Serializer<U>& obj)
{
	obj<<dt<<Q;
#ifndef NOBRANCH
	obj<<cavgQ;
#endif
}

template <class T, class U>
void BasicObs<T,U>::unserialize(Serializer<U>& obj)
{
	obj>>dt>>Q;
#ifndef NOBRANCH
	cavgQ.clear();
	obj>>cavgQ;
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


