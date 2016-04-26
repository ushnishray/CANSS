/*
 * RSOSObs.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
 */

#include "dmc.h"
#include "RSOSWalkerState.h"
#include "RSOSObs.h"

namespace measures {

template <class T,class U>
void RSOSObs<T,U>::measure() {
	Q.x += this->state.dQ.x;
	Q.y += this->state.dQ.y;
	Q.z += this->state.dQ.z;
	N += this->state.Rcurr->size();
}

template <class T, class U>
void RSOSObs<T,U>::writeViaIndex(int idx) {
#if DEBUG >= 2
	fprintf(this->log,"RSOSObs Write\n");
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

	H *= iZ;
	H2 *= iZ;
	double He = sqrt((H2-H*H)*iZ);

	gN *= iZ;
	N2 *= iZ;
	double Ne = sqrt((N2-gN*gN)*iZ);

	wif<<t*Zcount<<" "<<afE<<" "<<afE2;

	wif<<" "<<setfill(' ')<<Q.x<<" "<<Qxe<<" "<<setfill(' ')<<Q.y<<" "<<Qye<<" "<<setfill(' ')<<Q.z<<" "<<Qze
			<<" "<<V.x<<" "<<Vxe<<" "<<V.y<<" "<<Vye<<" "<<V.z<<" "<<Vze<<" "<<H<<" "<<He<<" "<<gN<<" "<<Ne;

#ifndef NOBRANCH
	//Average Trajectories
	double offset = log(this->totalWalkers);
	int lsize = cavgQ.size();
	double avgqx = 0.0, avgqx2 = 0.0;
	double avgqxe = 0.0, avgqx2e = 0.0;
	double avgqy = 0.0, avgqy2 = 0.0;
	double avgqye = 0.0, avgqy2e = 0.0;
	double avgqz = 0.0, avgqz2 = 0.0;
	double avgqze = 0.0, avgqz2e = 0.0;

	double avgn = 0.0, avgne = 0.0;
	double avgh = 0.0, avghe = 0.0;

	for(int i = 0;i<lsize;i++)
	{
		double lqx = cavgQ[i].x/totalWalkers;
		double lqx2 = (cavgQ2[i].x/totalWalkers - lqx*lqx)*it;
		avgqx2 += lqx2;
		avgqx2e += lqx2*lqx2;
		lqx *= it;
		avgqx += lqx;
		avgqxe += lqx*lqx;

		double lqy = cavgQ[i].y/totalWalkers;
		double lqy2 = (cavgQ2[i].y/totalWalkers - lqy*lqy)*it;
		avgqy2 += lqy2;
		avgqy2e += lqy2*lqy2;
		lqy *= it;
		avgqy += lqy;
		avgqye += lqy*lqy;

		double lqz = cavgQ[i].z/totalWalkers;
		double lqz2 = (cavgQ2[i].z/totalWalkers - lqz*lqz)*it;
		avgqz2 += lqz2;
		avgqz2e += lqz2*lqz2;
		lqz *= it;
		avgqz += lqz;
		avgqze += lqz*lqz;

		double lqn = avgN[i]/totalWalkers;
		avgn += lqn;
		avgne += lqn*lqn;

		double lqh = avgH[i]/totalWalkers;
		avgh += lqh;
		avghe += lqh*lqh;

	}

	double norm = 1.0/lsize;

	avgqx *= norm;
	avgqxe = sqrt((avgqxe*norm - avgqx*avgqx)/(lsize-1));
	avgqx2 *= norm;
	avgqx2e = sqrt((avgqx2e*norm - avgqx2*avgqx2)/(lsize-1));

	avgqy *= norm;
	avgqye = sqrt((avgqye*norm - avgqy*avgqy)/(lsize-1));
	avgqy2 *= norm;
	avgqy2e = sqrt((avgqy2e*norm - avgqy2*avgqy2)/(lsize-1));

	avgqz *= norm;
	avgqze = sqrt((avgqze*norm - avgqz*avgqz)/(lsize-1));
	avgqz2 *= norm;
	avgqz2e = sqrt((avgqz2e*norm - avgqz2*avgqz2)/(lsize-1));

	avgh  *= norm;
	avghe = sqrt((avghe*norm - avgh*avgh)/(lsize-1));
	avgn  *= norm;
	avgne = sqrt((avgne*norm - avgn*avgn)/(lsize-1));


	wif<<" # "<<avgqx<<" "<<avgqxe
			<<" "<<avgqy<<" "<<avgqye
			<<" "<<avgqz<<" "<<avgqze
			<<" "<<avgqx2<<" "<<avgqx2e
			<<" "<<avgqy2<<" "<<avgqy2e
			<<" "<<avgqz2<<" "<<avgqz2e
			<<" "<<avgh<<" "<<avghe
			<<" "<<avgn<<" "<<avgne
			<<endl;
#else
	wif<<endl;
#endif
	wif.close();

	clear();
}

template <class T, class U>
void RSOSObs<T,U>::clear()
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
	avgH.clear();
	avgH2.clear();
	avgN.clear();
	avgN2.clear();
#endif
}

template <class T, class U>
void RSOSObs<T,U>::branchGather(void* p)
{
#ifndef NOBRANCH
	RSOSObs<T,U>* obj = (RSOSObs<T,U>*)p;
	Q.x = Q.y = Q.z = 0;
	N = 0;

	//Averages
	if(obj->cavgQ.size()<cavgQ.size())
	{
		int ls = cavgQ.size();
		obj->cavgQ.resize(ls);
		obj->cavgQ2.resize(ls);
		obj->avgH.resize(ls);
		obj->avgH2.resize(ls);
		obj->avgN.resize(ls);
		obj->avgN2.resize(ls);
	}

	for(int i = 0;i<cavgQ.size();i++)
	{
		obj->cavgQ[i].x += cavgQ[i].x;
		obj->cavgQ2[i].x += cavgQ[i].x*cavgQ[i].x;
		obj->avgH[i] += avgH[i];
		obj->avgH2[i] += avgH[i]*avgH[i];
		obj->avgN[i] += avgN[i];
		obj->avgN2[i] += avgN[i]*avgN[i];
	}
#endif
}

template <class T, class U>
void RSOSObs<T,U>::gather(void* p)
{
	//Store into global accumulator
	RSOSObs<T,U>* obj = (RSOSObs<T,U>*)p;

	RSOSWalkerState<T,U>& ws = (dynamic_cast<RSOSWalkerState<T,U>&>(this->state));

	//Calculate avg height
	int lth = 0;
	for(int i=0;i<ws.L;i++)
		lth += ws.height[i];
	double avgth = lth*1.0/ws.L;

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

	temp = this->state.weight; temp.multUpdate(avgth);
	obj->eH.add(temp);
	temp = this->state.weight; temp.multUpdate(N);
	obj->eN.add(temp);

#ifdef NOBRANCH
	//Reset walker observables for next gather event
	Q.x = Q.y = Q.z = 0;
	N = 0;
#else
	cavgQ.push_back(Q);
	avgN.push_back(N);
	avgH.push_back(avgth);
#endif
}


template <class T, class U>
void RSOSObs<T,U>::copy(void* p)
{
	RSOSObs<T,U>* obj = (RSOSObs<T,U>*)p;
	this->ltime = obj->ltime;
	this->Q.x = obj->Q.x;
	this->Q.y = obj->Q.y;
	this->Q.z = obj->Q.z;
	this->N = obj->N;
#ifndef NOBRANCH
	this->cavgQ = obj->cavgQ;
	this->avgN = obj->avgN;
	this->avgH = obj->avgH;
#endif
}

template <class T, class U>
Observable<T,U>* RSOSObs<T,U>::duplicate(core::WalkerState<T,U>& ws)
{
	RSOSObs<T,U>* newo = new RSOSObs<T,U>(this->processId,this->procCount,this->totalWalkers,ws,
			this->baseFileName,this->log,this->dt);
	newo->ltime = this->ltime;
	newo->Q.x = this->Q.x;
	newo->Q.y = this->Q.y;
	newo->Q.z = this->Q.z;
	newo->N = this->N;
#ifndef NOBRANCH
	newo->cavgQ = this->cavgQ;
	newo->avgH = this->avgH;
	newo->avgN = this->avgN;
#endif
	return newo;
}

template <class T, class U>
void RSOSObs<T,U>::display()
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
int RSOSObs<T,U>::parallelSend()
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
	fprintf(this->log,"RSOSObs Received notification from master\n");
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

	eH.mpiSend(0);
	eH.resetValue();
	eN.mpiSend(0);
	eN.resetValue();

#ifndef NOBRANCH
	int lsize = this->cavgQ.size();
	MPI_Send(&lsize,1,MPI_INT,0,tag,MPI_COMM_WORLD);
	MPI_Send(this->cavgQ.data(),lsize*3,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(this->cavgQ2.data(),lsize*3,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(this->avgH.data(),lsize,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(this->avgH2.data(),lsize,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	MPI_Send(this->avgN.data(),lsize,MPI_INT,0,tag,MPI_COMM_WORLD);
	MPI_Send(this->avgN2.data(),lsize,MPI_INT,0,tag,MPI_COMM_WORLD);

	cavgQ.clear();
	cavgQ2.clear();
	avgH.clear();
	avgH2.clear();
	avgN.clear();
	avgN2.clear();
#endif


#if DEBUG >= 2
	fprintf(this->log,"RSOSObs Transfer Complete\n");
	fflush(this->log);
#endif

	return SUCCESS;
}

template <class T, class U>
int RSOSObs<T,U>::parallelReceive()
{
	if(!this->MPIEnabled)
		return NOTALLOWED;

#ifndef NOBRANCH
	cavgQ.clear();
	cavgQ2.clear();
#endif

	//Wait for all processes to get here
	MPI_Barrier(MPI_COMM_WORLD);
	for(int procId=1;procId<this->procCount;procId++)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&recv,1,MPI_INT,procId,tag,MPI_COMM_WORLD);
#if DEBUG >= 2
		fprintf(this->log,"RSOSObs Sending notification to process:%d\n",procId);
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

		this->eH.mpiReceive(procId);
		this->eN.mpiReceive(procId);

#ifndef NOBRANCH
		int lsize;
		MPI_Recv(&lsize,1,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		vector<vect<double>> tempcQ(lsize);
		vector<vect<double>> tempcQ2(lsize);
		if(this->cavgQ.size()<lsize)
		{
			cavgQ.resize(lsize);
			cavgQ2.resize(lsize);
			avgH.resize(lsize);
			avgH2.resize(lsize);
			avgN.resize(lsize);
			avgN2.resize(lsize);
		}

		MPI_Recv(tempcQ.data(),lsize*3,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		MPI_Recv(tempcQ2.data(),lsize*3,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);

		vector<double> tempH(lsize), tempH2(lsize);
		vector<int> tempN(lsize), tempN2(lsize);
		MPI_Recv(tempH.data(),lsize,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		MPI_Recv(tempH2.data(),lsize,MPI_DOUBLE,procId,tag,MPI_COMM_WORLD,&stat);
		MPI_Recv(tempN.data(),lsize,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);
		MPI_Recv(tempN2.data(),lsize,MPI_INT,procId,tag,MPI_COMM_WORLD,&stat);

		for(int i = 0;i<lsize;i++)
		{
			cavgQ[i].x += tempcQ[i].x;
			cavgQ2[i].x += tempcQ2[i].x;
			avgH[i] += tempH[i];
			avgH2[i] += tempH2[i];
			avgN[i] += tempN[i];
			avgN2[i] += tempN2[i];
		}


#endif

#if DEBUG >= 2
		fprintf(this->log,"RSOSObs Finished receiving from process:%d\n",procId);
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

	double lqN = (eN/freeEnergy).value()*it;
	double lqH = (eH/freeEnergy).value()*it;
	this->gN += lqN;
	this->N2 += lqN*lqN;
	this->H += lqH;
	this->H2 += lqH*lqH;

#ifndef NOBRANCH
	//Average Trajectory
	int lsize = cavgQ.size();
	double ait = it/lsize;
	double avgqx = 0.0, avgqx2 = 0.0;
	double avgqy = 0.0, avgqy2 = 0.0;
	double avgqz = 0.0, avgqz2 = 0.0;

	double avgh = 0.0, avgh2 = 0.0;
	double avgn = 0.0, avgn2 = 0.0;

#if 1
	ofstream dif(this->baseFileName + "A",std::ofstream::out);
#endif
	for(int i = 0;i<lsize;i++)
	{
		double lqx = cavgQ[i].x/this->totalWalkers;
		avgqx += lqx;
		double lqx2 = (cavgQ2[i].x/totalWalkers - lqx*lqx);
		avgqx2 += lqx2;

		double lqy = cavgQ[i].y/this->totalWalkers;
		avgqy += lqy;
		double lqy2 = (cavgQ2[i].y/totalWalkers - lqy*lqy);
		avgqy2 += lqy2;

		double lqz = cavgQ[i].z/this->totalWalkers;
		avgqx += lqz;
		double lqz2 = (cavgQ2[i].z/totalWalkers - lqz*lqz);
		avgqz2 += lqz2;

		double lh = (avgH[i]/totalWalkers);
		avgh += lh;
		double lh2 = (avgH2[i]/totalWalkers - lh*lh);
		avgh2 += lh2;

		double ln = (avgN[i]/totalWalkers);
		avgn += ln;
		double ln2 = (avgN2[i]/totalWalkers - ln*ln);
		avgn2 += ln2;
#if 1
		dif<<i<<" "<<lqx*it<<" "<<lqx2*it<<" "<<lqy*it<<" "<<lqy2*it<<" "<<lqz*it<<" "<<lqz2*it<<
				" "<<lh*it<<" "<<ln*it<<endl;
#endif
	}
#if 1
	dif<<"--------------------------\n";
	dif.close();
#endif

	avgqx2 *= ait;
	avgqx *= ait;
	avgqy2 *= ait;
	avgqy *= ait;
	avgqz2 *= ait;
	avgqz *= ait;
	avgh2 *= ait;
	avgh *= ait;
	avgn2 *= ait;
	avgn *= ait;

#endif

	lqx*=it;
	lqy*=it;
	lqz*=it;
	ofstream wif(this->baseFileName + "E",std::ofstream::app);
	wif<<it<<" "<<temp<<" "<<lqx<<" "<<lvx<<" "<<lqy<<" "<<lvy<<" "<<lqz<<" "<<lvz<<endl;
	wif.close();
#ifndef NOBRANCH
	ofstream aif(this->baseFileName + "D",std::ofstream::app);
	aif<<(this->ltime*this->dt)*lsize<<" "<<avgqx<<" "<<avgqx2<<" "<<avgqy<<" "<<avgqy2<<" "<<avgqz<<" "<<avgqz2
			<<" "<<avgh<<" "<<avgh2<<" "<<avgn<<" "<<avgn2<<endl;
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
	eH.resetValue();
	eN.resetValue();

#ifndef NOBRANCH
#endif

}

template <class T, class U>
void RSOSObs<T,U>::serialize(Serializer<U>& obj)
{
	obj<<dt<<Q;
#ifndef NOBRANCH
	obj<<cavgQ<<avgH<<avgN;
#endif
}

template <class T, class U>
void RSOSObs<T,U>::unserialize(Serializer<U>& obj)
{
	obj>>dt>>Q;
#ifndef NOBRANCH
	cavgQ.clear();
	obj>>cavgQ>>avgH>>avgN;
#endif
}
///////////////////////////////////////////////////////////////////////////

template void RSOSObs<int,stringstream>::measure();
template void RSOSObs<int,stringstream>::writeViaIndex(int idx);
template void RSOSObs<int,stringstream>::clear();
template void RSOSObs<int,stringstream>::gather(void* p);
template Observable<int,stringstream>* RSOSObs<int,stringstream>::duplicate(core::WalkerState<int,stringstream>&);
template void RSOSObs<int,stringstream>::copy(void* p);
template int RSOSObs<int,stringstream>::parallelSend();
template int RSOSObs<int,stringstream>::parallelReceive();
template void RSOSObs<int,stringstream>::serialize(Serializer<stringstream>&);
template void RSOSObs<int,stringstream>::unserialize(Serializer<stringstream>&);

template void RSOSObs<float,stringstream>::measure();
template void RSOSObs<float,stringstream>::writeViaIndex(int idx);
template void RSOSObs<float,stringstream>::clear();
template void RSOSObs<float,stringstream>::gather(void* p);
template Observable<float,stringstream>* RSOSObs<float,stringstream>::duplicate(core::WalkerState<float,stringstream>&);
template void RSOSObs<float,stringstream>::copy(void* p);
template int RSOSObs<float,stringstream>::parallelSend();
template int RSOSObs<float,stringstream>::parallelReceive();
template void RSOSObs<float,stringstream>::serialize(Serializer<stringstream>&);
template void RSOSObs<float,stringstream>::unserialize(Serializer<stringstream>&);
} /* namespace measures */


