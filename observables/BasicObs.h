/*
 * BasicObs.h
 *
 *  Created on: Aug 20, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
 */

#ifndef BasicObs_H_
#define BasicObs_H_

#include "Observable.h"

namespace measures {

template <class T, class U>
class BasicObs: public measures::Observable<T,U>, public measures::MPIObservable {
public:
	double dt;
	vect<double> Q, Q2;
	int Zcount;
	unsigned int ltime;
	Weight& freeEnergy;

#ifdef NOBRANCH
	//Just needed for collection
	Weight Qx,Qy,Qz;
	Weight Qx2,Qy2,Qz2;
	Weight lfE;
#endif

	BasicObs(core::WalkerState<T,U>& _state, string bsf, FILE* log) : Observable<T,U>(_state,bsf,log),freeEnergy(*(new Weight(_state.weight)))
	{
		dt = 0.0;
		Zcount = 0;
		ltime = 0;
		Q.x = 0.0; Q2.x = 0.0;
		Q.y = 0.0; Q2.y = 0.0;
		Q.z = 0.0; Q2.z = 0.0;
	}

	BasicObs(core::WalkerState<T,U>& _state, string bsf, FILE* log, double _dt) : Observable<T,U>(_state,bsf,log),freeEnergy(*(new Weight(_state.weight)))
	{
		dt = _dt;
		Zcount = 0;
		ltime = 0;
		Q.x = 0.0; Q2.x = 0.0;
		Q.y = 0.0; Q2.y = 0.0;
		Q.z = 0.0; Q2.z = 0.0;
	}

	BasicObs(int pId,int nprocs, core::WalkerState<T,U>& _state, string bsf, FILE* log, double _dt) : MPIObservable(pId,nprocs),Observable<T,U>(_state,bsf,log),
			freeEnergy(*(new Weight(_state.weight)))
	{
		dt = _dt;
		Zcount = 0;
		ltime = 0;
		Q.x = 0.0; Q2.x = 0.0;
		Q.y = 0.0; Q2.y = 0.0;
		Q.z = 0.0; Q2.z = 0.0;
	}

	~BasicObs()
	{
		delete &freeEnergy;
	}

	void display();
	void measure();
	void writeViaIndex(int idx);
	void gather(void*);
	void clear();
	Observable<T,U>* duplicate(core::WalkerState<T,U>&);
	void copy(void*);
	////////////////////////////////////////////////////////////////////////////////////////////
	//MPI Comm
	////////////////////////////////////////////////////////////////////////////////////////////

	int parallelSend(); //To be called by slaves
	int parallelReceive(); //To be called by master

	////////////////////////////////////////////////////////////////////////////////////////////
	//Serialization
	////////////////////////////////////////////////////////////////////////////////////////////
	virtual void serialize(Serializer<U>&);
	virtual void unserialize(Serializer<U>&);

};
} /* namespace measures */
#endif /* BasicObs_H_ */
