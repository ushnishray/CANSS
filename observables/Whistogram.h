/*
 * Whistogram.h
 *
 *  Created on: Aug 20, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
 */

#ifndef Whistogram_H_
#define Whistogram_H_

#include "Observable.h"

namespace measures {

template <class T, class U>
class Whistogram: public measures::Observable<T,U>, public measures::MPIObservable {
public:
	double dt;
	Weight& localWeight; 
	unsigned int ltime;

	//For gathering
	vector<double> Wcollection;
	vector<double> Wacollection;

	Whistogram(core::WalkerState<T,U>& _state, string bsf, FILE* log) : Observable<T,U>(_state,bsf,log), localWeight(*(new Weight(_state.weight)))
	{
		dt = 0.0;
		ltime = 0;
	}

	Whistogram(core::WalkerState<T,U>& _state, string bsf, FILE* log, double _dt) : Observable<T,U>(_state,bsf,log), localWeight(*(new Weight(_state.weight)))
	{
		dt = _dt;
		ltime = 0;
	}

	Whistogram(int pId,int nprocs, core::WalkerState<T,U>& _state, string bsf, FILE* log, double _dt) : MPIObservable(pId,nprocs),Observable<T,U>(_state,bsf,log),localWeight(*(new Weight(_state.weight)))
	{
		dt = _dt;
		ltime = 0;
	}

	~Whistogram()
	{
		Wcollection.clear();
		Wacollection.clear();
		delete &localWeight;	
	}

	void measure();
	void writeViaIndex(int idx);
	void gather(void*);
	void branchGather(void*);
	void clear();
	Observable<T,U>* duplicate(core::WalkerState<T,U>&);
	void copy(void*);
	void display();
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
#endif /* Whistogram_H_ */
