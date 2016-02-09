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

template <class T>
class Whistogram: public measures::Observable<T>, public measures::MPIObservable {
public:
	double dt;
	Weight& localWeight; 
	long ltime;

	//For gathering
	vector<double> Wcollection;

	Whistogram(core::WalkerState<T>& _state, string bsf, FILE* log, double _dt) : Observable<T>(_state,bsf,log), localWeight(*(new Weight(_state.weight)))
	{
		dt = _dt;
		ltime = 0;
	}

	Whistogram(int pId,int nprocs, core::WalkerState<T>& _state, string bsf, FILE* log, double _dt) : MPIObservable(pId,nprocs),Observable<T>(_state,bsf,log),localWeight(*(new Weight(_state.weight)))
	{
		dt = _dt;
		ltime = 0;
	}

	~Whistogram()
	{
		Wcollection.clear();
		delete &localWeight;	
	}

	void measure();
	void writeViaIndex(int idx);
	void gather(void*);
	void clear();
	Observable<T>* duplicate(core::WalkerState<T>&);
	void copy(void*);
	////////////////////////////////////////////////////////////////////////////////////////////
	//MPI Comm
	////////////////////////////////////////////////////////////////////////////////////////////

	int parallelSend(); //To be called by slaves
	int parallelReceive(); //To be called by master

};
} /* namespace measures */
#endif /* Whistogram_H_ */
