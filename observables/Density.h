/*
 * Density.h
 *
 *  Created on: Aug 20, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
 */

#ifndef Density_H_
#define Density_H_

#include "Observable.h"

namespace measures {

template <class T, class U>
class Density: public measures::Observable<T,U>, public measures::MPIObservable {
public:

	vectToValue<int> rho;
	int Zcount;

	Density(core::WalkerState<T,U>& _state, string bsf, FILE* log) : Observable<T,U>(_state,bsf,log)
	{
		Zcount = 0;
	}

	Density(int pId,int nprocs, core::WalkerState<T,U>& _state, string bsf, FILE* log) : MPIObservable(pId,nprocs),Observable<T,U>(_state,bsf,log)
	{
		Zcount = 0;
	}

	~Density()
	{
		rho.clear();
	}

	void measure();
	void writeViaIndex(int idx);
	void gather(void*);
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
#endif /* Density_H_ */
