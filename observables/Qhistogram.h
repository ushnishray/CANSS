/*
 * Qhistogram.h
 *
 *  Created on: Aug 20, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
 */

#ifndef Qhistogram_H_
#define Qhistogram_H_

#include "Observable.h"

namespace measures {

template <class T, class U>
class Qhistogram: public measures::Observable<T,U>, public measures::MPIObservable {
public:
	double dt;
	vect<double> Q;
	unsigned int ltime;

	//For gathering
	vector<vect<double>> Qcollection;

	Qhistogram(core::WalkerState<T,U>& _state, string bsf, FILE* log) : Observable<T,U>(_state,bsf,log)
	{
		dt = 0.0;
		ltime = 0;
		Q.x = 0.0;
		Q.y = 0.0;
		Q.z = 0.0;
	}

	Qhistogram(core::WalkerState<T,U>& _state, string bsf, FILE* log, double _dt) : Observable<T,U>(_state,bsf,log)
	{
		dt = _dt;
		ltime = 0;
		Q.x = 0.0;
		Q.y = 0.0;
		Q.z = 0.0;
	}

	Qhistogram(int pId,int nprocs, core::WalkerState<T,U>& _state, string bsf, FILE* log, double _dt) : MPIObservable(pId,nprocs),Observable<T,U>(_state,bsf,log)
	{
		dt = _dt;
		ltime = 0;
		Q.x = 0.0;
		Q.y = 0.0;
		Q.z = 0.0;
	}

	~Qhistogram()
	{
		Qcollection.clear();
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
#endif /* Qhistogram_H_ */
