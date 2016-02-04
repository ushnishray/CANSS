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

template <class T>
class BasicObs: public measures::Observable<T>, public measures::MPIObservable {
public:
	double dt;
	vect<double> Q;
	int Zcount;
	long ltime;
	Weight& freeEnergy;
	Weight Qx,Qy,Qz;
	Weight Q2x,Q2y,Q2z;

	BasicObs(core::WalkerState<T>& _state, string bsf, FILE* log, double _dt) : Observable<T>(_state,bsf,log),freeEnergy(*(new Weight(_state.weight)))
	{
		dt = _dt;
		Zcount = 0;
		ltime = 0;
		Q.x = 0.0;
		Q.y = 0.0;
		Q.z = 0.0;
	}

	BasicObs(int pId,int nprocs, core::WalkerState<T>& _state, string bsf, FILE* log, double _dt) : MPIObservable(pId,nprocs),Observable<T>(_state,bsf,log),
			freeEnergy(*(new Weight(_state.weight)))
	{
		dt = _dt;
		Zcount = 0;
		ltime = 0;
		Q.x = 0.0;
		Q.y = 0.0;
		Q.z = 0.0;
	}

	~BasicObs()
	{
		delete &freeEnergy;
	}

	void measure();
	void writeViaIndex(int idx);
	void gather(void*);
	void clear();
	Observable<T>* duplicate(core::WalkerState<T>&);

	////////////////////////////////////////////////////////////////////////////////////////////
	//MPI Comm
	////////////////////////////////////////////////////////////////////////////////////////////

	int parallelSend(); //To be called by slaves
	int parallelReceive(); //To be called by master

};
} /* namespace measures */
#endif /* BasicObs_H_ */
