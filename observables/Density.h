/*
 * Density.h
 *
 *  Created on: Aug 20, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.

 This program and the accompanying materials are made available under explicit agreement
 between Ushnish Ray and the end user. You may not redistribute the code and
 accompanying material to anyone.

 On the event that the software is used to generate data that is used implicitly or explicitly
 for research purposes, proper acknowledgment must provided  in the citations section of
 publications.

 This software is cannot be used for commercial purposes in any way whatsoever.
 */

#ifndef Density_H_
#define Density_H_

#include "Observable.h"

namespace measures {

template <class T>
class Density: public measures::Observable<T>, public measures::MPIObservable {
public:

	vectToValue<int> rho;
	int Zcount;

	Density(core::WalkerState<T>& _state, string bsf, FILE* log) : Observable<T>(_state,bsf,log)
	{
		Zcount = 0;
	}

	Density(int pId,int nprocs, core::WalkerState<T>& _state, string bsf, FILE* log) : MPIObservable(pId,nprocs),Observable<T>(_state,bsf,log)
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
	Observable<T>* duplicate();

	////////////////////////////////////////////////////////////////////////////////////////////
	//MPI Comm
	////////////////////////////////////////////////////////////////////////////////////////////

	int parallelSend(); //To be called by slaves
	int parallelReceive(); //To be called by master

};
} /* namespace measures */
#endif /* Density_H_ */
