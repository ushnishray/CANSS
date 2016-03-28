/*
 * AutoCorrFull.h
 *
 *  Created on: Mar 22, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */

#ifndef AutoCorrFull_H_
#define AutoCorrFull_H_

#include "Observable.h"

namespace measures{

template <class T, class U>
class AutoCorrFull: public measures::Observable<T,U>, public measures::MPIObservable {
public:

	vector<double> Wcollection;

	//For gathering
	int lsize;
	vector<double *> ac;


	AutoCorrFull(core::WalkerState<T,U>& _state, string bsf, FILE* log) : Observable<T,U>(_state,bsf,log)
	{
		lsize = 0;
	}

	AutoCorrFull(int pId,int nprocs, core::WalkerState<T,U>& _state, string bsf, FILE* log) : MPIObservable(pId,nprocs),Observable<T,U>(_state,bsf,log)
	{
		lsize = 0;
	}

	~AutoCorrFull()
	{
		clear();
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

}


#endif /* AutoCorrFull_H_ */
