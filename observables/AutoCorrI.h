/*
 * AutoCorrI.h
 *
 *  Created on: Mar 22, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */

#ifndef AutoCorrI_H_
#define AutoCorrI_H_

#include "Observable.h"

namespace measures{

template <class T, class U>
class AutoCorrI: public measures::Observable<T,U>, public measures::MPIObservable {
public:

	//For local gathering
	vector<Weight> lg;

	//For gathering
	int hsize;
	vector<vector<Weight>*> ac;

	AutoCorrI(core::WalkerState<T,U>& _state, string bsf, FILE* log) : Observable<T,U>(_state,bsf,log)
	{
		hsize = 0;
	}

	AutoCorrI(int pId,int nprocs,int tw, core::WalkerState<T,U>& _state, string bsf, FILE* log) : MPIObservable(pId,nprocs,tw),Observable<T,U>(_state,bsf,log)
	{
		hsize = 0;
	}

	~AutoCorrI()
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


#endif /* AutoCorrI_H_ */
