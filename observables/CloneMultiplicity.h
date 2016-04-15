/*
 * CloneMultiplicity.h
 *
 *  Created on: Aug 20, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
 */

#ifndef CloneMultiplicity_H_
#define CloneMultiplicity_H_

#include "Observable.h"

namespace measures {

template <class T, class U>
class CloneMultiplicity: public measures::Observable<T,U>, public measures::MPIObservable {
public:
	vector<set<int,gcmpr<int>>> idc;
	
	CloneMultiplicity(core::WalkerState<T,U>& _state, string bsf, FILE* log) : Observable<T,U>(_state,bsf,log)
	{
	}

	CloneMultiplicity(int pId,int nprocs,int tw, core::WalkerState<T,U>& _state, string bsf, FILE* log) : MPIObservable(pId,nprocs,tw),Observable<T,U>(_state,bsf,log)
	{
	}

	~CloneMultiplicity()
	{
		idc.clear();
	}

	void display();
	void measure();
	void writeViaIndex(int idx);
	void gather(void*);
	void branchGather(void*);
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
#endif /* CloneMultiplicity_H_ */
