/*
 * AutoCorr.h
 *
 *  Created on: Mar 22, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */

#ifndef AUTOCORR_H_
#define AUTOCORR_H_

#include "Observable.h"

namespace measures{

template <class T, class U>
class AutoCorr: public measures::Observable<T,U>, public measures::MPIObservable {
public:

	vector<double> Wcollection;

	//For gathering
	int Zcount;
	bool allocated;
	int lsize;
	double* autocorr;
	double* autocorr2;

	AutoCorr(core::WalkerState<T,U>& _state, string bsf, FILE* log) : Observable<T,U>(_state,bsf,log)
	{
		Zcount = 0;
		allocated = false;
		lsize = 0;
	}

	AutoCorr(int pId,int nprocs, core::WalkerState<T,U>& _state, string bsf, FILE* log) : MPIObservable(pId,nprocs),Observable<T,U>(_state,bsf,log)
	{
		Zcount = 0;
		allocated = false;
		lsize = 0;
	}

	~AutoCorr()
	{
		if(allocated)
		{
			lsize = 0;
			delete[] autocorr;
			delete[] autocorr2;
		}

		Wcollection.clear();
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


#endif /* AUTOCORR_H_ */
