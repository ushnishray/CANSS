/*
 * observable.h
 *
 *  Created on: Aug 19, 2014
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

#ifndef OBSERVABLE_H_
#define OBSERVABLE_H_

namespace measures{

template <class T>
class Observable
{
protected:
	core::WalkerState<T>& state;
	string baseFileName;
	FILE* log;

public:

	Observable(core::WalkerState<T>& _state, string bfn, FILE* _log) : state(_state)
	{
		baseFileName = bfn;
		log = _log;
	}

	virtual ~Observable() {}

	virtual void measure() = 0;
	virtual void writeViaIndex(int idx) = 0;
	virtual void clear() = 0;
	virtual void gather(void*) = 0;

	virtual void display() {}

	virtual Observable* duplicate() = 0; //So that walkers can self copy
};

}

#endif /* OBSERVABLE_H_ */
