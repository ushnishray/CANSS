/*
 * observable.h
 *
 *  Created on: Aug 19, 2014
 *      Author: ushnish
 *
 	Copyright (c) 2014 Ushnish Ray
	All rights reserved.
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

	virtual void copy(void*) = 0; //Copy observables
	virtual Observable* duplicate(core::WalkerState<T>& ws) = 0; //So that walkers can self copy

	template<typename U>
	friend class Serializer;
};

}

#endif /* OBSERVABLE_H_ */
