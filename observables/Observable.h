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

template <class T,typename U>
class Observable: public Serializable<U>
{
protected:
	core::WalkerState<T,U>& state;
	string baseFileName;
	FILE* log;

public:

	Observable(core::WalkerState<T,U>& _state, string bfn, FILE* _log) : state(_state)
	{
		baseFileName = bfn;
		log = _log;
	}

	virtual ~Observable() {}

	virtual void measure() = 0;
	virtual void writeViaIndex(int idx) = 0;
	virtual void clear() = 0;
	virtual void gather(void*) = 0;
	virtual void branchGather(void*) {}

	virtual void display() {}

	virtual void copy(void*) = 0; //Copy observables
	virtual Observable* duplicate(core::WalkerState<T,U>& ws) = 0; //So that walkers can self copy

	template<typename X>
	friend class Serializer;
};

}

#endif /* OBSERVABLE_H_ */
