/*
 * Walker.cpp
 *
 *  Created on: Aug 22, 2014
 *      Author: ushnish
 *
 	Copyright (c) 2014 Ushnish Ray
	All rights reserved.
 */


#include "dmc.h"

namespace core {

template <class T, class U>
void Walker<T,U>::measure()
{
	for(int i=0;i<this->observablesCollection.size();i++)
		(this->observablesCollection)[i]->measure();
}

template <class T, class U>
void Walker<T,U>::display()
{
	fprintf(state.out,"======================================================\n");
	fprintf(state.out,"Walker Display\n");
	fprintf(state.out,"======================================================\n");
	state.display();
//	for(int i=0;i<observablesCollection.size();i++)
//		observablesCollection[i]->display();
	fprintf(state.out,"======================================================\n");
}

template <class T, class U>
void Walker<T,U>::copy(Walker<T,U>& w)
{
	this->state.copy(w.state);
	//All observables needs to point to the walker state
	for(int i=0;i<w.observablesCollection.size();i++)
		observablesCollection[i]->copy(w.observablesCollection[i]);
}

template <class T, class U>
Walker<T,U>* Walker<T,U>::duplicate()
{
	vector<measures::Observable<T,U>*>* newob = new vector<measures::Observable<T,U>*>;

	WalkerState<T,U>* ws = state.duplicate();
	//All observables needs to point to the walker state
	for(int i=0;i<observablesCollection.size();i++)
		newob->push_back(observablesCollection[i]->duplicate(*ws));

	return (new Walker<T,U>(rgenref,*ws,*newob));
}

template <class T, class U>
void Walker<T,U>::serialize(Serializer<U>& obj)
{
	int l = observablesCollection.size();
	state.serialize(obj);
	obj<<l;

	for(int i=0;i<l;i++)
		observablesCollection[i]->serialize(obj);
}

template <class T, class U>
void Walker<T,U>::unserialize(Serializer<U>& obj)
{
	int l;
	//Replace don't add
	state.unserialize(obj);
	obj>>l;

//	fprintf(this->state.out,"Walker size read: %d\n",l);
//	fflush(this->state.out);

	for(int i = 0;i<l;i++)
		observablesCollection[i]->unserialize(obj);
}
//////////////////////////////////////////////////////////////////////
template void Walker<int,stringstream>::measure();
template void Walker<float,stringstream>::measure();
template Walker<int,stringstream>* Walker<int,stringstream>::duplicate();
template Walker<float,stringstream>* Walker<float,stringstream>::duplicate();
template void Walker<int,stringstream>::display();
template void Walker<float,stringstream>::display();
template void Walker<int,stringstream>::copy(Walker<int,stringstream>&);
template void Walker<float,stringstream>::copy(Walker<float,stringstream>&);
}
