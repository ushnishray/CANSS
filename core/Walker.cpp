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

template <class T>
void Walker<T>::measure()
{
	for(int i=0;i<this->observablesCollection.size();i++)
		(this->observablesCollection)[i]->measure();
}

template <class T>
void Walker<T>::display()
{
	fprintf(state.out,"======================================================\n");
	fprintf(state.out,"Walker Display\n");
	fprintf(state.out,"======================================================\n");
	state.display();
	for(int i=0;i<observablesCollection.size();i++)
		observablesCollection[i]->display();
	fprintf(state.out,"======================================================\n");
}

template <class T>
void Walker<T>::copy(Walker<T>& w)
{
	this->state.copy(w.state);
	//All observables needs to point to the walker state
	for(int i=0;i<w.observablesCollection.size();i++)
		observablesCollection[i]->copy(w.observablesCollection[i]);
}

template <class T>
Walker<T>* Walker<T>::duplicate()
{
	vector<measures::Observable<T>*>* newob = new vector<measures::Observable<T>*>;

	WalkerState<T>* ws = state.duplicate();
	//All observables needs to point to the walker state
	for(int i=0;i<observablesCollection.size();i++)
		newob->push_back(observablesCollection[i]->duplicate(*ws));

	return (new Walker<T>(rgenref,*ws,*newob));
}

//////////////////////////////////////////////////////////////////////
template void Walker<int>::measure();
template Walker<int>* Walker<int>::duplicate();
template void Walker<int>::display();
template void Walker<int>::copy(Walker<int>&);
}
