/*
 * Walker.h
 *
 *  Created on: Aug 21, 2014
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

#ifndef WALKER_H_
#define WALKER_H_

namespace core
{
	template <class T>
	class Walker
	{
	public:

		gsl_rng* rgenref;

		WalkerState<T>& state;
		vector<measures::Observable<T>*>& observablesCollection;

		Walker(gsl_rng* r, WalkerState<T>& _state, vector<measures::Observable<T>*>& _obc):
			rgenref(r),state(_state),observablesCollection(_obc)
		{ }

		//Copy Constructor
		Walker(Walker& w):
			state(*(w.state.duplicate())),observablesCollection(*(new vector<measures::Observable<T>*>))
		{
			rgenref = w.rgenref;

			for(int i=0;i<w.observablesCollection.size();i++)
				observablesCollection.push_back(w.observablesCollection[i]->duplicate());
		}

		~Walker()
		{
			delete &state;
			delete &observablesCollection;
		}
		
		void reset()
		{
			state.reset();		
		}
		
		bool operator() (Walker& a, Walker& b)
		{
			return (a.state.weight.logValue()<b.state.weight.logValue());
		}

		virtual Walker* duplicate();
		virtual void measure();


	};
}

#endif /* WALKER_H_ */
