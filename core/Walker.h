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
	template <class T, typename U>
	class Walker:public Serializable<U>
	{
	public:

		gsl_rng* rgenref;

		WalkerState<T,U>& state;
		vector<measures::Observable<T,U>*>& observablesCollection;

		Walker(gsl_rng* r, WalkerState<T,U>& _state, vector<measures::Observable<T,U>*>& _obc):
			rgenref(r),state(_state),observablesCollection(_obc)
		{ }

		//Copy Constructor
		Walker(Walker& w):
			state(*(w.state.duplicate())),observablesCollection(*(new vector<measures::Observable<T,U>*>))
		{
			rgenref = w.rgenref;

			for(int i=0;i<w.observablesCollection.size();i++)
				observablesCollection.push_back(w.observablesCollection[i]->duplicate());
		}

		~Walker()
		{
			delete &state;
			observablesCollection.clear();
			delete &observablesCollection;
		}
		
		void display();

		void reset()
		{
			state.reset();
		}

		void clear()
		{
			for(int i = 0;i<observablesCollection.size();i++)
				observablesCollection[i]->clear();
		}
		
		bool operator() (Walker& a, Walker& b)
		{
			return (a.state.weight.logValue()<b.state.weight.logValue());
		}

		void* getWalkerState()
		{
			return (&state);
		}

		virtual Walker* duplicate();
		virtual void measure();
		virtual void copy(Walker& w);

		template<typename X>
		friend class Serializer;

		virtual void serialize(Serializer<U>&);

		virtual void unserialize(Serializer<U>&);
	};
}

#endif /* WALKER_H_ */
