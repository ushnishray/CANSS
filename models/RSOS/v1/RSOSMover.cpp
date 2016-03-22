/*
 * BHMain.cpp
 *
 *  Created on: Aug 26, 2015
 *      Author: ushnish
 *
 	Copyright (c) 2014 Ushnish Ray
	All rights reserved.
 */

#include "dmc.h"
#include "RunParameters.h"
#include "RSOSMover.h"

template <class T, class U>
void RSOSMover<T,U>::initialize(Walker<T,U>* w)
{
	w->state.Rcurr->clear();

	w->state.particleCount = this->rp.L;

	for(int i = 0;i<this->rp.L;i++)
	{
		vect<float> vv(i,0,0);
		(*w->state.Rcurr)[vv] = 1;

		RSOSWalkerState<T,U>& ws = *(RSOSWalkerState<T,U>*(w->getWalkerState()));
		ws.height[i] = 1;
	}
	w->state.ltime = 0;
	w->state.reset();

//	fprintf(this->debugLog,"Initialized.\n");
//	fflush(this->debugLog);
}

template <class T, class U>
void RSOSMover<T,U>::move(Walker<T,U>* w)
{
	RSOSWalkerState<T,U>& ws = *(RSOSWalkerState<T,U>*(w->getWalkerState()));

	int site = gsl_rng_uniform_int(w->rgenref,rp.L);
	float rd = gsl_rng_uniform(w->rgenref);
	char species = (gsl_rng_uniform_int(w->regenref,2) == 0) ? 1:-1;

	int ht = ws.height[site];
	vect<int> siteloc(site,ht,0);

	if(rd<0.5)
	{
		//Attempt to add
		PtclMap<int>::iterator it = w->state.Rcurr->find(siteloc);

		if(it != w->state.Rcurr->end()) //no particle at current location so can add
		{
			//nearest neighbors
			vect<int> s1((site-1)%rp.L,ht,0);
			PtclMap<int>::iterator it1 = w->state.Rcurr->find(s1);
			char s1v = 0;
			if(it1 != w->state.Rcurr->end())
				s1v = it1->second;

			vect<int> s2((site+1)%rp.L,ht,0);
			it1 = w->state.Rcurr->find(s2);
			char s2v = 0;
			if(it1 != w->state.Rcurr->end())
				s2v = it1->second;

			vect<int> s3(site,ht-1,0);
			it1 = w->state.Rcurr->find(s3);
			char s3v = 0;
			if(it1 != w->state.Rcurr->end())
				s3v = it1->second;

			double cost = (-rp.trans.J*(s1v*species + s2v*species + s3v*species) + rp.trans.mu);
			double wt = exp(cost);
			if(gsl_rng_uniform(w->rgenref) < wt && abs(ws.height[(site-1)%rp.L]-(ht+1))<= 1 && abs(ws.height[(site+1)%rp.L]-(ht+1))<= 1)
			{
				w->state.Rcurr[siteloc] = species;
				ws.height[site]++;
				w->state.dQ.x = cost;
				w->state.dweight = wt;
				w->state.weight.multUpdate(w->state.dweight);
				w->state.particleCount++;
			}
		}
	}
	else
	{
		//Attempt to remove
		PtclMap<int>::iterator it = w->state.Rcurr->find(siteloc);

		if(it == w->state.Rcurr->end()) //particle found at current location so can remove
		{
			//nearest neighbors
			vect<int> s1((site-1)%rp.L,ht,0);
			PtclMap<int>::iterator it1 = w->state.Rcurr->find(s1);
			char s1v = 0;
			if(it1 != w->state.Rcurr->end())
				s1v = it1->second;

			vect<int> s2((site+1)%rp.L,ht,0);
			it1 = w->state.Rcurr->find(s2);
			char s2v = 0;
			if(it1 != w->state.Rcurr->end())
				s2v = it1->second;

			vect<int> s3(site,ht-1,0);
			it1 = w->state.Rcurr->find(s3);
			char s3v = 0;
			if(it1 != w->state.Rcurr->end())
				s3v = it1->second;

			double cost = -(-rp.trans.J*(s1v*species + s2v*species + s3v*species) + rp.trans.mu);
			double wt = exp(cost);
			if(gsl_rng_uniform(w->rgenref) < wt  && abs(ws.height[(site-1)%rp.L]-(ht-1))<= 1 && abs(ws.height[(site+1)%rp.L]-(ht-1))<= 1)
			{
				w->state.Rcurr->erase(it);
				w->state.particleCount--;
				ws.height[site]--;
				w->state.dQ.x = cost;
				w->state.dweight = wt;
				w->state.weight.multUpdate(w->state.dweight);
			}
		}
	}

}

/////////////////////////////////////////////////
template void RSOSMover<int,stringstream>::initialize(Walker<int,stringstream>*);
template void RSOSMover<int,stringstream>::move(Walker<int,stringstream>* w);
