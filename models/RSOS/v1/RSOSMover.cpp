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
#include "RSOSWalkerState.h"
#include "RSOSMover.h"

template <class T, class U>
void RSOSMover<T,U>::initialize(Walker<T,U>* w)
{
	w->state.Rcurr->clear();

	w->state.particleCount = this->rp.L;

	for(int i = 0;i<this->rp.L;i++)
	{
		vect<int> vv(i,1,0);
		(*w->state.Rcurr)[vv] = 1;

		RSOSWalkerState<T,U>& ws = (dynamic_cast<RSOSWalkerState<T,U>&>(w->state));
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
	RSOSWalkerState<T,U>& ws = (dynamic_cast<RSOSWalkerState<T,U>&>(w->state));

	int site = gsl_rng_uniform_int(w->rgenref,rp.L);
	float rd = gsl_rng_uniform(w->rgenref);

	w->state.ltime++;
	//Must do this for zero-move!
	w->state.dweight = 1.0;	
	w->state.dQ.x = 0;

	if(rd<0.5)
	{
		int ht = ws.height[site]+1;
		vect<int> siteloc(site,ht,0);

		int species = (gsl_rng_uniform_int(w->rgenref,2) == 0) ? 1:-1;

		//Attempt to add
		{
			//nearest neighbors
			vect<int> s1((site-1+rp.L)%rp.L,ht,0);
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

			double cost = (-rp.trans.J*species*(s1v + s2v + s3v)) - rp.trans.mu;
			double wt = exp(-cost);
			if(gsl_rng_uniform(w->rgenref) < wt && abs(ws.height[(site-1+rp.L)%rp.L]-(ht))<= 1 && abs(ws.height[(site+1)%rp.L]-(ht))<= 1)
			{
				(*w->state.Rcurr)[siteloc] = species;
				ws.height[site]++;
				w->state.dQ.x = species; //+1/-1 depending on spin; total spin increases by current spin being added
				w->state.dweight = (species == 1) ? rp.trans.rmc : rp.trans.lmc;
				w->state.weight.multUpdate(w->state.dweight);
				w->state.particleCount++;
			}
		}
	}
	else if(ws.height[site] > 1)
	{
		int ht = ws.height[site];
		vect<int> siteloc(site,ht,0);
		//Attempt to remove
		PtclMap<int>::iterator it = w->state.Rcurr->find(siteloc);

		if(it != w->state.Rcurr->end()) //particle found at current location so can remove
		{
			int species = it->second; //Get species

			//nearest neighbors
			vect<int> s1((site-1+rp.L)%rp.L,ht,0);
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

			double cost = rp.trans.J*(s1v*species + s2v*species + s3v*species) + rp.trans.mu;
			double wt = exp(-cost);
			if(gsl_rng_uniform(w->rgenref) < wt  && abs(ws.height[(site-1+rp.L)%rp.L]-(ht-1))<= 1 && abs(ws.height[(site+1)%rp.L]-(ht-1))<= 1)
			{
				w->state.Rcurr->erase(it);
				w->state.particleCount--;
				ws.height[site]--;
				w->state.dQ.x = -species; //Since we are removing total spin reduces by current spin
				w->state.dweight = (species == -1) ? rp.trans.rmc : rp.trans.lmc;
				w->state.weight.multUpdate(w->state.dweight);
			}
		}
	}

//	fprintf(w->state.out,"Moved\n");
//	fflush(w->state.out);
}

/////////////////////////////////////////////////
template void RSOSMover<int,stringstream>::initialize(Walker<int,stringstream>*);
template void RSOSMover<int,stringstream>::move(Walker<int,stringstream>* w);
