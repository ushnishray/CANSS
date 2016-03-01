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
#include "SNSMover.h"

template <class T, class U>
void SNSMover<T,U>::initialize(Walker<T,U>* w)
{
	w->state.Rcurr->clear();
	w->state.particleCount = 0;

	int initptclcount = this->rp.L/2*(this->rp.trans.linsert+this->rp.trans.rinsert);
	int pp = 0;
	while(pp<initptclcount)
	{
		vect<int> vv(gsl_rng_uniform_int(w->rgenref,this->rp.L),0,0);
		if(w->state.Rcurr->find(vv)== w->state.Rcurr->end())
		{
			(*w->state.Rcurr)[vv] = pp++;
			w->state.particleCount++;
		}
	}

	//Using the fact that there at most L particles possible
	//for(int i=pp;i<rp.L;i++)
	//	pidxavailable.push_back(i);
	pidx = pp;

	w->state.ltime = 0;
	w->state.reset();
	//w->state.display();
}

template <class T, class U>
void SNSMover<T,U>::move(Walker<T,U>* w)
{
	vect<int> siteloc(gsl_rng_uniform_int(w->rgenref,rp.L+2),0,0);

	double rd = gsl_rng_uniform(w->rgenref);
	int newhop = 0;

	bool ptcl = false;
	PtclMap<int>::iterator it = w->state.Rcurr->find(siteloc);
	if(it!=w->state.Rcurr->end())
		ptcl = true;

	double tweight = 1.0;

	if(siteloc.x == rp.L)
	{
		double num = rp.trans.lremove*rp.trans.lmc;
		tweight = num + rp.trans.linsert*rp.trans.rmc;

		if(rd<num/tweight)
		{
			double ar = gsl_rng_uniform(w->rgenref);
			PtclMap<int>::iterator it = w->state.Rcurr->find(vect<int>(0,0,0));
			if(it!=w->state.Rcurr->end()) //Particle present
			{
				w->state.Rcurr->erase(it);
				w->state.particleCount--;
				newhop--;
			}
		}
		else
		{
			double ar = gsl_rng_uniform(w->rgenref);
			PtclMap<int>::iterator it = w->state.Rcurr->find(vect<int>(0,0,0));
			if(it==w->state.Rcurr->end()) //Particle not present
			{
				(*w->state.Rcurr)[vect<int>(0,0,0)] = pidx++;
				w->state.particleCount++;
				newhop++;
			}
		}
	}
	else if(siteloc.x == rp.L+1)
	{
		double num = rp.trans.rremove*rp.trans.rmc;
		tweight = num + rp.trans.rinsert*rp.trans.lmc;

		if(rd<num/tweight)
		{
			double ar = gsl_rng_uniform(w->rgenref);
			PtclMap<int>::iterator it = w->state.Rcurr->find(vect<int>(rp.L-1,0,0));
			if(it!=w->state.Rcurr->end()) //Particle present
			{
				w->state.Rcurr->erase(it);
				w->state.particleCount--;
				newhop++;
			}
		}
		else
		{
			double ar = gsl_rng_uniform(w->rgenref);
			PtclMap<int>::iterator it = w->state.Rcurr->find(vect<int>(rp.L-1,0,0));
			if(it==w->state.Rcurr->end()) //Particle not present
			{
				(*w->state.Rcurr)[vect<int>(rp.L-1,0,0)] = pidx++;
				w->state.particleCount++;
				newhop--;
			}
		}
	}
	else if(ptcl)
	{
		double num = rp.trans.lmove*rp.trans.lmc;
		tweight = num + rp.trans.rmove*rp.trans.rmc;

		if(rd<num/tweight)
		{
			double ar = gsl_rng_uniform(w->rgenref);
			PtclMap<int>::iterator itl = w->state.Rcurr->find(vect<int>(siteloc.x-1,0,0));
			if(itl==w->state.Rcurr->end() && siteloc.x != 0) //Particle not found at x-1
			{
				newhop--;
				int ptclnum = itl->second;
				w->state.Rcurr->erase(it);
				(*w->state.Rcurr)[vect<int>(siteloc.x-1,0,0)] = ptclnum;
			}
		}
		else
		{
			double ar = gsl_rng_uniform(w->rgenref);
			PtclMap<int>::iterator itr = w->state.Rcurr->find(vect<int>(siteloc.x+1,0,0));
			if(itr==w->state.Rcurr->end() && siteloc.x != rp.L-1)
			{
				newhop++;
				int ptclnum = itr->second;
				w->state.Rcurr->erase(it);
				(*w->state.Rcurr)[vect<int>(siteloc.x+1,0,0)] = ptclnum;
			}
		}
	}
	else if(w->state.particleCount>0)
	{
//		fprintf(this->debugLog,"%f %f %d %d\n",lweight,rweight,lhop,rhop);
		double num = rp.trans.lmove*rp.trans.lmc;
		tweight = num + rp.trans.rmove*rp.trans.rmc;

		if(rd<num/tweight)
		{
			double ar = gsl_rng_uniform(w->rgenref);
			PtclMap<int>::iterator itl = w->state.Rcurr->find(vect<int>(siteloc.x+1,0,0));
			if(itl!=w->state.Rcurr->end())
			{
				newhop--;
				int ptclnum = itl->second;
				w->state.Rcurr->erase(itl);
				(*w->state.Rcurr)[siteloc] = ptclnum;
			}
		}
		else
		{
			double ar = gsl_rng_uniform(w->rgenref);
			PtclMap<int>::iterator itr = w->state.Rcurr->find(vect<int>(siteloc.x-1,0,0));
			if(itr!=w->state.Rcurr->end())
			{
				newhop++;
				int ptclnum = itr->second;
				w->state.Rcurr->erase(itr);
				(*w->state.Rcurr)[siteloc] = ptclnum;
			}
		}
	}

	w->state.dQ.x = newhop;
	w->state.ltime++;

//	fprintf(this->debugLog,"%d ================================\n",w->state.ltime);
//	w->state.weight.display(this->debugLog);
//	w->state.Q.x.addUpdate(newhop);
	w->state.dweight = tweight;
	w->state.weight.multUpdate(w->state.dweight);
//	w->state.weight.display(this->debugLog);
//	fprintf(this->debugLog,"=================================\n");
}

/////////////////////////////////////////////////
template void SNSMover<int,stringstream>::initialize(Walker<int,stringstream>*);
template void SNSMover<int,stringstream>::move(Walker<int,stringstream>* w);
