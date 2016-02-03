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

template <class T>
void SNSMover<T>::initialize(Walker<T>* w)
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

	w->state.ltime = 0;
	w->state.weight.resetValue();

	//w->state.display();
}

template <class T>
void SNSMover<T>::move(Walker<T>* w)
{
	vect<int> siteloc(gsl_rng_uniform_int(w->rgenref,rp.L+2),0,0);

	double rd = gsl_rng_uniform(w->rgenref);
	double tweight = 1.0;
	int newhop = 0;

	bool ptcl = false;
	PtclMap<int>::iterator it = w->state.Rcurr->find(siteloc);
	if(it!=w->state.Rcurr->end())
		ptcl = true;

	if(siteloc.x == rp.L)
	{
		PtclMap<int>::iterator it = w->state.Rcurr->find(vect<int>(0,0,0));
		if(it!=w->state.Rcurr->end())
		{
			double temp = rp.trans.lremove*rp.trans.lmc;
			tweight = 1.0-rp.trans.lremove + temp;
			if(rd<temp/tweight)
			{
				w->state.Rcurr->erase(it);
				w->state.particleCount--;
				newhop--;
			}
		}
		else
		{
			double temp = rp.trans.linsert*rp.trans.rmc;
			tweight = 1.0-rp.trans.linsert + temp;
			if(rd<temp/tweight)
			{
				(*w->state.Rcurr)[vect<int>(0,0,0)] = w->state.particleCount++;
				newhop++;
			}
		}
	}
	else if(siteloc.x == rp.L+1)
	{
		PtclMap<int>::iterator it = w->state.Rcurr->find(vect<int>(rp.L-1,0,0));
		if(it!=w->state.Rcurr->end())
		{
			double temp = rp.trans.rremove*rp.trans.rmc;
			tweight = 1.0-rp.trans.rremove + temp;
			if(rd<temp/tweight)
			{
				w->state.Rcurr->erase(it);
				w->state.particleCount--;
				newhop++;
			}
		}
		else
		{
			double temp = rp.trans.rinsert*rp.trans.lmc;
			tweight = 1.0-rp.trans.rinsert + temp;
			if(rd<temp/tweight)
			{
				(*w->state.Rcurr)[vect<int>(rp.L-1,0,0)] = w->state.particleCount++;
				newhop--;
			}
		}
	}
	else if(ptcl)
	{
		double lweight = 0.0;
		int lhop = 0;
		PtclMap<int>::iterator itl = w->state.Rcurr->find(vect<int>(siteloc.x-1,0,0));
		if(itl!=w->state.Rcurr->end() || siteloc.x == 0) //Particle found at x-1 or at left edge
			lweight = rp.trans.lmove;
		else
		{
			lhop = -1;
			lweight = rp.trans.lmove*rp.trans.lmc;
		}

		double rweight = 0.0;
		int rhop = 0;
		PtclMap<int>::iterator itr = w->state.Rcurr->find(vect<int>(siteloc.x+1,0,0));
		if(itr!=w->state.Rcurr->end() || siteloc.x == rp.L-1) //Particle found at x+1 or at right edge
			rweight = rp.trans.rmove;
		else
		{
			rhop = 1;
			rweight = rp.trans.rmove*rp.trans.rmc;
		}
		tweight = lweight+rweight;

		if(rd<lweight/tweight && lhop!=0)
		{
			newhop--;
			int ptclnum = itl->second;
			w->state.Rcurr->erase(it);
			(*w->state.Rcurr)[vect<int>(siteloc.x-1,0,0)] = ptclnum;
		}
		else if(rd>=lweight/tweight && rhop!=0)
		{
			newhop++;
			int ptclnum = itr->second;
			w->state.Rcurr->erase(it);
			(*w->state.Rcurr)[vect<int>(siteloc.x+1,0,0)] = ptclnum;
		}
	}
	else if(w->state.particleCount>0)
	{
		double lweight = 0.0;
		int lhop = 0;
		PtclMap<int>::iterator itl = w->state.Rcurr->find(vect<int>(siteloc.x+1,0,0));
		if(itl!=w->state.Rcurr->end())
		{
			lhop = -1;
			lweight = rp.trans.lmove*rp.trans.lmc;
		}
		else
			lweight = rp.trans.lmove;

		double rweight = 0.0;
		int rhop = 0;
		PtclMap<int>::iterator itr = w->state.Rcurr->find(vect<int>(siteloc.x-1,0,0));
		if(itr!=w->state.Rcurr->end())
		{
			rweight = rp.trans.rmove*rp.trans.rmc;
			rhop = 1;
		}
		else
			rweight = rp.trans.rmove;

		tweight = lweight+rweight;

//		fprintf(this->debugLog,"%f %f %d %d\n",lweight,rweight,lhop,rhop);
		if(rd<lweight/tweight && lhop!=0)
		{
			newhop--;
			int ptclnum = itl->second;
			w->state.Rcurr->erase(itl);
			(*w->state.Rcurr)[siteloc] = ptclnum;
		}
		else if(rd>=lweight/tweight && rhop!=0)
		{
			newhop++;
			int ptclnum = itr->second;
			w->state.Rcurr->erase(itr);
			(*w->state.Rcurr)[siteloc] = ptclnum;
		}
	}

	w->state.dQ.x = newhop;
	w->state.ltime++;

//	fprintf(this->debugLog,"%d ================================\n",w->state.ltime);
//	w->state.weight.display(this->debugLog);
	w->state.weight.update(tweight);
//	w->state.weight.display(this->debugLog);
//	fprintf(this->debugLog,"=================================\n");
}

/////////////////////////////////////////////////
template void SNSMover<int>::initialize(Walker<int>*);
template void SNSMover<int>::move(Walker<int>* w);
