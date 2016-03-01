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
#include "DBMMover.h"

template <class T, class U>
void DBMMover<T,U>::initialize(Walker<T,U>* w)
{
	w->state.Rcurr->clear();

	w->state.particleCount = 1;
	vect<float> vv(gsl_rng_uniform(w->rgenref)*this->rp.L,0.0,0.0);
	(*w->state.Rcurr)[vv] = 0;
	w->state.ltime = 0;
	w->state.reset();

//	fprintf(this->debugLog,"Initialized.\n");
//	fflush(this->debugLog);
}

template <class T, class U>
void DBMMover<T,U>::move(Walker<T,U>* w)
{
	PtclMap<float>::iterator it = w->state.Rcurr->begin();

	vect<float> vv = it->first;
	float x0 = vv.x;
	float noise = sqrt(2.0*rp.trans.dt)*gsl_ran_ugaussian(w->rgenref);
	vv.x = rp.trans.dt*(rp.trans.ga + rp.trans.V0*amp*sin(amp*vv.x)) + noise;

	w->state.dQ.x = vv.x*rp.trans.ga;
	w->state.dweight = exp(-rp.beta*w->state.dQ.x);
	w->state.weight.multUpdate(w->state.dweight);
	w->state.ltime++;

	//Update coordinates by first removing particle
	w->state.Rcurr->erase(it);
	vv.x += x0;

//	fprintf(this->debugLog,"%d, %10.6e -> %10.6e, W: %10.6e\n",w->state.ltime,vv.x,x0,w->state.dweight);
//	fflush(this->debugLog);

	while(vv.x<0.0)
		vv.x += rp.L;

	while(vv.x>=rp.L)
		vv.x -= rp.L;

	//Insert particle back
	(*w->state.Rcurr)[vv] = 0;
}

/////////////////////////////////////////////////
template void DBMMover<float,stringstream>::initialize(Walker<float,stringstream>*);
template void DBMMover<float,stringstream>::move(Walker<float,stringstream>* w);
