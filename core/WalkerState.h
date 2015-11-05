/*
 * WalkerState.h
 *
 *  Created on: Aug 20, 2015
 *      Author: ushnish
 *
 	Copyright (c) 2015 Ushnish Ray
	All rights reserved.
 */

#ifndef WALKERSTATE_H_
#define WALKERSTATE_H_

namespace core{

template <class T>
class WalkerState
{

public:
	int DIM;
	int particleCount;
	PtclMap<T>* Rcurr;
	vect<T> dQ;
	long ltime;
	Weight& weight;
	FILE* out;

	WalkerState(int _dim, Weight& _w, FILE* ot):weight(_w)
	{
		DIM = _dim;
		particleCount = 0;
		Rcurr = new PtclMap<T>;
		dQ.x = (T) 0.0;
		dQ.y = (T) 0.0;
		dQ.z = (T) 0.0;
		ltime = 0;
		out = ot;
	}

	WalkerState(int _dim,int _N,Weight& w):weight(w)
	{
		DIM = _dim;
		particleCount = _N;

		Rcurr = new PtclMap<T>;
		for(int i=0;i<_N;i++)
			(*Rcurr)[vect<T>(0.0,0.0,(T) i)] = i;

		dQ.x = (T) 0.0;
		dQ.y = (T) 0.0;
		dQ.z = (T) 0.0;
		ltime = 0;
	}

	WalkerState(int _dim,int _N, PtclMap<T>& Rcpy, vect<T> _dQ, long _time, Weight& w):weight(w)
	{
		DIM = _dim;
		particleCount = _N;

		Rcurr = new PtclMap<T>(Rcpy);

		dQ.x = _dQ.x;
		dQ.y = _dQ.y;
		dQ.z = _dQ.z;

		ltime = _time;
	}

	WalkerState(const WalkerState& ws):weight(*(new Weight(ws.weight)))
	{
		DIM = ws.DIM;
		particleCount = ws.particleCount;

		Rcurr = new PtclMap<T>(ws.Rcpy);

		dQ.x = ws.dQ.x;
		dQ.y = ws.dQ.y;
		dQ.z = ws.dQ.z;

		ltime = ws.ltime;
	}

	~WalkerState()
	{
		delete Rcurr;
		delete &weight;
	}

	WalkerState* duplicate()
	{
		WalkerState* a = new WalkerState(DIM,particleCount,*Rcurr,dQ,ltime,*(new Weight(weight)));
		return a;
	}

	void display()
	{
		fprintf(out,"Dimension: %d\n",DIM);
		fprintf(out,"Particle Count: %d\n",particleCount);
		fprintf(out,"dQ = (%d,%d,%d)\n",(int) dQ.x,(int) dQ.y,(int) dQ.z);
		fprintf(out,"time = %d\n",ltime);
		weight.display(out);
		fprintf(out,"State:\n");
		PtclMap<int>::iterator it;
		for(it=Rcurr->begin();it!=Rcurr->end();++it)
			fprintf(out,"(%d,%d,%d) = %d\n",it->first.x,it->first.y,it->first.z,it->second);
		fprintf(out,"---------------------------------------------------\n");
		fflush(out);
	}
};

}

#endif /* WALKERSTATE_H_ */
