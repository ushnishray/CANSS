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

template <class T, typename U>
class WalkerState: public Serializable<U>
{

public:
	int DIM;
	int particleCount;
	PtclMap<T>* Rcurr;
	vect<T> dQ;
//	vect<Weight> Q;
	unsigned int ltime;
	double dweight;
	Weight weight;

	FILE* out;

	WalkerState(int _dim, Weight& _w, FILE* ot)
	{
		DIM = _dim;
		particleCount = 0;
		Rcurr = new PtclMap<T>;
		dQ.x = (T) 0.0;
		dQ.y = (T) 0.0;
		dQ.z = (T) 0.0;
		dweight = 0.0;
		ltime = 0;
		weight.copy(_w);
		out = ot;
	}

	WalkerState(int _dim,int _N,Weight& w, FILE* ot)
	{
		DIM = _dim;
		particleCount = _N;

		Rcurr = new PtclMap<T>;
		for(int i=0;i<_N;i++)
			(*Rcurr)[vect<T>(0.0,0.0,(T) i)] = i;

		dQ.x = (T) 0.0;
		dQ.y = (T) 0.0;
		dQ.z = (T) 0.0;
		dweight = 0.0;
		ltime = 0;
		weight.copy(w);

		out = ot;
	}

//	WalkerState(int _dim,int _N, PtclMap<T>& Rcpy, vect<T> _dQ,vect<Weight> _Q, long _time, Weight& w)
	WalkerState(int _dim,int _N, PtclMap<T>& Rcpy, vect<T> _dQ, long _time, Weight& w)
	{
		DIM = _dim;
		particleCount = _N;

		Rcurr = new PtclMap<T>(Rcpy);

		dQ.x = _dQ.x;
		dQ.y = _dQ.y;
		dQ.z = _dQ.z;
		dweight = 0.0;
//		Q.x = _Q.x;
//		Q.y = _Q.y;
//		Q.z = _Q.z;

		ltime = _time;
		weight.copy(w);
	}

	WalkerState(const WalkerState& ws)
	{
		DIM = ws.DIM;
		particleCount = ws.particleCount;

		Rcurr = new PtclMap<T>(ws.Rcpy);

		dQ.x = ws.dQ.x;
		dQ.y = ws.dQ.y;
		dQ.z = ws.dQ.z;
		dweight = 0.0;
//		Q.x = ws.Q.x;
//		Q.y = ws.Q.y;
//		Q.z = ws.Q.z;

		ltime = ws.ltime;
		weight.copy(ws.weight);
	}

	~WalkerState()
	{
		Rcurr->clear();
		delete Rcurr;
	}

	void copy(WalkerState<T,U>& w)
	{
		DIM = w.DIM;
		particleCount = w.particleCount;
		Rcurr->clear();
		for(typename PtclMap<T>::iterator it = w.Rcurr->begin(); it!=w.Rcurr->end();++it)
				(*Rcurr)[it->first] = it->second;
		dQ = w.dQ;
//		Q = w.Q;
		ltime = w.ltime;
		dweight = w.dweight;
		weight.copy(w.weight);
		out = w.out;
	}

	WalkerState* duplicate()
	{
//		WalkerState* a = new WalkerState(DIM,particleCount,*Rcurr,dQ,Q,ltime,weight);
		WalkerState* a = new WalkerState(DIM,particleCount,*Rcurr,dQ,ltime,weight);
		a->dweight = this->dweight;
		a->out = this->out;
		return a;
	}

	void display()
	{
		fprintf(out,"Dimension: %d\n",DIM);
		fprintf(out,"Particle Count: %d\n",particleCount);
		fprintf(out,"dQ = (%d,%d,%d)\n",(int) dQ.x,(int) dQ.y,(int) dQ.z);
//		fprintf(out,"Q = (%10.6e,%10.6e,%10.6e)\n",Q.x.logValue(),Q.y.logValue(),Q.z.logValue());
		fprintf(out,"time = %d\n",ltime);
		weight.display(out);
		fprintf(out,"State:\n");
		PtclMap<int>::iterator it;
		for(it=Rcurr->begin();it!=Rcurr->end();++it)
		{
			fprintf(out,"(%d,%d,%d) = %d\n",it->first.x,it->first.y,it->first.z,it->second);
		}
		fprintf(out,"---------------------------------------------------\n");
		fflush(out);
	}

	void reset()
	{
		dQ.x = (T) 0.0;
		dQ.y = (T) 0.0;
		dQ.z = (T) 0.0;
//		Q.x = 0.0;
//		Q.y = 0.0;
//		Q.z = 0.0;
		ltime = 0;
		dweight = 0.0;
		weight.resetValue();	
	}

	virtual void serialize(Serializer<U>& obj)
	{
		obj<<DIM<<dQ<<ltime<<particleCount<<Rcurr<<weight<<dweight;
	}

	virtual void unserialize(Serializer<U>& obj)
	{
		obj>>DIM>>dQ>>ltime>>particleCount>>Rcurr>>weight>>dweight;
	}
};

}

#endif /* WALKERSTATE_H_ */
