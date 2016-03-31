/*
 * BasicObs.h
 *
 *  Created on: Aug 20, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
 */

#ifndef BasicObs_H_
#define BasicObs_H_

#include "Observable.h"

namespace measures {

template <class T, class U>
class BasicObs: public measures::Observable<T,U>, public measures::MPIObservable {
public:
	double dt;
	vect<double> Q;

	//For global collection into processes
	unsigned int ltime;

	//for pend
	Weight Qx,Qy,Qz;
	Weight Qx2,Qy2,Qz2;
	Weight freeEnergy;

#ifndef NOBRANCH
	//Need collector also (for pavg)
	Weight Qax,Qay,Qaz;
	Weight Qax2,Qay2,Qaz2;
	Weight freeEnergya;
#endif

	//For global collection into master
	int Zcount;
	vect<double> Q2;
	vect<double> V;
	vect<double> V2;
	double fe, fe2;
#ifndef NOBRANCH
	vect<double> Qa;
	vect<double> Qa2;
	vect<double> Va;
	vect<double> Va2;
	double fea, fea2;
#endif

	BasicObs(core::WalkerState<T,U>& _state, string bsf, FILE* log) : Observable<T,U>(_state,bsf,log)
	{
		dt = 0.0;
		Zcount = 0;
		ltime = 0;
		Q.x = 0.0; Q2.x = 0.0;
		Q.y = 0.0; Q2.y = 0.0;
		Q.z = 0.0; Q2.z = 0.0;
		fe = fe2 = 0.0;
#ifndef NOBRANCH
		Qa.x = 0.0; Qa2.x = 0.0;
		Qa.y = 0.0; Qa2.y = 0.0;
		Qa.z = 0.0; Qa2.z = 0.0;
		fea = fea2 = 0.0;
#endif
	}

	BasicObs(core::WalkerState<T,U>& _state, string bsf, FILE* log, double _dt) : Observable<T,U>(_state,bsf,log)
	{
		dt = _dt;
		Zcount = 0;
		ltime = 0;
		Q.x = 0.0; Q2.x = 0.0;
		Q.y = 0.0; Q2.y = 0.0;
		Q.z = 0.0; Q2.z = 0.0;
		fe = fe2 = 0.0;
#ifndef NOBRANCH
		Qa.x = 0.0; Qa2.x = 0.0;
		Qa.y = 0.0; Qa2.y = 0.0;
		Qa.z = 0.0; Qa2.z = 0.0;
		fea = fea2 = 0.0;
#endif
	}

	BasicObs(int pId,int nprocs, core::WalkerState<T,U>& _state, string bsf, FILE* log, double _dt) : MPIObservable(pId,nprocs),Observable<T,U>(_state,bsf,log)
	{
		dt = _dt;
		Zcount = 0;
		ltime = 0;
		Q.x = 0.0; Q2.x = 0.0;
		Q.y = 0.0; Q2.y = 0.0;
		Q.z = 0.0; Q2.z = 0.0;
		fe = fe2 = 0.0;
#ifndef NOBRANCH
		Qa.x = 0.0; Qa2.x = 0.0;
		Qa.y = 0.0; Qa2.y = 0.0;
		Qa.z = 0.0; Qa2.z = 0.0;
		fea = fea2 = 0.0;
#endif
	}

	~BasicObs()
	{
	}

	void display();
	void measure();
	void writeViaIndex(int idx);
	void gather(void*);
	void branchGather(void*);
	void clear();
	Observable<T,U>* duplicate(core::WalkerState<T,U>&);
	void copy(void*);
	////////////////////////////////////////////////////////////////////////////////////////////
	//MPI Comm
	////////////////////////////////////////////////////////////////////////////////////////////

	int parallelSend(); //To be called by slaves
	int parallelReceive(); //To be called by master

	////////////////////////////////////////////////////////////////////////////////////////////
	//Serialization
	////////////////////////////////////////////////////////////////////////////////////////////
	virtual void serialize(Serializer<U>&);
	virtual void unserialize(Serializer<U>&);

};
} /* namespace measures */
#endif /* BasicObs_H_ */
