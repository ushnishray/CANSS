/*
 * RSOSObs.h
 *
 *  Created on: Aug 20, 2014
 *      Author: ushnish
 *
 Copyright (c) 2014 Ushnish Ray
 All rights reserved.
 */

#ifndef RSOSObs_H_
#define RSOSObs_H_

#include "Observable.h"

namespace measures {

template <class T, class U>
class RSOSObs: public measures::Observable<T,U>, public measures::MPIObservable {
public:
	double dt;
	vect<double> Q;
	int N;

	//For global collection into processes
	unsigned int ltime;

	//for pend
	Weight Qx;
	Weight Qx2;
	Weight freeEnergy;
	Weight eH,eN;

	//For global collection into master
	int Zcount;
	vect<double> Q2;
	vect<double> V;
	vect<double> V2;
	double fe, fe2;
	double H,H2;
	double gN,N2;

#ifndef NOBRANCH
	//For averaging
	vector<vect<double>> cavgQ;
	vector<double> avgH;
	vector<int> avgN;

	//Global gathering
	vector<vect<double>> cavgQ2;
	vector<double> avgH2;
	vector<int> avgN2;
#endif

	RSOSObs(core::WalkerState<T,U>& _state, string bsf, FILE* log) : Observable<T,U>(_state,bsf,log)
	{

		dt = 0.0;
		Zcount = 0;
		ltime = 0;
		Q.x = 0.0; Q2.x = 0.0;
		fe = fe2 = 0.0;
#ifndef NOBRANCH
#endif
	}

	RSOSObs(core::WalkerState<T,U>& _state, string bsf, FILE* log, double _dt) : Observable<T,U>(_state,bsf,log)
	{
		dt = _dt;
		Zcount = 0;
		ltime = 0;
		Q.x = 0.0; Q2.x = 0.0;
		fe = fe2 = 0.0;
#ifndef NOBRANCH
#endif
	}

	RSOSObs(int pId,int nprocs,int tw, core::WalkerState<T,U>& _state, string bsf, FILE* log, double _dt) : MPIObservable(pId,nprocs,tw),Observable<T,U>(_state,bsf,log)
	{
		dt = _dt;
		Zcount = 0;
		ltime = 0;
		Q.x = 0.0; Q2.x = 0.0;
		fe = fe2 = 0.0;
#ifndef NOBRANCH
#endif
	}

	~RSOSObs()
	{
#ifndef NOBRANCH
		cavgQ.clear();
		cavgQ2.clear();
		avgH.clear();
		avgN.clear();
#endif
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
#endif /* RSOSObs_H_ */
