/*
 * runParams.h
 *
 *  Created on: Aug 25, 2014
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

#ifndef RUNPARAMS_H_
#define RUNPARAMS_H_

struct TransWeight
{
	float linsert,lremove;
	float rinsert,rremove;
	float lmove,rmove;

	float lmc,rmc;
};

struct RunParameters
{
	/////////////////////////////////////////
	int dimension;
	int L;
	float beta;
	/////////////////////////////////////////
	int walkerCount;
	int maxWalkerCount;
	int divisor;
	long double minBranchWeight;
	long double maxBranchWeight;
	/////////////////////////////////////////
	int rinitseed;
	/////////////////////////////////////////
	string logFile;
	/////////////////////////////////////////
	int observableCount;
	string* observableType;
	string* observableName;
	/////////////////////////////////////////
	int bins;
	int nSteps;
	/////////////////////////////////////////
	TransWeight trans;
	string moverName;
	/////////////////////////////////////////

	void display()
	{
		cout<<"===========================================================\n";
		cout<<"Run Parameters\n";
		cout<<"===========================================================\n";
		cout<<"Length: "<<L<<endl;
		cout<<"Dimesion: "<<dimension<<endl;
		cout<<"beta: "<<beta<<endl;
		cout<<"-----------------------------------------------------------\n";
		cout<<"Walker Count: "<<walkerCount<<endl;
		cout<<"Max. Walker Count: "<<maxWalkerCount<<endl;
		cout<<"Divisor: "<<divisor<<endl;
		cout<<"Min. Branch Weight: "<<minBranchWeight<<endl;
		cout<<"Max. Branch Weight: "<<maxBranchWeight<<endl;
		cout<<"-----------------------------------------------------------\n";
		cout<<"Random Seed: "<<rinitseed<<endl;
		cout<<"-----------------------------------------------------------\n";
		cout<<"Log file: "<<logFile<<endl;
		cout<<"-----------------------------------------------------------\n";
		cout<<"Observable Count: "<<observableCount<<endl;
		for(int i=0;i<observableCount;i++)
			cout<<observableType[i]<<" "<<observableName[i]<<"\n";
		cout<<"-----------------------------------------------------------\n";
		cout<<"Bins: "<<bins<<endl;
		cout<<"Steps: "<<nSteps<<endl;
		cout<<"===========================================================\n\n";

		cout<<"===========================================================\n";
		cout<<"Model Parameters\n";
		cout<<"===========================================================\n";
		cout<<"Mover name: "<<moverName<<endl;
		cout<<endl;
		cout<<"Left Insert: "<<trans.linsert<<endl;
		cout<<"Left Remove: "<<trans.lremove<<endl;
		cout<<"Right Insert: "<<trans.rinsert<<endl;
		cout<<"Right Remove: "<<trans.rremove<<endl;
		cout<<"Left Move: "<<trans.lmove<<endl;
		cout<<"Right Move: "<<trans.rmove<<endl;
		cout<<"===========================================================\n";
	}

	int loadFile(string baseFile)
	{
		ifstream bf(baseFile);
		///////////////////////////////////////////////////////////////

		if(!bf)
			return FILENOTFOUND;
		bf>>L;
		bf>>dimension;
#if (DIMENSION == 3)
		if(dimension != 3)
			return DIMERROR;
#elif (DIMENSION == 2)
		if(dimension != 2)
			return DIMERROR;
#elif (DIMENSION == 1)
		if(dimension != 1)
			return DIMERROR;
#endif
		bf>>beta;

		bf>>walkerCount;
		bf>>maxWalkerCount;
		bf>>divisor;
		bf>>minBranchWeight;
		bf>>maxBranchWeight;

		bf>>rinitseed;
		bf>>logFile;

		bf>>observableCount;
		observableType = new string[observableCount];
		observableName = new string[observableCount];
		for(int i=0;i<observableCount;i++)
			bf>>observableType[i]>>observableName[i];

		bf>>bins;
		bf>>nSteps; nSteps *= (L+2); //So that time is nSteps rather than sweeps

		///////////////////////////////////////////////////////////////
		bf>>moverName;
		bf>>trans.linsert;
		bf>>trans.lremove;
		bf>>trans.rinsert;
		bf>>trans.rremove;
		bf>>trans.lmove;
		bf>>trans.rmove;
		bf.close();

		trans.lmc = exp(beta);
		trans.rmc = exp(-beta);
		return SUCCESS;
	}

};

#endif /* RUNPARAMS_H_ */
