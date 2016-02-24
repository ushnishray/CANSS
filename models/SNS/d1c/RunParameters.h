/*
 * runParams.h
 *
 *  Created on: Aug 25, 2014
 *      Author: ushnish
 *
 	Copyright (c) 2014 Ushnish Ray
	All rights reserved.
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
	int branchStep;
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
	int eSteps;
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
		cout<<"Walker Count (per process): "<<walkerCount<<endl;
		cout<<"Max. Walker Count (per process): "<<maxWalkerCount<<endl;
		cout<<"Branch interval (sweeps): "<<branchStep<<endl;
		cout<<"Log(Min. Branch Weight): "<<minBranchWeight<<endl;
		cout<<"Log(Max. Branch Weight): "<<maxBranchWeight<<endl;
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
		cout<<"E-Steps (max time) (sweeps): "<<eSteps/(L+2)<<" "<<nSteps<<endl;
		cout<<"D-Steps (max time) (sweeps): "<<nSteps/(L+2)<<" "<<nSteps<<endl;
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
		bf>>branchStep;
		bf>>minBranchWeight; minBranchWeight = log(minBranchWeight);
		bf>>maxBranchWeight; maxBranchWeight = log(maxBranchWeight);

		bf>>rinitseed;
		bf>>logFile;

		bf>>observableCount;
		observableType = new string[observableCount];
		observableName = new string[observableCount];
		for(int i=0;i<observableCount;i++)
			bf>>observableType[i]>>observableName[i];

		bf>>bins;
		bf>>eSteps;
		bf>>nSteps;
		eSteps *= (L+2); 
		nSteps *= (L+2); //So that nSteps provided is time rather than sweeps
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
