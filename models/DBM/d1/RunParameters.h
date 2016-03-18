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
	float ga;
	float dt;
	float V0;

	void display()
	{
		cout<<"===========================================================\n";
		cout<<"Model Parameters\n";
		cout<<"===========================================================\n";
		cout<<"Gamma : "<<ga<<endl;
		cout<<"dt: "<<dt<<endl;
		cout<<"V0: "<<V0<<endl;
		cout<<"===========================================================\n";
	}

	void load(ifstream &bf)
	{
		bf>>ga;
		bf>>dt;
		bf>>V0;
	}
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
		cout<<"E-Steps (max time) (sweeps): "<<eSteps*trans.dt<<" "<<nSteps<<endl;
		cout<<"D-Steps (max time) (sweeps): "<<nSteps*trans.dt<<" "<<nSteps<<endl;
		cout<<"===========================================================\n\n";
		cout<<"Mover name: "<<moverName<<endl;
		cout<<"===========================================================\n\n";
		cout<<endl;

		trans.display();
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
		///////////////////////////////////////////////////////////////
		bf>>moverName;
		trans.load(bf);

		//We are specifying final time
		eSteps /= trans.dt;
#ifdef NOBRANCH
		branchStep /= trans.dt;
#else
		nSteps /= trans.dt;
#endif

		return SUCCESS;
	}

};

#endif /* RUNPARAMS_H_ */
