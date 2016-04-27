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
	float J;
	float mu;
	float dt;

	float lmc,rmc;

	void display()
	{
		cout<<"===========================================================\n";
		cout<<"Model Parameters\n";
		cout<<"===========================================================\n";
		cout<<"J : "<<J<<endl;
		cout<<"mu: "<<mu<<endl;
		cout<<"===========================================================\n";
	}

	void load(ifstream &bf)
	{
		bf>>J;
		bf>>mu;
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
	int branchStep;
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
		cout<<"E-Steps (sweeps): "<<eSteps<<endl;
		cout<<"D-Steps (max time) (sweeps): "<<nSteps*trans.dt<<" "<<nSteps<<endl;
		cout<<"Integration Time (time) (sweeps): "<<branchStep*trans.dt<<" "<<branchStep<<endl;
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

		float temp;
		bf>>temp;
		///////////////////////////////////////////////////////////////
		bf>>moverName;
		trans.load(bf);
		trans.dt = 1.0/L;

		eSteps = temp/trans.dt*eSteps;
		branchStep = temp/trans.dt;

		trans.lmc = exp(beta);
		trans.rmc = exp(-beta);

		return SUCCESS;
	}

};

#endif /* RUNPARAMS_H_ */
