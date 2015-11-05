/*
 * base.cpp
 *
 *  Created on: Aug 19, 2014
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


#include "dmc.h"
#include "RunParameters.h"
#include "SNSMover.h"

using namespace std;
using namespace __gnu_cxx;


int setup(int rank, string baseSpecFile)
{
	RunParameters runParams;

	int status = runParams.loadFile(baseSpecFile);
	if(status == FILENOTFOUND)
		return FAIL;
	else if(status == DIMERROR)
		return DIMERROR;
	if(rank == 0)
		runParams.display();


	//Output Log
	stringstream sp;
	sp << rank;
	string logf = runParams.logFile + "_" + sp.str();
	FILE* log = fopen(logf.c_str(),"a");

	//COMM WORLD
	int totalProcs;
	MPI_Comm_size(MPI_COMM_WORLD,&totalProcs);

	//Time scaler
	double dt = 1.0/(runParams.L+2);

	// All ranks need a collector
	// For rank 0 the global collector and the mpi collector are used exclusively to receive
	// For slave ranks, global collector and mpi-collector are used for local assembly (i.e. per node)

	Weight* freeEnergy = new Weight(0.0,runParams.divisor,runParams.maxBranchWeight,runParams.minBranchWeight);
	WalkerState<int>* gwstate = new WalkerState<int>(runParams.dimension,*freeEnergy,log);
	vector<Observable<int>*> globalObs;
	//Parallel Observables
	vector<MPIObservable*> mpiGlobalObs;

	fprintf(log,"Setting up observables.\n");
	fflush(log);

	for(int i=0;i<runParams.observableCount;i++)
	{
		if(runParams.observableType[i].compare("Basic")==0)
		{
			BasicObs<int>* oo = new BasicObs<int>(rank,totalProcs,*gwstate,runParams.observableName[i],log,dt);
			globalObs.push_back(oo);
			mpiGlobalObs.push_back(oo);
		}
		else if(runParams.observableType[i].compare("Density")==0)
		{
			Density<int>* oo = new Density<int>(rank,totalProcs,*gwstate,runParams.observableName[i],log);
			globalObs.push_back(oo);
			mpiGlobalObs.push_back(oo);
		}
	}

	fprintf(log,"Observables setup.\n");
	fflush(log);

	//Other observables
	vector<core::Walker<int>*> walkerCollection;

	//Setup MPI Runner
	MPIBasicRunner<int> *brunner;

	if(rank == 0) //Reserve Master
	{
		fprintf(log,"Starting master run.\n");
		fflush(log);

		brunner = new MPIBasicRunner<int>(log,totalProcs,runParams.bins,runParams.nSteps,globalObs,mpiGlobalObs);
		brunner->masterRun();
	}
	else
	{

		gsl_rng_env_setup();
		gsl_rng* rgenref = gsl_rng_alloc(gsl_rng_mt19937);
		int procSeed;
		if(runParams.rinitseed == -1)
			procSeed = (int) time(NULL) + totalProcs*rank;
		else
			procSeed = (int) totalProcs*rank + runParams.rinitseed;
		gsl_rng_set(rgenref,procSeed);

		fprintf(log,"Randomized offset for seed: %d\n",procSeed);

		int localWalkerCount = runParams.walkerCount;

		fprintf(log,"Number of walkers for current process: %d\n",localWalkerCount);
		fflush(log);

		//Setup Mover
		SNSMover<int>* mov;
		if(runParams.moverName.compare("SNSMover")==0)
			mov = new SNSMover<int>(log,runParams);
		else
			return FAIL;

		//Total number of walkers
		long totalnwalkers = totalProcs*runParams.walkerCount;
		double initweight = 1.0/totalnwalkers;

		//setup walker threads
		for(int w=0;w<localWalkerCount;w++)
		{
			//////////////////////////////////////////////////////////////////////////////////////////
			//Prepare Walker State Objects
			//////////////////////////////////////////////////////////////////////////////////////////

			//State File
			Weight* wt = new Weight(initweight,runParams.divisor,runParams.maxBranchWeight,runParams.minBranchWeight);
			WalkerState<int>* wstate = new WalkerState<int>(runParams.dimension,*wt,log);

			//Observable Files
			vector<Observable<int>*>* localObs = new vector<Observable<int>*>;
			for(int i=0;i<runParams.observableCount;i++)
			{
				if(runParams.observableType[i].compare("Basic")==0)
					localObs->push_back(new BasicObs<int>(*wstate,runParams.observableName[i],log,dt));
				if(runParams.observableType[i].compare("Density")==0)
					localObs->push_back(new Density<int>(*wstate,runParams.observableName[i],log));
			}

			//Figure out seed
			//Total no. of samplers = totalProcs*maxWalkerCount
			int localSeed = runParams.maxWalkerCount*rank + w + procSeed;
			fprintf(log,"Seed: %d\n",localSeed);
			fflush(log);


			Walker<int>* lwalker = new Walker<int>(w,localSeed,*wstate,*localObs);
			walkerCollection.push_back(lwalker);
			//////////////////////////////////////////////////////////////////////////////////////////
		}

		brunner = new MPIBasicRunner<int>(log,totalProcs,mov,runParams.bins,runParams.nSteps,
				runParams.minBranchWeight,runParams.maxBranchWeight,&walkerCollection,globalObs,mpiGlobalObs);
		fprintf(log,"Starting spawn runs.\n");

		//Now run
		brunner->run();
	}

	//This is important. Wait for master to finish
	MPI_Barrier(MPI_COMM_WORLD);
	fclose(log);
}

