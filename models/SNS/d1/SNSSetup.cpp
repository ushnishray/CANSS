/*
 * base.cpp
 *
 *  Created on: Aug 19, 2014
 *      Author: ushnish
 *
 	Copyright (c) 2014 Ushnish Ray
	All rights reserved.
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
		else if(runParams.observableType[i].compare("Qhistogram")==0)
		{
			Qhistogram<int>* oo = new Qhistogram<int>(rank,totalProcs,*gwstate,runParams.observableName[i],log,dt);
			globalObs.push_back(oo);
			mpiGlobalObs.push_back(oo);
		}
		else if(runParams.observableType[i].compare("Whistogram")==0)
		{
			Whistogram<int>* oo = new Whistogram<int>(rank,totalProcs,*gwstate,runParams.observableName[i],log,dt);
			globalObs.push_back(oo);
			mpiGlobalObs.push_back(oo);
		}
	}

	fprintf(log,"Observables setup.\n");
	fflush(log);

	//Other observables
	NumMap<Walker<int>> walkerCollection;

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
		//Setup Random Seed for process
		gsl_rng_env_setup();
		int procSeed;
		if(runParams.rinitseed == -1)
			runParams.rinitseed = (int) time(NULL);

		gsl_rng* rgenref = gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(rgenref,runParams.rinitseed + rank);
		fprintf(log,"Seed: %d\n",runParams.rinitseed + rank);
		fprintf(log,"Number of walkers for current process: %d\n",runParams.walkerCount);
		fflush(log);

		//Setup Mover
		SNSMover<int>* mov;
		if(runParams.moverName.compare("SNSMover")==0)
			mov = new SNSMover<int>(log,runParams);
		else
			return FAIL;

		//Total number of walkers
		long totalnwalkers = (totalProcs-1)*runParams.walkerCount;
		double initweight = 1.0/totalnwalkers;

		//setup walker threads
		for(int w=0;w<runParams.walkerCount;w++)
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
				else if(runParams.observableType[i].compare("Density")==0)
					localObs->push_back(new Density<int>(*wstate,runParams.observableName[i],log));
				else if(runParams.observableType[i].compare("Qhistogram")==0)
					localObs->push_back(new Qhistogram<int>(*wstate,runParams.observableName[i],log,dt));
				else if(runParams.observableType[i].compare("Whistogram")==0)
					localObs->push_back(new Whistogram<int>(*wstate,runParams.observableName[i],log,dt));
			}

			Walker<int>* lwalker = new Walker<int>(rgenref,*wstate,*localObs);
			walkerCollection[w] = lwalker;
			//////////////////////////////////////////////////////////////////////////////////////////
		}

		brunner = new MPIBasicRunner<int>(log,totalProcs,mov,runParams.bins,runParams.nSteps,runParams.maxWalkerCount,
				runParams.maxBranchWeight,runParams.minBranchWeight,&walkerCollection,globalObs,mpiGlobalObs);
		fprintf(log,"Starting spawn runs.\n");

		//Now run
		brunner->run();
	}

	//This is important. Wait for master to finish
	MPI_Barrier(MPI_COMM_WORLD);
	fclose(log);
}

