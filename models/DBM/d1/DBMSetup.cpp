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
#include "DBMMover.h"

using namespace std;
using namespace __gnu_cxx;


int setup(int rank, string baseSpecFile, int argc, char* argv[])
{
	RunParameters runParams;

	int status = runParams.loadFile(baseSpecFile);
	if(argc==3)
	{
		runParams.beta = atof(argv[2]);
	}

	if(status == FILENOTFOUND)
		return FAIL;
	else if(status == DIMERROR)
		return DIMERROR;
	if(rank == 0)
	{
		cout<<"Supplied Argument count: "<<argc<<endl;
		runParams.display();
	}

	//Output Log
	stringstream sp;
	sp << rank;
	string logf = runParams.logFile + "_" + sp.str();
	FILE* log = fopen(logf.c_str(),"a");

	//COMM WORLD
	int totalProcs;
	MPI_Comm_size(MPI_COMM_WORLD,&totalProcs);

	//Time scaler
	double dt = runParams.trans.dt;

	// All ranks need a collector
	// For rank 0 the global collector and the mpi collector are used exclusively to receive
	// For slave ranks, global collector and mpi-collector are used for local assembly (i.e. per node)

	Weight* freeEnergy = new Weight(0.0);
	WalkerState<float,stringstream>* gwstate = new WalkerState<float,stringstream>(runParams.dimension,*freeEnergy,log);
	vector<Observable<float,stringstream>*> globalObs;
	//Parallel Observables
	vector<MPIObservable*> mpiGlobalObs;

	fprintf(log,"Setting up observables.\n");
	fflush(log);

	for(int i=0;i<runParams.observableCount;i++)
	{
		if(runParams.observableType[i].compare("Basic")==0)
		{
			BasicObs<float,stringstream>* oo = new BasicObs<float,stringstream>(rank,totalProcs,*gwstate,runParams.observableName[i],log,dt);
			globalObs.push_back(oo);
			mpiGlobalObs.push_back(oo);
		}
		else if(runParams.observableType[i].compare("Density")==0)
		{
			Density<float,stringstream>* oo = new Density<float,stringstream>(rank,totalProcs,*gwstate,runParams.observableName[i],log);
			globalObs.push_back(oo);
			mpiGlobalObs.push_back(oo);
		}
		else if(runParams.observableType[i].compare("Qhistogram")==0)
		{
			Qhistogram<float,stringstream>* oo = new Qhistogram<float,stringstream>(rank,totalProcs,*gwstate,runParams.observableName[i],log,dt);
			globalObs.push_back(oo);
			mpiGlobalObs.push_back(oo);
		}
		else if(runParams.observableType[i].compare("Whistogram")==0)
		{
			Whistogram<float,stringstream>* oo = new Whistogram<float,stringstream>(rank,totalProcs,*gwstate,runParams.observableName[i],log,dt);
			globalObs.push_back(oo);
			mpiGlobalObs.push_back(oo);
		}
		else if(runParams.observableType[i].compare("AutoCorr")==0)
		{
			AutoCorr<float,stringstream>* oo = new AutoCorr<float,stringstream>(rank,totalProcs,*gwstate,runParams.observableName[i],log);
			globalObs.push_back(oo);
			mpiGlobalObs.push_back(oo);
		}
	}

	fprintf(log,"Observables setup.\n");
	fflush(log);

	//Other observables
	NumMap<Walker<float,stringstream>> walkerCollection;

	//Setup MPI Runner
	MPIBasicRunner<float,stringstream> *brunner;

	if(rank == 0) //Reserve Master
	{
		fprintf(log,"Starting master run.\n");
		fflush(log);

		gsl_rng_env_setup();
		int procSeed;
		if(runParams.rinitseed == -1)
			runParams.rinitseed = (int) time(NULL);
		gsl_rng* rgenref = gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(rgenref,runParams.rinitseed + rank);
		fprintf(log,"Seed: %d\n",runParams.rinitseed + rank);

		brunner = new MPIBasicRunner<float,stringstream>(log,totalProcs,
				runParams.eSteps,runParams.bins,runParams.nSteps,runParams.branchStep,
				runParams.walkerCount,globalObs,mpiGlobalObs,rgenref);
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
		DBMMover<float,stringstream>* mov;
		if(runParams.moverName.compare("DBMMover")==0)
			mov = new DBMMover<float,stringstream>(log,runParams);
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
			Weight* wt = new Weight(initweight);
			WalkerState<float,stringstream>* wstate = new WalkerState<float,stringstream>(runParams.dimension,*wt,log);

			//Observable Files
			vector<Observable<float,stringstream>*>* localObs = new vector<Observable<float,stringstream>*>;
			for(int i=0;i<runParams.observableCount;i++)
			{
				if(runParams.observableType[i].compare("Basic")==0)
					localObs->push_back(new BasicObs<float,stringstream>(*wstate,runParams.observableName[i],log,dt));
				else if(runParams.observableType[i].compare("Density")==0)
					localObs->push_back(new Density<float,stringstream>(*wstate,runParams.observableName[i],log));
				else if(runParams.observableType[i].compare("Qhistogram")==0)
					localObs->push_back(new Qhistogram<float,stringstream>(*wstate,runParams.observableName[i],log,dt));
				else if(runParams.observableType[i].compare("Whistogram")==0)
					localObs->push_back(new Whistogram<float,stringstream>(*wstate,runParams.observableName[i],log,dt));
				else if(runParams.observableType[i].compare("AutoCorr")==0)
					localObs->push_back(new AutoCorr<float,stringstream>(*wstate,runParams.observableName[i],log));
			}

			Walker<float,stringstream>* lwalker = new Walker<float,stringstream>(rgenref,*wstate,*localObs);
			walkerCollection[w] = lwalker;
			//////////////////////////////////////////////////////////////////////////////////////////
		}

		brunner = new MPIBasicRunner<float,stringstream>(log,totalProcs,mov,runParams.eSteps,runParams.bins,
				runParams.nSteps,runParams.branchStep,
				runParams.maxWalkerCount, runParams.maxBranchWeight,runParams.minBranchWeight,
				&walkerCollection,globalObs,mpiGlobalObs);
		fprintf(log,"Starting spawn runs.\n");
		fflush(log);

		//Now run
		brunner->run();
	}

	//This is important. Wait for master to finish
	MPI_Barrier(MPI_COMM_WORLD);
	fclose(log);
}

