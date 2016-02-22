/*
 * branchLocal.cpp
 *
 *  Created on: Feb 20, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */

#include "dmc.h"

template <class T, class U>
void MPIBasicRunner<T,U>::branchLocal(float maxpercent)
{
	int* idx = new int[walkers.walkerCount];
	int* ni = new int[walkers.walkerCount];
	int* newPops = new int[walkers.walkerCount];

	//Send to master
	Weight lw;
	for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
		lw.add(it->second->state.weight);

	//Now compute new population of walkers based on local weight
	int i = 0, newtotalwalkers = 0;
#if DEBUG >= 4
	fprintf(this->log,"++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif
	for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
	{
		idx[i] = it->first;
		float probi = (it->second->state.weight/lw).value();
		ni[i] = probi*walkers.walkerCount + gsl_rng_uniform(it->second->rgenref);
#if DEBUG >= 4
		fprintf(this->log,"Walker %d, new population %d [prob %10.6e]\n",idx[i],ni[i],probi);
#endif
		newtotalwalkers += ni[i];
		i++;
	}
#if DEBUG >= 4
	fprintf(this->log,"++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif

	memcpy(newPops,ni,sizeof(int)*walkers.walkerCount);

#if DEBUG >= 1
	fprintf(this->log,"################################################################\n");
	for(int j=0;j<walkers.walkerCount;j++)
			fprintf(this->log,"Walker: %d Population [Old] [New]: %d %d\n",j,ni[j],newPops[j]);

	fprintf(this->log,"################################################################\n");
	fprintf(this->log,"New population: %d Expected: %d\n",newtotalwalkers,walkers.walkerCount);
	fflush(this->log);
#endif

#if 0

	if(newtotalwalkers>walkers.walkerCount)
	{
		//Find smallest probability walkers and reduce them to keep total
		//population constant

		int tsize = tlist.size()-1;
		for(int i = 0,j = 0;i<newtotalwalkers-walkers.walkerCount;i++)
		{
			int si = tlist[tsize-j].first;

			while(newPops[si]<1)
			{
				j++;
				//cycle if needed;
				j = (j>tsize) ? 0: j;
				si = tlist[tsize-j].first;
			}

			newPops[si]--;
		}
	}
	else if(newtotalwalkers<walkers.walkerCount)
	{
		//Find largest probability walkers and reduce them to keep total
		//population constant
		for(int i = 0;i<walkers.walkerCount-newtotalwalkers;i++)
		{
			int si = tlist[i].first;
			newPops[si]++;
		}
	}
#else

	//Randomly copy or destroy if not right number of walkers
	int incr = (newtotalwalkers>walkers.walkerCount) ? -1: 1;
	int diff = abs(newtotalwalkers-walkers.walkerCount);
	for(int i = 0;i<diff;i++)
	{
		int ridx = 0;
		do{
			ridx = gsl_rng_uniform_int(walkers[0]->rgenref,walkers.walkerCount);
		} while(incr<0 && newPops[ridx]<=0);
		newPops[ridx] += incr;
	}
#endif

#if DEBUG >= 1
	fprintf(this->log,"################################################################\n");
	for(int j=0;j<walkers.walkerCount;j++)
		fprintf(this->log,"Walker: %d Population [Old] [New]: %d %d\n",j,ni[j],newPops[j]);
		fprintf(this->log,"################################################################\n");
	fflush(this->log);
#endif

	memcpy(ni,newPops,sizeof(int)*walkers.walkerCount);
	delete[] newPops;

	int* zvals = new int[walkers.walkerCount];

	//Compute local redistribution
	int ccount = 0, zc=0;
	i = 0;
	for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
	{
		//Keep track of zero values
		if(ni[i]==0)
			zvals[zc++] = idx[i];
		ccount += ni[i++];
	}

	//Copy over
	int p = 0;
	for(int i=0;i<walkers.walkerCount && p<zc;i++)
	{
		while(ni[i]>1)
		{
			(*walkers.walkerCollection)[zvals[p++]]->copy(*walkers[idx[i]]);
			ni[i]--;
		}
	}

	//Ensure that all 0 occupancies are copied into
	assert(p == zc);

	delete[] idx;
	delete[] ni;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
template void MPIBasicRunner<int,stringstream>::branchLocal(float);
