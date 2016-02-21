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
	float* wProbs = new float[walkers.walkerCount];


	//Send to master
	Weight lw;
	for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
		lw.add(it->second->state.weight);

	//Now compute new population of walkers based on local weight
	int i = 0;
#if DEBUG >= 4
	fprintf(this->log,"++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif
	for(typename NumMap<Walker<T,U>>::iterator it = walkers.walkerCollection->begin();it!=walkers.walkerCollection->end();++it)
	{
		idx[i] = it->first;
		float probi = (it->second->state.weight/lw).value();
		wProbs[i] = probi;
		ni[i] = probi*walkers.walkerCount + gsl_rng_uniform(it->second->rgenref);
#if DEBUG >= 4
		fprintf(this->log,"Walker %d, new population %d [prob %10.6e]\n",idx[i],ni[i],probi);
#endif
		i++;
	}
#if DEBUG >= 4
	fprintf(this->log,"++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif

	memcpy(newPops,ni,sizeof(int)*walkers.walkerCount);

	//Now that we have total populations adjust them
	//Find extras and then redistribute among walkers that have 0 population with largest weight
	int extras = 0;
	double totalp = 0.0;
	vector<pair<int,float>> plist, tlist;

	int maxpopallowed = maxpercent*walkers.walkerCount;
	for(int j=0;j<walkers.walkerCount;j++)
	{
		if(ni[j]>maxpopallowed)
		{
			extras += ni[j] - maxpopallowed;
			newPops[j] = maxpopallowed;
		}
		else if(ni[j] == 0)
		{
			plist.push_back(pair<int,float>(j,wProbs[j]));
			totalp += wProbs[j];
		}

		tlist.push_back(pair<int,float>(j,wProbs[j]));
	}


	struct compobj
	{
		bool operator()(pair<int,float> a, pair<int,float> b)
		{
			return (a.second > b.second);
		}
	} co;
	sort(plist.begin(),plist.end(),co);
	for(int i = 0;i<extras;i++)
	{
		int si = plist[i].first;
		newPops[si] = 1;
	}

	///////////////////////////////////////////////////////////////////////////////
	// Need to do a check that total population is exactly right
	// This should be small - happens because of fractional occupation
	///////////////////////////////////////////////////////////////////////////////
	sort(tlist.begin(),tlist.end(),co);
//	fprintf(this->log,"################################################################\n");
	int newtotalwalkers = 0;
	for(int i = 0;i<tlist.size();i++)
	{
		int si = tlist[i].first;
//		float prb = tlist[i].second;

		newtotalwalkers += newPops[si];
//		fprintf(this->log,"Sorted Process: %d Walker: %d New pop: %d Prob: %10.6e\n",si,sj,newPops[si][sj],prb);
	}
//	fprintf(this->log,"################################################################\n");
//	fflush(this->log);

#if DEBUG >= 1
	fprintf(this->log,"################################################################\n");
	for(int j=0;j<walkers.walkerCount;j++)
			fprintf(this->log,"Walker: %d Population [Old] [New]: %d %d %10.6e\n",j,ni[j],newPops[j],wProbs[j]);

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

	}
#endif

#if DEBUG >= 1
	fprintf(this->log,"################################################################\n");
	for(int j=0;j<walkers.walkerCount;j++)
		fprintf(this->log,"Walker: %d Population [Old] [New]: %d %d %10.6e\n",j,ni[j],newPops[j],wProbs[j]);
		fprintf(this->log,"################################################################\n");
	fflush(this->log);
#endif

	memcpy(ni,newPops,sizeof(int)*walkers.walkerCount);

	delete[] newPops;
	delete[] wProbs;
	plist.clear();
	tlist.clear();

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

	delete[] idx;
	delete[] ni;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
template void MPIBasicRunner<int,stringstream>::branchLocal(float);
