/*
 * BHMain.cpp
 *
 *  Created on: Aug 26, 2014
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

#define MPI
using namespace std;

extern int setup(int rank, string baseSpecFile);

int main(int argc,char* argv[])
{
#ifdef MPI
	if(argc<2)
	{
		cout<<"Format is: <Run mode> <filenamelist> [<prefix file>]"<<endl;
		return -1;
	}

	int rank;
	int prov;

	MPI_Status Stat;
	MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE, &prov);

	if(prov == MPI_THREAD_SINGLE)
		cout<<"Thread Support: SINGLE"<<endl;
	else if(prov == MPI_THREAD_FUNNELED)
		cout<<"Thread Support: FUNNEL"<<endl;
	else if(prov == MPI_THREAD_SERIALIZED)
		cout<<"Thread Support: SERIALIZED"<<endl;
	else if(prov == MPI_THREAD_MULTIPLE)
		cout<<"Thread Support: MULTIPLE"<<endl;
	else
		cout<<"Thread Support: ERROR"<<endl;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	string bfile(argv[1]);
	int status;
	status = setup(rank,bfile);

	MPI_Finalize();
#else


#endif
}
