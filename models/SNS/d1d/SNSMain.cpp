/*
 * BHMain.cpp
 *
 *  Created on: Aug 26, 2014
 *      Author: ushnish
 *
 	Copyright (c) 2014 Ushnish Ray
	All rights reserved.
 */

#include "dmc.h"

#define MPI
using namespace std;

extern int setup(int rank, string baseSpecFile, int argc, char* argv[]);

int main(int argc,char* argv[])
{
#ifdef MPI
	if(argc<2)
	{
		cout<<"Format is: <filenamelist>"<<endl;
		return -1;
	}

	int rank;
	int prov;

	MPI_Status Stat;
	MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE, &prov);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0)
	{
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
	}

	string bfile(argv[1]);
	int status;
	status = setup(rank,bfile,argc,argv);

	MPI_Finalize();
#else


#endif
}
