/*
 * serializeTester.cpp
 *
 *  Created on: Feb 11, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */

#include <iostream>
#include <sstream>
#include <iomanip>

#include "dmc.h"
#include "Serializer.h"
#include "DMCSNSSerializer.h"


using namespace std;

void test1()
{
	//Simple test of input and output of primitives

	Serializer<stringstream> bs;
	unsigned int b = 13;
	bs <<"Hello World"<<b;
	cout<<"WS: "<<bs.printable(bs.str())<<endl;

	//Create character buffer
	bs.seekg(0,bs.end);
	int length = bs.tellg();
	cout<<"Size write: "<<length<<endl;
	bs.seekg(0,bs.beg);
	char* data = new char[length];
	bs.read(data,length);

	for(int i = 0;i<length;i++)
		cout<<i<<" "<<data[i]<<endl;

	//Write to HDD
	ofstream of("Test.txt",ios::binary);
	of.write(data,length);
	of.close();
	delete[] data;

	//Read back to another stream
	cout<<"\nReading back:\n";

	ifstream f("Test.txt",ios::binary);
	f.seekg(0,f.end);
	int ll = f.tellg();
	cout<<"Size read: "<<ll<<endl;

	char* buffer = new char[ll];
	f.seekg(0,f.beg);
	f.read(buffer,ll);
	f.close();


	for(int i = 0;i<length;i++)
		cout<<i<<" "<<buffer[i]<<endl;


	Serializer<stringstream> ns;
	ns.write(buffer,ll);
	cout<<"RS: "<<ns.printable(ns.str())<<endl;

	string a; unsigned int bb;
	ns>>a>>bb;
	cout<<a<<bb<<endl;
}

#if 1
void test2()
{
	gsl_rng_env_setup();
	gsl_rng* rgenref = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rgenref,1);

	FILE* logf = fopen("out.txt","a");

#if 1
	//Write part
	double dt = 0.001;
	Serializer<stringstream> s1;
	//Testing composite object
	Weight* wt = new Weight(10.0);
	WalkerState<int>& wstate = *(new WalkerState<int>(1,*wt,logf));

	/*
	Serializer<stringstream> ns;
	ns<<wstate;
	ns.seekg(0,ns.end);
	int size = ns.tellg();
	cout<<"Size: "<<size<<endl;
	char* data = new char[size];
	ns.seekg(0,ns.beg);
	ns.read(data,size);
	ofstream f("data.out",ios::binary);
	f.write(data,size);
	f.close();
	cout<<"Written\n";
	delete[] data;
	wstate.display();
	delete &wstate;

	ifstream fr("data.out",ios::binary);
	fr.seekg(0,fr.end);
	int ll = fr.tellg();
	cout<<"Size to read: "<<ll<<endl;
	char* buffer = new char[ll];
	fr.seekg(0,fr.beg);
	fr.read(buffer,ll);
	fr.close();
	cout<<"Read data\n";

	Serializer<stringstream> rs;
	rs.write(buffer,ll);
	WalkerState<int>& wstate2 = *(new WalkerState<int>(1,*(new Weight(1.0)),logf));
	cout<<"Attempting to un-stream\n";
	rs>>wstate2;
	wstate2.display();
	delete &wstate2;
	*/

	//Observable Files
	vector<Observable<int>*>* localObs = new vector<Observable<int>*>;
	localObs->push_back(new BasicObs<int>(wstate,"Basic.txt",logf,dt));
	localObs->push_back(new Density<int>(wstate,"Density.txt",logf));
	localObs->push_back(new Qhistogram<int>(wstate,"Qhist.txt",logf,dt));
	localObs->push_back(new Whistogram<int>(wstate,"Whist.txt",logf,dt));
	Walker<int>* lwalker = new Walker<int>(rgenref,wstate,*localObs);

	lwalker->display();

	DMCSNSSerializer<stringstream> ns;
	ns<<*lwalker;
	ns.seekg(0,ns.end);
	int size = ns.tellg();
	ns.seekg(0,ns.beg);//Restore head
	cout<<"Write Size: "<<size<<endl;
	char* data = new char[size];
	ns.read(data,size);
	ofstream f("data.out",ios::binary);
	f.write(data,size);
	cout<<"Written\n";
	delete[] data;

#else
	//Read Part

	Weight* wt = new Weight(0.0);
	WalkerState<int>& wstate = *(new WalkerState<int>(1,*wt,logf));
	vector<Observable<int>*>* localObs = new vector<Observable<int>*>;
	Walker<int>* lwalker = new Walker<int>(rgenref,wstate,*localObs);

	//Read from data file
	ifstream f("data.out",ios::binary);
	f.seekg(0,f.end);
	int size = f.tellg();
	f.seekg(0,f.beg);
	char* data = new char[size];
	f.read(data,size);
	cout<<"Read Size: "<<size<<endl;

	Serializer<stringstream> ns;
	ns.write(data,size);
	delete[] data;

	ns>>*lwalker;
	lwalker->display();

#endif
	fclose(logf);
}
#endif

int main(int argc, char* argv[])
{
	//	Weight* wt = new Weight(0.0);
	//	core::WalkerState<int> w(1,10,*wt);
	//	bs << w;

//	test1();
	test2();

}


