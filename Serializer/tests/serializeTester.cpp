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

	DMCSNSSerializer<stringstream> bs;

	bs << "Hello World "<<13;
	cout<<"WS: "<<bs.printable(bs.str())<<endl;

	//Create character buffer
	bs.seekg(0,bs.end);
	int length = bs.tellg();
	cout<<"Size write: "<<length<<endl;
	bs.seekg(0,bs.beg);
	char* data = new char[length];
	bs.read(data,length);

	/*
	for(int i = 0;i<length;i++)
		cout<<i<<" "<<data[i]<<endl;
	*/

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

	/*
	for(int i = 0;i<length;i++)
		cout<<i<<" "<<buffer[i]<<endl;
	*/

	DMCSNSSerializer<stringstream> ns;
	ns.write(buffer,ll);
	cout<<"RS: "<<ns.printable(ns.str())<<endl;

	string a; int b;
	ns>>a>>b;
	cout<<a<<b<<endl;
}

void test2()
{
	//Testing composite object
}

int main(int argc, char* argv[])
{
	//	Weight* wt = new Weight(0.0);
	//	core::WalkerState<int> w(1,10,*wt);
	//	bs << w;


}


