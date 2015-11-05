/*
 * flagEnums.h
 *
 *  Created on: Apr 21, 2014
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

#ifndef BDATATYPES_H_
#define BDATATYPES_H_

#include "mainIncs.h"

enum bConditions {HW, PBC};
///////////////////////////////////////////////////////////////////////


struct eqnum{
	bool operator()(int num1, int num2)
	const
	{
		return num1 == num2;
	}
};

struct gcmpr{
	bool operator()(const int a, const int b)
	{
		return a<b;
	}
};

typedef unordered_map<int,double,hash<int>,eqnum> Row;
typedef map<long,double,gcmpr> sortedRow;
typedef map<int,int,gcmpr> IntMap;


typedef unordered_map<int,sortedRow,hash<int>,eqnum> orderedMap2d;
typedef unordered_map<int,Row,hash<int>,eqnum> hashmap2d;

typedef unordered_map<int,hashmap2d,hash<int>,eqnum> hashmap3d;

typedef hashmap2d map2d;
typedef hashmap3d map3d;

template <class T>
struct vect
{
	T x,y,z;

	vect() {x=y=z=0;}
	vect(T _x, T _y, T _z) : x(_x),y(_y),z(_z) {}

	bool operator<(vect<T> b)
	{
		return ((z == b.z) ? ((y == b.y) ? (x<b.x) : y<b.y) : z<b.z);
	}

	bool operator==(vect<T> b)
	{
		return ((z==b.z) && (y==b.y) && (x==b.x));
	}
};

template <class T>
struct vectComp
{
	bool operator()(vect<T> a, vect<T> b)
	{
		return a<b;
	}
};

template <class T>
struct vectorPairComp
{
	bool operator()(std::pair<vect<T>,vect<T>> p1, pair<vect<T>,vect<T>> p2)
	{
		bool (p1.first == p2.first) ? p1.second<p2.second : p1.first<p2.first;
	}
};

template <class T>
using vectToValue = map<vect<T>,double,vectComp<T>>;


template <class T>
using vectorMap2d = map<vect<T>,vectToValue<T>,vectComp<T>>;

template <class T>
using PtclMap = map<vect<T>,int,vectComp<T>>;

/*
struct ival
{
	int endLevel;
	float value;

	ival(): endLevel(0), value(0.0) {}
	ival(int l, float v): endLevel(l), value(v) { }
	ival(const ival& obj)
	{
		endLevel = obj.endLevel;
		value = obj.value;
	}
};

typedef unordered_map<int,ival,hash<int>,eqnum> Interval;
typedef unordered_map<int,Interval,hash<int>, eqnum>  SiteToInterval;
typedef unordered_map<int,SiteToInterval,hash<int>,eqnum> PairSiteToInterval;
typedef PairSiteToInterval GijI;
*/

#endif /* FLAGENUMS_H_ */
