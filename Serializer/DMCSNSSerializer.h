/*
 * DMCSNSSerializer.h
 *
 *  Created on: Feb 11, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */

#ifndef DMCSNSSERIALIZER_H_
#define DMCSNSSERIALIZER_H_

template <typename XStream>
class DMCSNSSerializer:public Serializer<XStream>
{
public:
	template<typename A>
	using B = Serializer<A>;

	//////////////////////////////////////////////////////////////////
	//Vect
	//////////////////////////////////////////////////////////////////

	template <class T>
	friend B<XStream>& operator<<(B<XStream>& bs, vect<T>& obj)
	{
		bs<<obj.x<<obj.y<<obj.z;
		return bs;
	}

	template <class T>
	friend B<XStream>& operator>>(B<XStream>& bs, vect<T>& obj)
	{
		bs>>obj.x>>obj.y>>obj.z;
		return bs;
	}

	//////////////////////////////////////////////////////////////////
	//Particle Map
	//////////////////////////////////////////////////////////////////

	template <class T>
	friend B<XStream>& operator<<(B<XStream>& bs, PtclMap<T>& obj)
	{
		bs<<obj.size();
		typename PtclMap<T>::iterator it;
		for(it=obj->begin();it!=obj->end();++it)
			bs<<(*it->first)<<it->second;
		return bs;
	}

	template <class T>
	friend B<XStream>& operator>>(B<XStream>& bs, PtclMap<T>& obj)
	{
		int tsize;
		bs>>tsize;
		for(int i = 0;i<tsize;i++)
		{
			vect<T> v;
			int val;
			bs>>v>>val;
			obj[v] = val;
		}
		return bs;
	}

	//////////////////////////////////////////////////////////////////
	//Vector
	//////////////////////////////////////////////////////////////////

	template <class T>
	friend B<XStream>& operator<<(B<XStream>& bs, vector<T>& obj)
	{
		bs>>obj.size();
		typename vect<T>::iterator it;
		for(it=obj->begin();it!=obj->end();++it)
			bs<<(*it);

		return bs;
	}

	template <class T>
	friend B<XStream>& operator>>(B<XStream>& bs, vector<T>& obj)
	{
		int ts;
		bs>>ts;
		for(int i=0;i<ts;i++)
		{
			T val;
			bs>>val;
			obj[i] = val;
		}
		return bs;
	}

	//////////////////////////////////////////////////////////////////
	//Vector-To-Value
	//////////////////////////////////////////////////////////////////

	template <class T>
	friend B<XStream>& operator<<(B<XStream>& bs, vectToValue<T>& obj)
	{
		bs>>obj.size();
		typename PtclMap<T>::iterator it;
		for(it=obj->begin();it!=obj->end();++it)
			bs<<(*it->first)<<it->second;
		return bs;
	}

	template <class T>
	friend B<XStream>& operator>>(B<XStream>& bs, vectToValue<T>& obj)
	{
		int tsize;
		bs>>tsize;
		for(int i = 0;i<tsize;i++)
		{
			vect<T> v;
			double val;
			bs>>v>>val;
			obj[v] = val;
		}
		return bs;
	}

	//////////////////////////////////////////////////////////////////
	//Weight
	//////////////////////////////////////////////////////////////////

	friend B<XStream>& operator<<(B<XStream>& bs, Weight& w)
	{
		bs<<w.divisor<<w.exponent<<w.initval<<w.maxvalue<<w.minvalue<<w.val;
		return bs;
	}

	friend B<XStream>& operator>>(B<XStream>& bs, Weight& w)
	{
		bs>>w.divisor>>w.exponent>>w.initval>>w.maxvalue>>w.minvalue>>w.val;
		return bs;
	}

	//////////////////////////////////////////////////////////////////
	//WalkerState
	//////////////////////////////////////////////////////////////////

	template <class T>
	friend B<XStream>& operator<<(B<XStream>& bs, WalkerState<T>& obj)
	{
		bs<<obj.DIM<<obj.Rcurr<<obj.Rcurr<<obj.dQ<<obj.ltime
				<<obj.particleCount<<obj.weight;
		return bs;
	}

	template <class T>
	friend B<XStream>& operator>>(B<XStream>& bs, WalkerState<T>& obj)
	{
		bs>>obj.DIM>>obj.Rcurr>>obj.Rcurr>>obj.dQ>>obj.ltime
				>>obj.particleCount>>obj.weight;
		return bs;
	}

	//////////////////////////////////////////////////////////////////
	//Observables: BasicObs
	//////////////////////////////////////////////////////////////////

	template <class T>
	friend B<XStream>& operator<<(B<XStream>& bs, BasicObs<T>& obj)
	{
		bs<<obj.dt<<obj.Q<<obj.Zcount<<obj.ltime<<obj.FreeEnergy<<obj.Qx<<obj.Qy<<obj.Qz<<obj.Q2x<<obj.Q2y<<obj.Q2z;
		return bs;
	}

	template <class T>
	friend B<XStream>& operator>>(B<XStream>& bs, BasicObs<T>& obj)
	{
		bs>>obj.dt>>obj.Q>>obj.Zcount>>obj.ltime>>obj.FreeEnergy>>obj.Qx>>obj.Qy>>obj.Qz>>obj.Q2x>>obj.Q2y>>obj.Q2z;
		return bs;
	}

	//////////////////////////////////////////////////////////////////
	//Observables: Density
	//////////////////////////////////////////////////////////////////

	template <class T>
	friend B<XStream>& operator>>(B<XStream>& bs, Density<T>& obj)
	{
		bs>>obj.rho>>obj.Zcount;
		return bs;
	}

	template <class T>
	friend B<XStream>& operator<<(B<XStream>& bs, Density<T>& obj)
	{
		bs<<obj.rho<<obj.Zcount;
		return bs;
	}

	//////////////////////////////////////////////////////////////////
	//Observables: Qhistogram
	//////////////////////////////////////////////////////////////////

	template <class T>
	friend B<XStream>& operator>>(B<XStream>& bs, Qhistogram<T>& obj)
	{
		bs>>obj.dt>>obj.Q>>obj.ltime>>obj.Qcollection;
		return bs;
	}

	template <class T>
	friend B<XStream>& operator<<(B<XStream>& bs, Qhistogram<T>& obj)
	{
		bs<<obj.dt<<obj.Q<<obj.ltime<<obj.Qcollection;
		return bs;
	}

	//////////////////////////////////////////////////////////////////
	//Observables: Whistogram
	//////////////////////////////////////////////////////////////////

	template <class T>
	friend B<XStream>& operator>>(B<XStream>& bs, Whistogram<T>& obj)
	{
		bs>>obj.dt>>obj.Q>>obj.ltime>>obj.Wcollection;
		return bs;
	}

	template <class T>
	friend B<XStream>& operator<<(B<XStream>& bs, Whistogram<T>& obj)
	{
		bs<<obj.dt<<obj.Q<<obj.ltime<<obj.Wcollection;
		return bs;
	}

	//////////////////////////////////////////////////////////////////
	//Walker
	//////////////////////////////////////////////////////////////////

	template <class T>
	friend B<XStream>& operator<<(B<XStream>& bs, Walker<T>& w)
	{
		int l = w.observablesCollection.size();
		bs<<w.state<<l;
		for(int i=0;i<l;i++)
			bs<<(*w.observablesCollection[i]);
		return bs;
	}

	template <class T>
	friend B<XStream>& operator>>(B<XStream>& bs, Walker<T>& w)
	{
		int l;
		bs>>w.state>>l;
		for(int i=0;i<l;i++)
			bs>>(*w.observablesCollection[i]);
		return bs;
	}

};

#endif /* DMCSNSSERIALIZER_H_ */
