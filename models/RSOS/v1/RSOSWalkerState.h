/*
 * RSOSWalkerState.h
 *
 *  Created on: Mar 22, 2016
 *      Author: ushnishray
 */

#ifndef MODELS_RSOS_V1_RSOSWALKERSTATE_H_
#define MODELS_RSOS_V1_RSOSWALKERSTATE_H_

template <class T, class U>
class RSOSWalkerState:public core::WalkerState<T,U>
{
public:
	int L;
	int* height;

	RSOSWalkerState(int _dim, Weight& _w, FILE* ot, int _L):WalkerState<T,U>(_dim,_w,ot)
	{
		L = _L;
		height = new int[L];
		for(int i = 0;i<L;i++)
			height[i] = 0;
	}

	RSOSWalkerState(int _dim, Weight& _w, FILE* ot, int _L, int wid):WalkerState<T,U>(_dim,_w,ot,wid)
	{
		L = _L;
		height = new int[L];
		for(int i = 0;i<L;i++)
			height[i] = 0;
	}

	RSOSWalkerState(WalkerState<T,U>& a, int _L)
	{
		core::WalkerState<T,U>::copy(a);
		L = _L;
		height = new int[L];
		for(int i = 0;i<L;i++)
			height[i] = 0;
	}

	~RSOSWalkerState()
	{
		delete[] height;
	}

	void copy(RSOSWalkerState<T,U>& w)
	{
		core::WalkerState<T,U>::copy((WalkerState<T,U>) w);
		//this->L = w.L;
		//delete[] height;
		//height = new int[L];
		memcpy(height,w.height,sizeof(int)*L);
	}

	RSOSWalkerState* duplicate()
	{
		RSOSWalkerState* a = new RSOSWalkerState(*this,this->L);
		memcpy(a,this->height,sizeof(int)*L);
		return a;
	}

	void reset()
	{
		this->dQ.x = (T) 0.0;
		this->dQ.y = (T) 0.0;
		this->dQ.z = (T) 0.0;

		this->ltime = 0;
		this->dweight = 0.0;
		this->weight.resetValue();
	}

	virtual void serialize(Serializer<U>& obj)
	{
		core::WalkerState<T,U>::serialize(obj);
		obj<<L;
		for(int i = 0;i<L;i++)
			obj<<height[i];
	}

	virtual void unserialize(Serializer<U>& obj)
	{
		core::WalkerState<T,U>::unserialize(obj);
		obj>>L;
		for(int i = 0;i<L;i++)
			obj>>height[i];
	}
};

#endif /* MODELS_RSOS_V1_RSOSWALKERSTATE_H_ */
