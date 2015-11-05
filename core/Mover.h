/*
 * walker.h
 *
 *  Created on: Aug 4, 2015
 *      Author: ushnish
 *
 Copyright (c) 2015 Ushnish Ray
 All rights reserved.
 */

#ifndef MOVER_H_
#define MOVER_H_

namespace core {

template <class T>
class Mover {
protected:

	double acceptr;

	FILE* debugLog;

public:

	Mover(FILE* _dlog)
	{
		acceptr = 0.0;
		debugLog = _dlog;
	}

	~Mover()
	{
	}

	virtual void initialize(Walker<T>*) = 0;
	virtual void move(Walker<T>*) = 0;

	double getAccept()
	{
		double retval = acceptr;
		acceptr = 0.0;
		return retval;
	}

	void setDebugFile(FILE* _dlog)
	{
		debugLog = _dlog;
	}
};

}
#endif /* MOVER_H_ */
