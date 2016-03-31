/*
 * walker.h
 *
 *  Created on: Aug 4, 2015
 *      Author: ushnish
 *
 Copyright (c) 2015 Ushnish Ray
 All rights reserved.
 */

#ifndef RSOSMOVER_H_
#define RSOSMOVER_H_

template <class T, class U>
class RSOSMover:public core::Mover<T,U>
{
protected:

	RunParameters& rp;
	unsigned int pidx;
	float amp;
public:

	RSOSMover(FILE* _dlog, RunParameters& _rp):core::Mover<T,U>(_dlog),rp(_rp)
	{
		pidx = 0;
		amp = 2.0*M_PI;
	}

	~RSOSMover()
	{}

	void initialize(Walker<T,U>*);
	void move(Walker<T,U>*);

};
#endif /* MOVER_H_ */
