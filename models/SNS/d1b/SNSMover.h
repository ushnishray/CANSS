/*
 * walker.h
 *
 *  Created on: Aug 4, 2015
 *      Author: ushnish
 *
 Copyright (c) 2015 Ushnish Ray
 All rights reserved.
 */

#ifndef SNSMOVER_H_
#define SNSMOVER_H_

template <class T, class U>
class SNSMover:public core::Mover<T,U>
{
protected:

	RunParameters& rp;
	vector<int> pidxavailable;

public:

	SNSMover(FILE* _dlog, RunParameters& _rp):core::Mover<T,U>(_dlog),rp(_rp)
	{}

	~SNSMover()
	{}

	void initialize(Walker<T,U>*);
	void move(Walker<T,U>*);

};
#endif /* MOVER_H_ */
