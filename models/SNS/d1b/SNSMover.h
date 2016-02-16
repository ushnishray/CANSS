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

template <class T>
class SNSMover:public core::Mover<T>
{
protected:

	RunParameters& rp;
	vector<int> pidxavailable;

public:

	SNSMover(FILE* _dlog, RunParameters& _rp):core::Mover<T>(_dlog),rp(_rp)
	{}

	~SNSMover()
	{}

	void initialize(Walker<T>*);
	void move(Walker<T>*);

};
#endif /* MOVER_H_ */
