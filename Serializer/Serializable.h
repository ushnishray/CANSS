/*
 * Serializable.h
 *
 *  Created on: Feb 11, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.

	This software is cannot be used for commercial purposes in any way whatsoever.
 */

#ifndef SERIALIZABLE_H_
#define SERIALIZABLE_H_

template <typename U>
class Serializable
{
public:

	virtual void serialize(Serializer<U>&) = 0;

	virtual void unserialize(Serializer<U>&) = 0;
};

#endif /* SERIALIZABLE_H_ */
