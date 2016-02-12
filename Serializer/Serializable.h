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
#include "Serializer.h"

class Serializable
{
public:


	void serializeWith(Serializer&) = 0;
};

#endif /* SERIALIZABLE_H_ */
