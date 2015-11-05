/*
 * Weight.h
 *
 *  Created on: Nov 2, 2015
 *      Author: Ushnish Ray
 */

#ifndef INCLUDE_WEIGHT_H_
#define INCLUDE_WEIGHT_H_


class Weight
{
private:
	int divisor;
	long double maxvalue,minvalue;

	double val;
	int exponent;

public:

	Weight(const Weight& w)
	{
		divisor = w.divisor;
		maxvalue = w.maxvalue;
		minvalue = w.minvalue;

		val = w.val;
		exponent = w.exponent;
	}

	Weight(long double maxv, long double minv)
	{
		val = 0.0;
		exponent = 0;
		divisor = DIVISOR;

		maxvalue = maxv;
		minvalue = minv;
	}

	Weight(long double maxv, long double minv, int _divisor)
	{
		val = 0.0;
		exponent = 0;
		divisor = _divisor;

		maxvalue = maxv;
		minvalue = minv;
	}

	Weight(long double a, int _divisor, long double maxv, long double minv)
	{
		maxvalue = maxv;
		minvalue = minv;
		divisor = _divisor;
		exponent = 0;

		if(fabs(a)>ZEROTOL)
		{
			while(fabs(a)>maxv)
			{
				a /= divisor;
				exponent++;
			}

			while(fabs(a)<minv)
			{
				a *= divisor;
				exponent--;
			}
		}
		else
			a = 0.0;

		val = a;
	}

	Weight(double _val, int _divisor, int _exponent, long double maxv, long double minv)
	{
		maxvalue = maxv;
		minvalue = minv;

		val = _val;
		divisor = _divisor;
		exponent = _exponent;
	}

	long double value()
	{
		return (val*pow(divisor,exponent));
	}

	void setValue(double v,int e)
	{
		val = v;
		exponent = e;
	}

	void setBounds(long double _maxvalue,long double _minvalue)
	{
		maxvalue = _maxvalue;
		minvalue = _minvalue;
	}

	Weight& operator =(const Weight& a)
	{
		val = a.val;
		divisor = a.divisor;
		exponent = a.exponent;
		return *this;
	}

	bool operator ==(const Weight& a)
	{
		return (val == a.val && divisor == a.divisor && exponent == a.exponent);
	}

	Weight operator *(const Weight& a)
	{
		double nv = val*a.val*pow((a.divisor/divisor),a.exponent);
		int e = exponent + a.exponent;

		while(fabs(nv)>maxvalue)
		{
			nv /= divisor;
			e++;
		}

		while(fabs(nv)<minvalue)
		{
			nv *= divisor;
			e--;
		}

		Weight nw = Weight(nv,divisor,e,maxvalue,minvalue);
		return nw;
	}

	void update(double a)
	{
		val *= a;
		if(fabs(val)>ZEROTOL)
		{
			while(fabs(val)>maxvalue)
			{
				val /= divisor;
				exponent++;
			}

			while(fabs(a)<minvalue)
			{
				val *= divisor;
				exponent--;
			}
		}
	}

	double logValue()
	{
		return log(val) + exponent*log(divisor);
	}

	void  add(const Weight& w)
	{
		int maxe = (exponent<w.exponent) ? w.exponent : exponent;
		val = val*pow(divisor,exponent-maxe) + w.val*pow(divisor,w.exponent-maxe);
		exponent = maxe;

		if(fabs(val)>ZEROTOL)
		{
			while(fabs(val)>maxvalue)
			{
				val /= divisor;
				exponent++;
			}

			while(fabs(val)<minvalue)
			{
				val *= divisor;
				exponent--;
			}
		}
	}

	void mpiReceive(int from)
	{
		int tag,recv;
		MPI_Status stat;

		double mval;
		int mexponent;

		MPI_Recv(&mval,1,MPI_DOUBLE,from,tag,MPI_COMM_WORLD,&stat);
		MPI_Recv(&mexponent,1,MPI_INT,from,tag,MPI_COMM_WORLD,&stat);

		Weight *w = new Weight(mval,divisor,mexponent,maxvalue,minvalue);
		this->add(*w);
		delete w;
	}

	void mpiSend(int to)
	{
		int tag,recv;
		MPI_Status stat;

		MPI_Send(&this->val,1,MPI_DOUBLE,to,tag,MPI_COMM_WORLD);
		MPI_Send(&this->exponent,1,MPI_INT,to,tag,MPI_COMM_WORLD);
	}

	void display(FILE* out)
	{
		fprintf(out,"%lf %d %d",val,divisor,exponent);
	}
};

#endif /* INCLUDE_WEIGHT_H_ */
