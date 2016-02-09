/*
 * Weight.h
 *
 *  Created on: Nov 2, 2015
 *      Author: Ushnish Ray
 */

#ifndef INCLUDE_WEIGHT_H_
#define INCLUDE_WEIGHT_H_

const double MAXV = 1.0e10;
const double MINV = 1.0e-10;
const int DIVISOR = 1e5;

class Weight
{
private:
	int divisor;
	long double maxvalue,minvalue;

	double val;
	int exponent;

	double initval;
public:

	Weight(const Weight& w)
	{
		divisor = w.divisor;
		maxvalue = w.maxvalue;
		minvalue = w.minvalue;

		val = w.val;
		exponent = w.exponent;
		initval = w.initval;
	}

	Weight(double a = 0.0)
	{
		initval = a;
		maxvalue = MAXV;
		minvalue = MINV;
		divisor = DIVISOR;
		exponent = 0;

		if(fabs(a)>DBL_EPSILON)
		{
			while(fabs(a)>maxvalue)
			{
				a /= divisor;
				exponent++;
			}

			while(fabs(a)<minvalue)
			{
				a *= divisor;
				exponent--;
			}
		}
		else
			a = 0.0;

		val = a;
	}

	Weight(long double maxv, long double minv)
	{
		initval = 0.0;
		val = 0.0;
		exponent = 0;
		divisor = DIVISOR;

		maxvalue = maxv;
		minvalue = minv;
	}

	Weight(long double maxv, long double minv, int _divisor)
	{
		initval = val = 0.0;
		exponent = 0;
		divisor = _divisor;

		maxvalue = maxv;
		minvalue = minv;
	}

	Weight(double a, int _divisor, long double maxv, long double minv)
	{
		initval = a;
		maxvalue = maxv;
		minvalue = minv;
		divisor = _divisor;
		exponent = 0;

		if(fabs(a)>DBL_EPSILON)
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
		
		initval = _val;
		val = _val;
		divisor = _divisor;
		exponent = _exponent;
	}

	long double value()
	{
		return (val*pow(divisor,exponent));
	}

	void resetValue()
	{		
		double a = initval;
		exponent = 0;

		if(fabs(a)>DBL_EPSILON)
		{
			while(fabs(a)>maxvalue)
			{
				a /= divisor;
				exponent++;
			}

			while(fabs(a)<minvalue)
			{
				a *= divisor;
				exponent--;
			}
		}
		else
			a = 0.0;

		val = a;
	}

	void setInitValue(double v)
	{
		initval = v;
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
		//initval = a.initval;
		val = a.val;
		divisor = a.divisor;
		exponent = a.exponent;
		return *this;
	}

	void copy(Weight& w)
	{
		divisor = w.divisor;
		maxvalue = w.maxvalue;
		minvalue = w.minvalue;

		val = w.val;
		exponent = w.exponent;
		initval = w.initval;
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

	Weight operator /(const Weight& a)
	{
		double nv = val/a.val*pow((divisor/a.divisor),a.exponent);
		int e = exponent - a.exponent;

		if(fabs(val)>DBL_EPSILON)
		{
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
		}

		Weight nw = Weight(nv,divisor,e,maxvalue,minvalue);
		return nw;
	}

	void update(double a)
	{
		val *= a;
		if(fabs(val)>DBL_EPSILON)
		{
			//printf("%10.6e %10.6Le %10.6Le %10.6e\n",DBL_EPSILON,minvalue,maxvalue,val);
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

	double logValue()
	{
		return log(val) + exponent*log(divisor);
	}	

	bool zerovalued()
	{
		return (exponent == 0 && fabs(val)<DBL_EPSILON);
	}

	void  add(const Weight& w)
	{
		if(zerovalued())
		{
			val = w.val;
			divisor = w.divisor;
			exponent = w.exponent;
		}
		else
		{
			int maxe = (exponent<w.exponent) ? w.exponent : exponent;
			val = val*pow(divisor,exponent-maxe) + w.val*pow(divisor,w.exponent-maxe);
			exponent = maxe;

			if(fabs(val)>DBL_EPSILON)
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
		fprintf(out,"%10.6e %d %d\n",val,divisor,exponent);
	}
};

#endif /* INCLUDE_WEIGHT_H_ */
