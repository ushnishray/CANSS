/*
 * Serializer.h
 *
 *  Created on: Feb 11, 2016
 *      Author: ushnish
 *
 	Copyright (c) 2016 Ushnish Ray
	All rights reserved.
 */

#ifndef SERIALIZER_H_
#define SERIALIZER_H_

#include <arpa/inet.h> // for htonl/s, ntoh/s
#include <endian.h> // for htonbe64, if you have it...

template <typename XStream>
class Serializer: public virtual XStream //It is a stream
{
private:
	//uint64_t hton(uint64_t n) { return htonbe64(n); }
	uint32_t mhton(uint32_t n) { return htonl(n); }
	uint16_t mhton(uint16_t n) { return htons(n); }

	// there are no "int" versions - this is ugly but effective...
	uint32_t mhton(int32_t n) { return htonl(n); }
	uint16_t mhton(int16_t n) { return htons(n); }

	//uint64_t ntoh(uint64_t n) { return betoh64(n); }
	uint32_t mntoh(uint32_t n) { return ntohl(n); }
	uint16_t mntoh(uint16_t n) { return ntohs(n); }

	//uint64_t ntoh(uint64_t n) { return betoh64(n); }
	uint32_t mntoh(int32_t n) { return ntohl(n); }
	uint16_t mntoh(int16_t n) { return ntohs(n); }

public:
	typedef Serializer This;

	//Handle primitives
	This& write(const char* s, std::streamsize n)
	{
		XStream::write(s, n);
		return *this;
	}

	This& read(char* s, std::streamsize n)
	{
		XStream::read(s, n);
		return *this;
	}

	template <typename T>
	This& rawwrite(const T& t)
	{
		XStream& p = static_cast<XStream&>(*this);
		size_t n = sizeof(t);
		p<<'['<<n<<']';
		write((const char*)&t,n);
		return *this;
	}

	template <typename T>
	This& rawread(T& t)
	{
		XStream& p = static_cast<XStream&>(*this);
		size_t n;
		char x;
		p >>x>>n>>x;
		read((char*)&t, n);
		return *this;
	}

	template <typename T>
	This& hton(T h)
	{
		T n = mhton(h);
		return rawwrite(n);
	}

	template <typename T>
	This& ntoh(T& h)
	{
		rawread(h);
		h = mntoh(h);
		return *this;
	}

	//Print Functions
	std::string printable(char c)
	{
	    std::ostringstream oss;
	    if (isprint(c))
	        oss << c;
	    else
	        oss << "\\x" << std::hex << std::setw(2) << std::setfill('0')<<(int)(uint8_t)c << std::dec;
	    return oss.str();
	}

	std::string printable(const std::string& s)
	{
	    std::string result;
	    for (std::string::const_iterator i = s.begin(); i != s.end(); ++i)
	        result += printable(*i);
	    return result;
	}

	// conversions for inbuilt & Standard-library types...
/*
	This& operator<<(bool x)
	{
		XStream& bs = *this;
		bs << (x ? 'T' : 'F');
		return bs;
	}

	This& operator>>(bool &x)
	{
		char v;
		static_cast<XStream&>(*this)>>v;
		x = (v == 'T');
		return (*this);
	}
*/

	virtual This& operator<<(int8_t x) { XStream& bs = *this; bs << x; return (*this);}
	virtual This& operator>>(int8_t& x) { XStream& bs = *this; bs >> x; return (*this);}

	virtual This& operator<<(uint8_t x) { XStream& bs = *this; bs << x; return (*this);}
	virtual This& operator>>(uint8_t& x) { XStream& bs = *this; bs >> x; return (*this);}

	virtual This& operator<<(uint16_t x){return (*this).hton(x);}
	virtual This& operator>>(uint16_t& x) { return (*this).ntoh(x);}

	virtual This& operator<<(int16_t x){return (*this).hton(x);}
	virtual This& operator>>(int16_t& x) { return (*this).ntoh(x);}

	virtual This& operator<<(int32_t x){return (*this).hton(x);}
	virtual This& operator>>(int32_t& x) { return (*this).ntoh(x);}

	virtual This& operator<<(uint32_t x){return (*this).hton(x);}
	virtual This& operator>>(uint32_t& x) { return (*this).ntoh(x);}

	virtual This& operator<<(double d)
	{
//		long int i64 = *(reinterpret_cast<long int *>(&d)); /* Ugly, but works */
//		int hiword = (static_cast<int>(i64 >> 32));
//	    int loword = (static_cast<int>(i64));
//
//	    return *this<<hiword<<loword;

	    return (*this).rawwrite(d);
	}

	virtual This& operator>>(double& d)
	{
//		int hiword,loword;
//		*this>>hiword>>loword;
//		long int i64 = static_cast<long int>(hiword);
//		i64 = i64 << 32 | loword;

//		d = *(reinterpret_cast<double *>(&i64));
//		return *this;
		return (*this).rawread(d);
	}

//	virtual This& operator<<(long double d){return (*this).rawwrite(d);}
//	virtual This& operator>>(long double& d) {return (*this).rawread(d);}

	virtual This& operator<<(const std::string& x)
	{
		unsigned int l = x.size();
		(*this) << l;
		return (*this).write(x.data(), x.size());
	}

	virtual This& operator>>(std::string& x)
	{
		unsigned int size;
		(*this) >> size;
		char* data = new char[size+1];
		This& rv = (*this).read(data,size);
		data[size] = 0;
		x = string(data);
		delete[] data;
		return rv;
	}

	virtual This& operator<<(const char* x)
	{
		unsigned int l = strlen(x);
		(*this)<<l;
		return (*this).write(x,l);
	}

	virtual This& operator>>(char* x)
	{
		unsigned int l;
		(*this)>>l;
		return (*this).read(x,l);
	}

	//****************************************************************

	//////////////////////////////////////////////////////////////////
	//Vect
	//////////////////////////////////////////////////////////////////

	template <class T>
	This& operator<<(vect<T>& obj)
	{
		*this<<obj.x<<obj.y<<obj.z;
		return *this;
	}

	template <class T>
	This& operator>>(vect<T>& obj)
	{
		*this>>obj.x>>obj.y>>obj.z;
		return *this;
	}

	//////////////////////////////////////////////////////////////////
	//Particle Map
	//////////////////////////////////////////////////////////////////

	/*
	template <class T>
	This& operator<<( PtclMap<T>& obj)
	{
		Serializer<XStream>& bs = *this;
		unsigned int size = obj.size();
		bs<<size;
		typename PtclMap<T>::iterator it;
		for(it=obj.begin();it!=obj.end();++it)
		{
			vect<T> a = it->first;
			(*this)<<a;
			bs<<it->second;
		}
		return *this;
	}

	template <class T>
	This& operator>>( PtclMap<T>& obj)
	{
		Serializer<XStream>& bs = *this;
		unsigned int size;
		bs>>size;
		for(unsigned int i = 0;i<size;i++)
		{
			vect<T> v;
			int val;
			(*this)>>v;
			bs>>val;
			obj[v] = val;
		}

		return *this;
	}*/

	template <class T>
	This& operator<<( PtclMap<T>* obj)
	{
		unsigned int tsize = obj->size();
		(*this)<<tsize;

		typename PtclMap<T>::iterator it;
		for(it=obj->begin();it!=obj->end();++it)
		{
			vect<T> v(it->first.x,it->first.y,it->first.z);
			(*this)<<v<<it->second;
		}
		return *this;
	}

	template <class T>
	This& operator>>( PtclMap<T>* obj)
	{
		unsigned int tsize;
		(*this)>>tsize;

		for(unsigned int i = 0;i<tsize;i++)
		{
			vect<T> v;
			int val;
			(*this)>>v>>val;
			(*obj)[v] = val;
		}

		return *this;
	}

	//////////////////////////////////////////////////////////////////
	//Vector
	//////////////////////////////////////////////////////////////////

	template <class T>
	This& operator<<( vector<T>& obj)
	{
		Serializer<XStream>& bs = *this;
		bs<<(unsigned int) obj.size();
		typename vector<T>::iterator it;
		for(it=obj.begin();it!=obj.end();++it)
			bs<<(*it);

		return bs;
	}

	template <class T>
	This& operator>>( vector<T>& obj)
	{
		This& bs = *this;
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
	This& operator<<( vectToValue<T>& obj)
	{
		Serializer<XStream>& bs = *this;
		bs<<(unsigned int) obj.size();

		typename vectToValue<T>::iterator it;
		for(it=obj.begin();it!=obj.end();++it)
		{
			vect<T> a = it->first;
			(*this)<<a;
			bs<<it->second;
		}
		return bs;
	}

	template <class T>
	This& operator>>( vectToValue<T>& obj)
	{
		Serializer<XStream>& bs = *this;
		int tsize;
		bs>>tsize;
		for(int i = 0;i<tsize;i++)
		{
			vect<T> v;
			double val;
			(*this)>>v;
			bs>>val;
			obj[v] = val;
		}
		return bs;
	}

	//////////////////////////////////////////////////////////////////
	//Weight
	//////////////////////////////////////////////////////////////////

	This& operator<<(Weight& w)
	{
		(*this)<<w.divisor<<w.exponent<<w.initval<<w.maxvalue<<w.minvalue<<w.val;
		return *this;
	}

	This& operator>>(Weight& w)
	{
		(*this)>>w.divisor>>w.exponent>>w.initval>>w.maxvalue>>w.minvalue>>w.val;
		return *this;
	}


#if 0
	//////////////////////////////////////////////////////////////////
	//WalkerState
	//////////////////////////////////////////////////////////////////

	template <class T>
	This& operator<<( WalkerState<T>& obj)
	{
		//(*this)<<obj.DIM<<obj.dQ<<obj.ltime<<obj.particleCount<<obj.weight<<(obj.Rcurr);
		//(*this)<<obj.DIM<<obj.dQ<<obj.ltime<<obj.particleCount<<obj.weight;
		(*this)<<obj.DIM<<obj.dQ<<obj.ltime<<obj.particleCount<<(obj.Rcurr)<<obj.weight;
		return (*this);
	}

	template <class T>
	This& operator>>( WalkerState<T>& obj)
	{
		obj.Rcurr->clear();
		//(*this)>>obj.DIM>>obj.dQ>>obj.ltime>>obj.particleCount>>obj.weight>>(obj.Rcurr);
		//(*this)>>obj.DIM>>obj.dQ>>obj.ltime>>obj.particleCount>>obj.weight;
		(*this)>>obj.DIM>>obj.dQ>>obj.ltime>>obj.particleCount>>(obj.Rcurr)>>obj.weight;
		return (*this);
	}

	//////////////////////////////////////////////////////////////////
	//Observables: BasicObs
	//////////////////////////////////////////////////////////////////

	template <class T>
	This& operator<<( BasicObs<T>& obj)
	{
		(*this)<<obj.dt<<obj.Zcount<<obj.ltime<<obj.Q
				<<obj.freeEnergy<<obj.Qx<<obj.Qy<<obj.Qz<<obj.Q2x<<obj.Q2y<<obj.Q2z;
		return *this;
	}

	template <class T>
	This& operator>>( BasicObs<T>& obj)
	{
		(*this)>>obj.dt>>obj.Zcount>>obj.ltime>>obj.Q
				>>obj.freeEnergy>>obj.Qx>>obj.Qy>>obj.Qz>>obj.Q2x>>obj.Q2y>>obj.Q2z;
		return *this;
	}

	//////////////////////////////////////////////////////////////////
	//Observables: Density
	//////////////////////////////////////////////////////////////////

	template <class T>
	This& operator>>( Density<T>& obj)
	{
		Serializer<XStream>& bs = *this;
		(*this)>>obj.rho;
		bs>>obj.Zcount;
		return *this;
	}

	template <class T>
	This& operator<<( Density<T>& obj)
	{
		Serializer<XStream>& bs = *this;
		(*this)<<obj.rho;
		bs<<obj.Zcount;
		return *this;
	}

	//////////////////////////////////////////////////////////////////
	//Observables: Qhistogram
	//////////////////////////////////////////////////////////////////

	template <class T>
	This& operator>>( Qhistogram<T>& obj)
	{
		Serializer<XStream>& bs = *this;
		bs>>obj.dt>>obj.ltime;
		(*this)>>obj.Q>>obj.Qcollection;
		return *this;
	}

	template <class T>
	This& operator<<( Qhistogram<T>& obj)
	{
		Serializer<XStream>& bs = *this;
		bs<<obj.dt<<obj.ltime;
		(*this)<<obj.Q<<obj.Qcollection;
		return *this;
	}

	//////////////////////////////////////////////////////////////////
	//Observables: Whistogram
	//////////////////////////////////////////////////////////////////

	template <class T>
	This& operator>>( Whistogram<T>& obj)
	{
		(*this)>>obj.dt>>obj.localWeight>>obj.ltime>>obj.Wcollection;
		return *this;
	}

	template <class T>
	This& operator<<( Whistogram<T>& obj)
	{
		(*this)<<obj.dt<<obj.localWeight<<obj.ltime<<obj.Wcollection;
		return *this;
	}
#endif

	//////////////////////////////////////////////////////////////////
	//Walker
	//////////////////////////////////////////////////////////////////
#if 0
	template <class T>
	This& operator<<( Walker<T>& w)
	{
		int l = w.observablesCollection.size();
		(*this)<<w.state<<l;

		//fprintf(w.state.out,"State and obs size %d\n",l);
		//fflush(w.state.out);
		//w.state.display();
		//fflush(w.state.out);

		for(int i=0;i<l;i++)
		{
			measures::Observable<T>* ob = w.observablesCollection[i];
			if( BasicObs<T>* obj = dynamic_cast<BasicObs<T>*>(ob))
				(*this)<<"BasicObs"<<obj->baseFileName<<*obj;
			else if( Density<T>* obj = dynamic_cast<Density<T>*>(ob))
				(*this)<<"Density"<<obj->baseFileName<<*obj;
			else if( Qhistogram<T>* obj = dynamic_cast<Qhistogram<T>*>(ob))
				(*this)<<"Qhistogram"<<obj->baseFileName<<*obj;
			else if( Whistogram<T>* obj = dynamic_cast<Whistogram<T>*>(ob))
				(*this)<<"Whistogram"<<obj->baseFileName<<*obj;
		}
		return *this;
	}

	template <class T>
	This& operator>>( Walker<T>& w)
	{
		//copy don't add

		vector<Observable<T>*>* localObs = new vector<Observable<T>*>;

		int l;
		(*this)>>w.state>>l;
		//fprintf(w.state.out,"State and obs size %d\n",l);
		//w.state.display();
		//fflush(w.state.out);

		for(int i=0;i<l;i++)
		{
			string type,bfname;
			(*this)>>type>>bfname;
			if(type.compare("BasicObs") == 0)
			{
				BasicObs<T>& obj = *(new BasicObs<T>(w.state,bfname,w.state.out));
				(*this)>>obj;
				localObs->push_back(&obj);
			}
			else if(type.compare("Density") == 0)
			{
				Density<T>& obj = *(new Density<T>(w.state,bfname,w.state.out));
				(*this)>>obj;
				localObs->push_back(&obj);
			}
			else if(type.compare("Qhistogram") == 0)
			{
				Qhistogram<T>& obj = *(new Qhistogram<T>(w.state,bfname,w.state.out));
				(*this)>>obj;
				localObs->push_back(&obj);
			}
			else if(type.compare("Whistogram") == 0)
			{
				Whistogram<T>& obj = *(new Whistogram<T>(w.state,bfname,w.state.out));
				(*this)>>obj;
				localObs->push_back(&obj);
			}
		}

		//Now copy i.e replace current data
		for(int i=0;i<w.observablesCollection.size();i++)
			w.observablesCollection[i]->copy((*localObs)[i]);
		delete localObs;

		return *this;
	}
#endif

	/*
	template <class T, class U>
	This& operator<<( core::Walker<T,U>& w)
	{
		w.serialize(*this);
	}

	template <class T, class U>
	This& operator>>( core::Walker<T,U>& w)
	{
		w.unserialize(*this);
	}
	*/

};

#endif /* SERIALIZER_H_ */
