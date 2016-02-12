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

// uint64_t hton(uint64_t n) { return htonbe64(n); }
uint32_t hton(uint32_t n) { return htonl(n); }
uint16_t hton(uint16_t n) { return htons(n); }

// there are no "int" versions - this is ugly but effective...
uint32_t hton(int32_t n) { return htonl(n); }
uint16_t hton(int16_t n) { return htons(n); }

// uint64_t ntoh(uint64_t n) { return betoh64(n); }
uint32_t ntoh(uint32_t n) { return ntohl(n); }
uint16_t ntoh(uint16_t n) { return ntohs(n); }

// uint64_t ntoh(uint64_t n) { return betoh64(n); }
uint32_t ntoh(int32_t n) { return ntohl(n); }
uint16_t ntoh(int16_t n) { return ntohs(n); }

template <typename XStream>
class Serializer : public XStream
{
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
		static_cast<XStream&>(*this) << '[' << sizeof t << ']';
		return write((const char*)&t, sizeof t);
	}

	template <typename T>
	This& rawread(T& t)
	{
		char x; int ss;
		static_cast<XStream&>(*this) >> x >> ss >> x;
		return read((char*)&t, ss);
	}

	template <typename T>
	This& hton(T h)
	{
		T n = ::hton(h);
		return rawwrite(n);
	}

	template <typename T>
	This& ntoh(T& h)
	{
		rawread(h);
		T v = h;
		h = ::ntoh(v);
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
	friend This& operator<<(This& bs, bool x) { return bs << (x ? 'T' : 'F'); }

	friend This& operator<<(This& bs, int8_t x) { return bs << x; }
	friend This& operator>>(This& bs, int8_t& x) { return bs >> x;}

	friend This& operator<<(This& bs, uint8_t x) { return bs << x; }
	friend This& operator>>(This& bs, uint8_t& x) { return bs >> x;}

	friend This& operator<<(This& bs, uint16_t x){return bs.hton(x);}
	friend This& operator>>(This& bs, uint16_t& x) { return bs.ntoh(x);}

	friend This& operator<<(This& bs, int16_t x){return bs.hton(x);}
	friend This& operator>>(This& bs, int16_t& x) { return bs.ntoh(x);}

	friend This& operator<<(This& bs, int32_t x){return bs.hton(x);}
	friend This& operator>>(This& bs, int32_t& x) { return bs.ntoh(x);}

	friend This& operator<<(This& bs, uint32_t x){return bs.hton(x);}
	friend This& operator>>(This& bs, uint32_t& x) { return bs.ntoh(x);}

	friend This& operator<<(This& bs, double d){return bs.rawwrite(d);}
	friend This& operator>>(This& bs, double& d) {return bs.rawread(d);}

	friend This& operator<<(This& bs, const std::string& x)
	{
		bs << x.size();
		return bs.write(x.data(), x.size());
	}

	friend This& operator>>(This& bs, std::string& x)
	{
		int size;
		bs >> size;
		char* data = new char[size+1];
		This& rv = bs.read(data,size);
		data[size] = 0;
		x = string(data);
		delete[] data;
		return rv;
	}

	friend This& operator<<(This& bs, const char* x)
	{
		int l = strlen(x);
		bs<<l;
		return bs.write(x,l);
	}
	friend This& operator>>(This& bs, char* x)
	{
		int l;
		bs>>l;
		return bs.read(x,l);
	}


};

#endif /* SERIALIZER_H_ */
