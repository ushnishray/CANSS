#include <iostream>

using namespace std;

struct Enc
{
	int& a;
	Enc(int &b):a(b) {}
};


class Walker
{
	Enc& y;
public:
	
	Walker(int& b):y(*(new Enc(b)))
	{
	}

	void display()
	{	
		cout<<y.a<<endl;
	}

};

int main(int argc,char *argv[])
{
	int& x = *(new int);

	x = 50;
	cout<<x<<endl;

	delete x;
}
