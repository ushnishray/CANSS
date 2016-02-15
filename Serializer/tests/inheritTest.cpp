#include <iostream>

using namespace std;

class Weight
{
template <typename U>
friend class S;

private:
	double val;

public:
	Weight() {val = 1.0;}

};

template <typename U>
class S
{
public:
	void fa(Weight& x)
	{
		cout<<x.val<<endl;	
	}
};

int main(int argc, char* argv[])
{
	Weight x;
	S<int> b;

	b.fa(x);
}
