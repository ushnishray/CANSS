#include <iostream>
#include <sstream>

using namespace std;

int main(int argc,char* argv[])
{
	stringstream s;
	s<<'['<<sizeof(double)<<']';
	double x = 3.142;
	s.write((const char*)&x,sizeof x);

	char c;	
	int a;
	s>>c>>a>>c;

	double y;
	s.read((char*)&y,sizeof y);
	cout<<"Output: "<<a<<" "<<y<<"\n";	


	long f = *(reinterpret_cast<long *>(&x));
	int high = static_cast<int>(f >> 32);
	int low = static_cast<int>(f);

	cout<<low<<" "<<high<<"\n";

	long g = (static_cast<long>(high) << 32) ;
	long q = (static_cast<long>(low) << 32) ; q = q >> 32;
	g = g | q;
	double h = *(reinterpret_cast<double *>(&g));
	cout<<h<<"\n";
	
	return 0;
}
