#include <iostream>
#include <vector>

using namespace std;

int main(int argc,char* argv[])
{
	int procCount = 11;
	//         0,1,2,3, 4,5,6, 7,8, 9,10		
	int A[] = {0,4,2,0,-2,0,3,-2,1,-4,-2};

	//Accummulate senders and receivers
	vector<int> sendProcs, sendProcCount, recvProcs, recvProcCount;
	for(int i = 1;i<procCount;i++)
	{
		int rem = A[i];
		//MPI_Recv(&rem,1,MPI_INT,i,tag,MPI_COMM_WORLD,&stat);
		
		if(rem>0)
		{
			sendProcs.push_back(i);
			sendProcCount.push_back(rem);
		}
		else if(rem<0)
		{
			recvProcs.push_back(i);
			recvProcCount.push_back(-rem);
		}
	}

	//Assemble messages for senders and receivers
	vector<int>* mproc = new vector<int>[sendProcs.size()+recvProcs.size()];
	vector<int>* mcount = new vector<int>[sendProcs.size()+recvProcs.size()];
	int sendp = 0,recvp = 0;

	while(sendp<sendProcs.size() && recvp<recvProcs.size())
	{
		if(recvProcCount[recvp]<sendProcCount[sendp])
		{
			//Assign receivers which processors to expect data from and data count
			mproc[recvp+sendProcs.size()].push_back(sendProcs[sendp]);
			mcount[recvp+sendProcs.size()].push_back(recvProcCount[recvp]);

			//Assign senders which processors to send data to and data count
			mproc[sendp].push_back(recvProcs[recvp]);
			mcount[sendp].push_back(recvProcCount[recvp]);

			sendProcCount[sendp] -= recvProcCount[recvp];
			recvp++;
		}
		else
		{
			//Assign receivers which processors to expect data from and data count
			mproc[recvp+sendProcs.size()].push_back(sendProcs[sendp]);
			mcount[recvp+sendProcs.size()].push_back(sendProcCount[sendp]);

			//Assign senders which processors to send data to and data count
			mproc[sendp].push_back(recvProcs[recvp]);
			mcount[sendp].push_back(sendProcCount[sendp]);

			recvProcCount[recvp] -= sendProcCount[sendp];
			//sendProcCount[sendp] = 0;

			sendp++;
			if(recvProcCount[recvp]==0)
				recvp++;
		}
	}

	for(int i = 0;i<sendProcs.size()+recvProcs.size();i++)
	{
		if(i<sendProcs.size())
		{
			cout<<"Processor: "<<sendProcs[i]<<endl;
			cout<<"Sender\n";
		}
		else
		{
			cout<<"Processor: "<<recvProcs[i-sendProcs.size()]<<endl;
			cout<<"Receiver\n";
		}

		for(int j=0;j<mproc[i].size();j++)
			cout<<"With: "<<mproc[i][j]<<" count: "<<mcount[i][j]<<endl;				
	}
	
	return 0;
}
