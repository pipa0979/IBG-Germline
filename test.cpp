#include<iostream>
#include<fstream>
using namespace std;

int main()
{

	cout<<"please input the file name"<<endl;
	string s,line;
	int count;
	cin>>s;
	ifstream stream;
	stream.open(s.c_str());
	int a=stream.tellg();
	while(getline(stream,line))
	{
		stream.seekg(a);
		cout<<a<<endl;
		count++;
		getline(stream,line);
		a=stream.tellg();

	}

	cout<<count<<endl;

}



