#include <iostream>
#include <fstream>
#include <cstdlib>


using namespace std;

int main(int argc, char** argv)
{
	int n;
	float a;
	ifstream f("matrixn.txt", ifstream::binary);
	f.read((char*)&n, sizeof(int));
	cout << n << endl;
	for (int i = 0; i < n * (n+1); i++)
	{
		f.read((char*)&a, sizeof(float));
		cout << a << " ";
		if (i % n == n-1) cout << endl;
	}
	f.close();
	return 0;
}

