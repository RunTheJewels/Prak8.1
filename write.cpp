#include <iostream>
#include <fstream>
#include <cstdlib>


using namespace std;

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		cout << "Usage: ./a.out [size of system]";
		return -1;
	}
	int n = atoi(argv[1]);
	float a;
	ofstream f("matrixn.txt", ofstream::binary);
	f.write((char*)&n, sizeof(int));
	srand(time(NULL));
	for (int i = 0; i < n * (n+1); i++)
	{
		a = (float) (rand()) / (float) (RAND_MAX) * 200.0 - 100.0;
		f.write((char*)&a, sizeof(float));
	}
	f.close();
	return 0;
}
