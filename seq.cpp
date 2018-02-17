#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>


const bool DEBUG = false;

const bool timer = true;

using namespace std;

void show_system(int n, float** A, float b[])
{
	for (int i = 0; i < n; i++) 
	{
		for (int j = 0; j < n; j++)
		{
			cout << (j!=0 and A[i][j]>=0?"+":"") << A[i][j] << "*x" << j << " ";
		}
		cout << "= " << b[i] << endl;
	}
}

int main(int argc, char** argv)
{
	// Инициализация системы
	char* filename;
	if (argc > 1)
	{
		filename = argv[1];
	}
	else
	{
		filename = new char[12];
		strcpy(filename, "matrixn.txt");
	}
	ifstream file(filename, ifstream::binary);
	int n;
	file.read((char*)&n, sizeof(int));
	float **A = new float*[n];
	for (int i = 0; i < n; i++)
	{
		A[i] = new float[n];
	}
	float tmpA[n][n];
	float b[n];
	float tmpb[n];
	for (int i = 0; i < n; i++) 
	{
		for (int j = 0; j < n; j++)
		{
			file.read((char*)&A[i][j], sizeof(float));
		}
	}
	for (int i = 0; i < n; i++)
	{
		file.read((char*)&b[i], sizeof(float));
	}
	file.close();

	show_system(n,A,b);

	clock_t t1_start, t1_end, t2_start, t2_end;
	t1_start = clock();
	// Основной цикл
	for (int i = 0; i < n-1; i++)
	{
		float x[n-i];
		float xA[n-i];
		float x_norm = 0;
		for (int j = 0; j < n-i; j++)
		{
			x[j] = A[j+i][i];
			x_norm += x[j]*x[j];
		}
		x_norm = sqrt(x_norm);
		x[0] -= x_norm;
		x_norm = 0;
		for (int j = 0; j < n-i; j++)
		{
			x_norm += x[j]*x[j];
		}
		x_norm = sqrt(x_norm);
		for (int j = 0; j < n-i; j++)
		{
			x[j] /= x_norm;
		}
		
		/*for (int j = i; j < n; j++)
		{
			for (int k = i; k < n; k++)
			{
				tmpA[j][k] = A[j][k];
				for (int l = i; l < n; l++)
				{
					tmpA[j][k] -= 2*x[j-i]*x[l-i]*A[l][k];
				}
			}
		}*/
		for (int j = i; j < n; j++)
		{
			xA[j-i] = 0;
			for (int k = i; k < n; k++)
			{
				xA[j-i] += x[k-i]*A[k][j];
			}
		}
		for (int j = i; j < n; j++)
		{
			for (int k = i; k < n; k++)
			{
				A[j][k] -= 2*x[j-i]*xA[k-i];
			}
		}
		/*for (int j = i; j < n; j++)
		{
			tmpb[j] = b[j];
			for (int k = i; k < n; k++)
			{
				tmpb[j] -= 2*x[j-i]*x[k-i]*b[k];
			}
		}*/
		float bx = 0;
		for (int j = i; j < n; j++)
		{
			bx += x[j-i]*b[j];
		}
		for (int j = i; j < n; j++)
		{
			b[j] -= 2*x[j-i]*bx;
		}
		/*for (int j = i; j < n; j++)
		{
			for (int k = i; k < n; k++)
			{
				A[j][k] = tmpA[j][k];
			}
		}
		for (int j = i; j < n; j++)
		{
			b[j] = tmpb[j];
		}*/
		if (DEBUG)
		{
			cout << endl << "------" << i << "------" << endl;
			show_system(n,A,b);
		}
	}
	t1_end = clock();
	t2_start = clock();
	// Обратный ход
	float ans[n];
	for (int i = n-1; i >=0; i--)
	{
		float right = 0;
		for (int j = n-1; j > i; j--)
		{
			right += A[i][j] * ans[j];
		}
		ans[i] = (b[i] - right) / A[i][i];
	}
	t2_end = clock();
	if (timer)
	for (int i = 0; i < n; i++)
	{
		cout << "x" << i << " = " << ans[i] << endl;
	}

	ofstream stats("stats", ofstream::app);
	stats << 1 << " " << n << " " << (float)(t1_end-t1_start)/CLOCKS_PER_SEC << " " << (float)(t2_end-t2_start)/CLOCKS_PER_SEC << endl;
	stats.close();
	return 0;
}

