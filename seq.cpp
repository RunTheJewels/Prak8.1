#include <iostream>
#include <fstream>
#include <math.h>


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
	ifstream file("matrixn.txt", ifstream::binary);
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


	// Основной цикл
	for (int i = 0; i < n-1; i++)
	{
		float x[n-i];
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
		for (int j = i; j < n; j++)
		{
			for (int k = i; k < n; k++)
			{
				tmpA[j][k] = A[j][k];
				for (int l = i; l < n; l++)
				{
					tmpA[j][k] -= 2*x[j-i]*x[l-i]*A[l][k];
				}
			}
		}
		for (int j = i; j < n; j++)
		{
			tmpb[j] = b[j];
			for (int k = i; k < n; k++)
			{
				tmpb[j] -= 2*x[j-i]*x[k-i]*b[k];
			}
		}
		for (int j = i; j < n; j++)
		{
			for (int k = i; k < n; k++)
			{
				A[j][k] = tmpA[j][k];
			}
		}
		for (int j = i; j < n; j++)
		{
			b[j] = tmpb[j];
		}
		cout << endl << "------" << i << "------" << endl;
		show_system(n,A,b);
	}
	
	// Обратный ход
	for (int i = n-1; i > 0; i--)
	{
		for (int j = i-1; j >= 0; j--)
		{
			float coef = A[j][i] / A[i][i];
			A[j][i] -= A[i][i] * coef;
			b[j] -= b[i] * coef;
		}
		cout << endl << "------ Обратный" << i << "------" << endl;
		show_system(n,A,b);
	}

	for (int i = 0; i < n; i++)
	{
		cout << "x" << i << " = " << b[i]/A[i][i] << endl;
	}

	return 0;
}

