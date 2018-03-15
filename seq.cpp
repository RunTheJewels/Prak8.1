#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>

const bool residual = true; // Считать ли невязку

const bool func_fill = false; // Заполнять таблицу функцией fill или же по файлу

const bool DEBUG = false; 

const bool timer = true; // Писать ли время выполнения в stats

using namespace std;

float fill(int i, int j)
{
	return i>j?i:j;
}



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
	int n;
	float **A;
	float *b;
	// Инициализация системы
	if (func_fill)
	{
		int n;
		cout << "Размер генерируемой системы: ";
		cin >> n;
		A = new float*[n];
		for (int i = 0; i < n; i++)
		{
			A[i] = new float[n];
		}
		b = new float[n];
		for (int i = 0; i < n; i++)
		{
			b[i] = 0;
			for (int j = 0; j < n; j++)
			{
				A[i][j] = fill(i,j);
				if (j % 2 == 0) b[i] += A[i][j];
			}
		}
	} else
	{
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
		n;
		file.read((char*)&n, sizeof(int));
		A = new float*[n];
		for (int i = 0; i < n; i++)
		{
			A[i] = new float[n];
		}
		b = new float[n];
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
	}
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
		float bx = 0;
		for (int j = i; j < n; j++)
		{
			bx += x[j-i]*b[j];
		}
		for (int j = i; j < n; j++)
		{
			b[j] -= 2*x[j-i]*bx;
		}
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

	if (residual)
	{
		float r = 0;
		float rtmp;
		for (int i = 0; i < n; i++)
		{
			rtmp = b[i];
			for (int j = 0; j < n; j++)
			{
				rtmp -= A[i][j] * ans[j];
			}
			r += rtmp * rtmp;
		}
		cout << "Невязка: " << sqrt(r) << endl;
	}	

	ofstream stats("stats", ofstream::app);
	stats << 1 << " " << n << " " << (float)(t1_end-t1_start)/CLOCKS_PER_SEC << " " << (float)(t2_end-t2_start)/CLOCKS_PER_SEC << endl;
	stats.close();
	return 0;
}

