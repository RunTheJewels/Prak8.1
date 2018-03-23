#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>
#include <cstdlib>


const bool func_fill = true; // Заполнять ли матрицу функцией

using namespace std;

const bool timer = true;

const bool residual = true; // Считать ли невязку

double fill(int i, int j)
{
	return i>j?i:j;
}

void print_by_cols(double **A, double *b, int size, int n, int rank)
{
	MPI_Barrier(MPI_COMM_WORLD);
	for (int p = 0; p < size; p++)
	{
		if (rank == p)
		{
			if (rank == 0)
			{
				cout << "Правая часть" << endl;
				for (int i = 0; i < n; i++)
				{
					cout << b[i] << " ";
				}
				cout << endl;
			}
			cout << "Процесс " << rank << endl;
			int k = 0;
			for (int i = rank; i < n; i += size, k++)
			{
				cout << "Столбец " << i << endl;
				for (int j = 0; j < n; j++)
				{
					cout << A[k][j] << " ";
				}
				cout << endl;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

int main(int argc, char** argv)
{
	int rank, size;
	MPI_Init(&argc, &argv);	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int n;
	int numcols;
	double **A;
	double *b;
	if (func_fill)
	{
		if (argc > 1)
		{
			n = atoi(argv[1]);
		} else
		{
			if (rank == 0) cout << "Размер системы должен быть указан\n";
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			return -1;
		}
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		numcols = n / size + ((rank < n % size)?1:0);
		A = new double*[numcols];
		for (int i = 0; i < numcols; i++)
		{
			A[i] = new double[n];
		}
		b = new double[n];
		for (int i = 0; i < n; i++)
		{
			double sum_part = 0;
			int j;
			int k = 0;
			for (j = rank; j < n; j += size, k++)
			{
				A[k][i] = fill(i,j);
				if (j % 2 == 0)
				{
					sum_part += A[k][i];
				}
			}
			MPI_Reduce(&sum_part, &b[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
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
		MPI_File file;
		MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
		MPI_File_read(file, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
		numcols = n / size + ((rank < n % size)?1:0);
		//cout << "Процесс " << rank << " берёт " << numcols << " столбцов" << endl;
		A = new double*[numcols];
		for (int i = 0; i < numcols; i++)
		{
			A[i] = new double[n];
		}
		// Чтение матрицы
		MPI_File_seek(file, sizeof(int)+rank*sizeof(double), MPI_SEEK_SET);
		for (int i = 0; i < n; i++)
		{
			int j;
			int k = 0;
			for (j = rank; j < n; j += size, k++)
			{
				MPI_File_read(file, &A[k][i], 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
				MPI_File_seek(file, (size-1)*sizeof(double), MPI_SEEK_CUR);
			}
			if (i != n-1)
			{
				MPI_File_seek(file, (n-j+rank)*sizeof(double), MPI_SEEK_CUR);
			}
		}
		b = new double[n];
		MPI_File_seek(file, sizeof(int)+n*n*sizeof(double), MPI_SEEK_SET);
		MPI_File_read(file, b, n, MPI_DOUBLE, MPI_STATUS_IGNORE);
	}

	// Вывод матрицы
	// print_by_cols(A, b, size, n, rank);
	float t1_start, t1_end, t2_start, t2_end;
	t1_start = MPI_Wtime();
	// Основная часть
	for (int i = 0; i < n-1; i++)
	{
		double x[n-i];
		int root = i % size;
		if (rank == root)
		{
			double x_norm = 0;
			for (int j = 0; j < n-i; j++)
			{
				x[j] = A[i/size][j+i];
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
		}
		MPI_Bcast(x, n-i, MPI_DOUBLE, root, MPI_COMM_WORLD);
		int k = 0;
		for (int j = rank; j < n; j += size, k++)
		{
			if (j < i) continue;
			double xA = 0;
			for (int l = i; l < n; l++)
			{
				xA += x[l-i] * A[k][l];
			}
			for (int l = i; l < n; l++)
			{
				A[k][l] -= 2*xA*x[l-i];
			}
		}
		if (rank == 0)
		{
			double bx = 0;
			for (int j = i; j < n; j++)
			{
				bx += x[j-i]*b[j];
			}
			for (int j = i; j < n; j++)
			{
				b[j] -= 2*x[j-i]*bx;
			}
		}
	}
	t1_end = MPI_Wtime();
	t2_start = MPI_Wtime();
	double ans[n];
	double mypart;
	double allparts;
	MPI_Bcast(b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int i = n-1; i >= 0; i--)
	{
		mypart = 0;
		int k = 0;
		for (int j = rank; j < n; j += size, k++)
		{
			if (j < i+1) continue;
			mypart += A[k][i] * ans[j];
		}
		MPI_Reduce(&mypart, &allparts, 1, MPI_DOUBLE, MPI_SUM, i%size, MPI_COMM_WORLD);
		if (rank == i%size)
		{
			ans[i] = (b[i] - allparts) / A[i/size][i];
		}
		MPI_Bcast(&ans[i], 1, MPI_DOUBLE, i%size, MPI_COMM_WORLD);
	}	
	t2_end = MPI_Wtime();
	if (rank == 0)
	{
		for (int i = 0; i < n; i++)
		{
			cout << "x" << i << "=" << ans[i] << endl; 
		}
	}
	if (residual)
	{
		double r = 0;
		double rtmp_part, rtmp;
		for (int i = 0; i < n; i++)
		{
			rtmp_part = 0;
			int k = 0;
			for (int j = rank; j < n; j += size, k++)
			{
				rtmp_part += A[k][i] * ans[j];
			}
			MPI_Reduce(&rtmp_part, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if (rank == 0)
			{
				r += (b[i]-rtmp) * (b[i]-rtmp);
			}
		}
		if (rank == 0)
		{
			cout << "Невязка: " << sqrt(r) << endl;
		}
	}
	if (timer)
	{
		float t1 = t1_end-t1_start, t1_max;
		float t2 = t2_end-t2_start, t2_max;
		MPI_Reduce(&t1, &t1_max, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&t2, &t2_max, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			ofstream stats("stats", ofstream::app);
			stats << size << " " << n << " " << t1_max<< " " << t2_max << endl;
			stats.close();
		}
	}
	MPI_Finalize();
	return 0;
}

