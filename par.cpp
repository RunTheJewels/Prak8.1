#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>


using namespace std;

const bool timer = true;

void print_by_cols(float **A, float *b, int size, int n, int rank)
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
	int n;
	MPI_File_read(file, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
	int numcols = n / size + ((rank < n % size)?1:0);
	//cout << "Процесс " << rank << " берёт " << numcols << " столбцов" << endl;
	float **A = new float*[numcols];
	for (int i = 0; i < numcols; i++)
	{
		A[i] = new float[n];
	}
	float tmpA[n];
	// Чтение матрицы
	MPI_File_seek(file, sizeof(int)+rank*sizeof(float), MPI_SEEK_SET);
	for (int i = 0; i < n; i++)
	{
		int j;
		int k = 0;
		for (j = rank; j < n; j += size, k++)
		{
			MPI_File_read(file, &A[k][i], 1, MPI_FLOAT, MPI_STATUS_IGNORE);
			MPI_File_seek(file, (size-1)*sizeof(float), MPI_SEEK_CUR);
		}
		if (i != n-1)
		{
			MPI_File_seek(file, (n-j+rank)*sizeof(float), MPI_SEEK_CUR);
		}
	}
	float *b;
	float *tmpb;
	b = new float[n];
	tmpb = new float[n];
	MPI_File_seek(file, sizeof(int)+n*n*sizeof(float), MPI_SEEK_SET);
	MPI_File_read(file, b, n, MPI_FLOAT, MPI_STATUS_IGNORE);


	// Вывод матрицы
	// print_by_cols(A, b, size, n, rank);
	double t1_start, t1_end, t2_start, t2_end;
	t1_start = MPI_Wtime();
	// Основная часть
	for (int i = 0; i < n-1; i++)
	{
		float x[n-i];
		int root = i % size;
		if (rank == root)
		{
			float x_norm = 0;
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
		MPI_Bcast(x, n-i, MPI_FLOAT, root, MPI_COMM_WORLD);
		int k = 0;
		for (int j = rank; j < n; j += size, k++)
		{
			if (j < i) continue;
			float xA = 0;
			for (int l = i; l < n; l++)
			{
				xA += x[l-i] * A[k][l];
			}
			for (int l = i; l < n; l++)
			{
				A[k][l] -= 2*xA*x[l-i];
			}
			/*for (int l = i; l < n; l++)
			{
				tmpA[l] = A[k][l];
				for (int m = i; m < n; m++)
				{
					tmpA[l] -= 2*x[l-i]*x[m-i]*A[k][m];
				}
			}
			for (int l = i; l < n; l++)
			{
				A[k][l] = tmpA[l];
			}*/
		}
		if (rank == 0)
		{
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
				b[j] = tmpb[j];
			}*/
		}
	}
	t1_end = MPI_Wtime();
	t2_start = MPI_Wtime();
	float ans[n];
	float mypart;
	float allparts;
	MPI_Bcast(b, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

	for (int i = n-1; i >= 0; i--)
	{
		mypart = 0;
		int k = 0;
		for (int j = rank; j < n; j += size, k++)
		{
			if (j < i+1) continue;
			mypart += A[k][i] * ans[j];
		}
		MPI_Reduce(&mypart, &allparts, 1, MPI_FLOAT, MPI_SUM, i%size, MPI_COMM_WORLD);
		if (rank == i%size)
		{
			ans[i] = (b[i] - allparts) / A[i/size][i];
		}
		MPI_Bcast(&ans[i], 1, MPI_FLOAT, i%size, MPI_COMM_WORLD);
	}	
	t2_end = MPI_Wtime();
	if (rank == 0)
	{
		for (int i = 0; i < n; i++)
		{
			cout << "x" << i << "=" << ans[i] << endl; 
		}
	}
	if (rank == 0 and timer)
	{
		ofstream stats("stats", ofstream::app);
		stats << size << " " << n << " " << t1_end-t1_start << " " << t2_end-t2_start << endl;
		stats.close();
	}
	// Вывод матрицы
/*	MPI_Barrier(MPI_COMM_WORLD);
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
	}*/
	MPI_Finalize();
	return 0;
}

