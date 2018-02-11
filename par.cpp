#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>


using namespace std;

int main(int argc, char** argv)
{
	int rank, size;
	MPI_Init(&argc, &argv);	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_File file;
	MPI_File_open(MPI_COMM_WORLD, "matrixn.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	int n;
	MPI_File_read_all(file, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_close(&file);
	cout << "Процесс " << rank << " прочитал значение " << n << endl;
	int start = n * rank / size;
	int end = n * (rank+1) / size -1;
	cout << "Процесс " << rank << " берёт столбцы с " << start << " до " << end << endl;
	MPI_Finalize();
	return 0;
}

