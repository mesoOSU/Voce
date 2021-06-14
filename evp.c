#define EVP_GLOBALS_ONCE
#include<assert.h>
#include<mpi.h>
#include "evp.h"
#undef EVP_GLOBALS_ONCE

int main(int argc, char **argv)
{
	// inputfile is needed as 2nd argument on the command line
	assert(argc>1);
	GetNameList(argv);

	/* MPI initialization */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NumPE);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank==0) PrintNameList();

	if(rank==0) OpenFiles();
	SetupJob();

	Evolution();

	DestroyJob();
	if(rank==0) CloseFiles();

	/* MPI Finalize */
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}

