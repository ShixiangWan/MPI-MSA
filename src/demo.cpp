#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

int main2(){

	cout << "OK";

	int my_rank, num_procs;
	int argc;
	char **argv;
	
	MPI_Init(&argc, &argv);
	MPI_Status status;
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	int test_num = 0;
	//MPI_Barrier(MPI_COMM_WORLD);     // 同步
	if (my_rank == 0){
		test_num = 1024;
		MPI_Send(&test_num, 1, MPI_INT, 1, 99, MPI_COMM_WORLD);   // MPI_Type_vector
		vector<string> test_string;
		test_string.push_back("ABCD");
		test_string.push_back("EFGH");
		//MPI_Send(&test_string, 1, MPI_Type_vector, 1, 99, MPI_COMM_WORLD);   // MPI_Type_vector
		printf("send: %d", &test_num);
	} else if (my_rank == 1) {
		MPI_Recv(&test_num, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
		printf("received: %d", test_num);
	}
	MPI_Finalize();
	return 0;
	// 从进程号为0的进程，广播1个MPI_INT类型的&n到当前进程
	//MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//
	//width = 1.0 / n;
	//sum = 0.0;
	//for (i = my_rank; i < n; i += num_procs){
	//	local = width * ((double)i - 0.5);
	//	sum += 4.0 / (1.0 + local * local);
	//}
	//mypi = width * sum;
	//MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//
	//if (my_rank == 0){
	//	printf("PI is %.20f\n", pi);
	//	stop = MPI_Wtime();
	//	printf("Time: %f on %s\n", stop - start, processor_name);
	//	fflush(stdout);
	//}
}
