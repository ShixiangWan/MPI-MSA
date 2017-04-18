#include <stdio.h>
#include <time.h>
#include "util.h"
#include "sp.h"
#include "center-star.h"
//#include "cuda-nw.h"
#include "nw.h"
//#include "omp.h"
#include "global.h"
#include <mpi.h>

using namespace std;


/**
* ����ȫ�ֱ���
* centerSeq �洢���Ĵ�
* seqs �洢����������
*/
string centerSeq;
vector<string> titles;
vector<string> seqs;    // ���д�
int maxLength;          // ��Ĵ��ĳ���
int centerSeqIdx;


void pre_compute();

/**
* ��path����fasta��ʽ�ļ���
* ��ɳ�ʼ����������������Ϣ
*/
void init(const char *path) {
	// ���������ַ���
	// centerSeq, ͼ�е����򣬾���������m
	// seqs[idx], ͼ�еĺ��򣬾���������n
	double start = clock();
	FastaSeqs fastaSeqs = readFastaFile(path);
	titles = fastaSeqs.titles;
	seqs = fastaSeqs.seqs;
	double end = clock();
	printf("Read Sequences, use time: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

	// �ҳ����Ĵ�
	start = clock();
	centerSeqIdx = findCenterSequence(seqs);
	end = clock();
	printf("Find the Center Sequence, use time: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

	centerSeq = seqs[centerSeqIdx];
	seqs.erase(seqs.begin() + centerSeqIdx);

	unsigned long sumLength = 0;
	maxLength = centerSeq.size();
	int minLength = centerSeq.size();
	for (int i = 0; i<seqs.size(); i++) {
		sumLength += seqs[i].size();
		if (maxLength < seqs[i].size())
			maxLength = seqs[i].size();
		if (minLength > seqs[i].size())
			minLength = seqs[i].size();
	}
	int avgLength = sumLength / seqs.size();

	MODE = CPU_ONLY;

	// ��������Ϣ
	printf("\n\n=========================================\n");
	printf("Sequences Size: %lu\n", seqs.size() + 1);
	printf("Max: %d, Min: %d, Avg: %d\n", maxLength, minLength, avgLength);
	printf("Center Sequence Index: %d\n", centerSeqIdx);
	//printf("Workload Ratio of GPU/CPU: %.2f:%d\n", (MODE == GPU_ONLY) ? 1 : WORKLOAD_RATIO, (MODE == GPU_ONLY) ? 0 : 1);
	//printf("Block Size: %d, Thread Size: %d\n", BLOCKS, THREADS);
	printf("=========================================\n\n");
}

/**
* ��MSA��������path�ļ���
* ����n������ƽ������m
* ������ո�����Ĵ����Ӷ�Ϊ:O(nm)
* ������ո�������������Ӷ�Ϊ:O(nm)
*/
void output(short *space, short *spaceForOther, const char* path) {
	vector<string> allAlignedStrs;

	int sWidth = centerSeq.size() + 1;      // space[] ��ÿ�������
	int soWidth = maxLength + 1;            // spaceForOther[] ��ÿ�������

	// �����д���ӵĿո���ܵ�һ��������
	// Ȼ������Ĵ�����ո�
	string alignedCenter(centerSeq);
	vector<int> spaceForCenter(centerSeq.size() + 1, 0);
	for (int pos = centerSeq.size(); pos >= 0; pos--) {
		int count = 0;
		for (int idx = 0; idx < seqs.size(); idx++)
			count = (space[idx*sWidth + pos] > count) ? space[idx*sWidth + pos] : count;
		spaceForCenter[pos] = count;
		if (spaceForCenter[pos] > 0)
			//printf("pos:%d, space:%d\n", pos, spaceForCenter[pos]);
			alignedCenter.insert(pos, spaceForCenter[pos], '-');
	}

	//printf("\n\n%s\n", alignedCenter.c_str());
	//allAlignedStrs.push_back(alignedCenter);

	for (int idx = 0; idx < seqs.size(); idx++) {
		int shift = 0;
		string alignedStr(seqs[idx]);
		// �Ȳ����Լ��ȶ�ʱ�Ŀո�
		for (int pos = seqs[idx].size(); pos >= 0; pos--) {
			if (spaceForOther[idx*soWidth + pos] > 0)
				alignedStr.insert(pos, spaceForOther[idx*soWidth + pos], '-');
		}
		// �ٲ����������ȶ�ʱ����Ŀո�
		for (int pos = 0; pos < spaceForCenter.size(); pos++) {
			int num = spaceForCenter[pos] - space[idx*sWidth + pos];
			if (num > 0) {
				alignedStr.insert(pos + shift, num, '-');
			}
			shift += spaceForCenter[pos];
		}
		//printf("%s\n", alignedStr.c_str());
		allAlignedStrs.push_back(alignedStr);
	}
	allAlignedStrs.insert(allAlignedStrs.begin() + centerSeqIdx, alignedCenter);

	// �����д���ļ�
	writeFastaFile(path, titles, allAlignedStrs);

}


int main(int argc, char *argv[]) {
	int my_rank, num_procs;
	int argc;
	char **argv;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	double start, stop;
	if (my_rank == 0) {
		start = MPI_Wtime();
		// �����û�����
		//int argvIdx = parseOptions(argc, argv);
		//if (argvIdx < 0) return 0;    // �������ѡ���ѡ���ʱ��ִ�г���
		//const char *inputPath = argv[argvIdx];
		//const char *outputPath = argv[argvIdx + 1];
		const char *inputPath = "E:\\mpi\\MPI-demo\\x64\\Debug\\protein.txt";
		const char *outputPath = "E:\\mpi\\MPI-demo\\x64\\Debug\\results.txt";
		MODE == CPU_ONLY;

		// �������д����ҳ����Ĵ�
		init(inputPath);

		// ��¼�ո������
		short *space = new short[seqs.size() * (centerSeq.size() + 1)];
		short *spaceForOther = new short[seqs.size() * (maxLength + 1)];
		int space_len = seqs.size() * (centerSeq.size() + 1);
		int spaceForOther_len = seqs.size() * (maxLength + 1);
		for (int i = 0; i < space_len; i++) {
			space[i] = 0;
		}
		for (int i = 0; i < spaceForOther_len; i++) {
			spaceForOther[i] = 0;
		}

		// MSA����
		double cpu_time;
		double start = clock();
		int workCount = 0;
	}

	// ��������ÿ��MPI���̣�ÿ�����̼���һ��forѭ������
	for (int i = my_rank; i < seqs.size(); i += num_procs){
		short **matrix = nw(centerSeq, seqs[my_rank]);
		backtrack(matrix, centerSeq, seqs[my_rank], my_rank, space, spaceForOther, maxLength);
		//printf("%d/%lu, sequence length:%lu\n", idx+1, seqs.size(), seqs[idx].size());
	}
	
	if (my_rank == 0) {
		// ������
		output(space, spaceForOther, outputPath);
		delete[] space;
		delete[] spaceForOther;
		stop = MPI_Wtime();
		printf("total time: %f\n", stop - start);
	}

	MPI_Finalize();
	return 0;
}


