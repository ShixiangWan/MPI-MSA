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
* 定义全局变量
* centerSeq 存储中心串
* seqs 存储所有其他串
*/
string centerSeq;
vector<string> titles;
vector<string> seqs;    // 所有串
int maxLength;          // 最长的串的长度
int centerSeqIdx;


void pre_compute();

/**
* 从path读如fasta格式文件，
* 完成初始化工作并输出相关信息
*/
void init(const char *path) {
	// 读入所有字符串
	// centerSeq, 图中的纵向，决定了行数m
	// seqs[idx], 图中的横向，决定了列数n
	double start = clock();
	FastaSeqs fastaSeqs = readFastaFile(path);
	titles = fastaSeqs.titles;
	seqs = fastaSeqs.seqs;
	double end = clock();
	printf("Read Sequences, use time: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

	// 找出中心串
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

	// 输出相关信息
	printf("\n\n=========================================\n");
	printf("Sequences Size: %lu\n", seqs.size() + 1);
	printf("Max: %d, Min: %d, Avg: %d\n", maxLength, minLength, avgLength);
	printf("Center Sequence Index: %d\n", centerSeqIdx);
	//printf("Workload Ratio of GPU/CPU: %.2f:%d\n", (MODE == GPU_ONLY) ? 1 : WORKLOAD_RATIO, (MODE == GPU_ONLY) ? 0 : 1);
	//printf("Block Size: %d, Thread Size: %d\n", BLOCKS, THREADS);
	printf("=========================================\n\n");
}

/**
* 将MSA结果输出到path文件中
* 共有n条串，平均长度m
* 构造带空格的中心串复杂度为:O(nm)
* 构造带空格的其他条串复杂度为:O(nm)
*/
void output(short *space, short *spaceForOther, const char* path) {
	vector<string> allAlignedStrs;

	int sWidth = centerSeq.size() + 1;      // space[] 的每条串宽度
	int soWidth = maxLength + 1;            // spaceForOther[] 的每条串宽度

	// 将所有串添加的空格汇总到一个数组中
	// 然后给中心串插入空格
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
		// 先插入自己比对时的空格
		for (int pos = seqs[idx].size(); pos >= 0; pos--) {
			if (spaceForOther[idx*soWidth + pos] > 0)
				alignedStr.insert(pos, spaceForOther[idx*soWidth + pos], '-');
		}
		// 再插入其他串比对时引入的空格
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

	// 将结果写入文件
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
		// 解析用户参数
		//int argvIdx = parseOptions(argc, argv);
		//if (argvIdx < 0) return 0;    // 输入错误选项或选项不够时不执行程序
		//const char *inputPath = argv[argvIdx];
		//const char *outputPath = argv[argvIdx + 1];
		const char *inputPath = "E:\\mpi\\MPI-demo\\x64\\Debug\\protein.txt";
		const char *outputPath = "E:\\mpi\\MPI-demo\\x64\\Debug\\results.txt";
		MODE == CPU_ONLY;

		// 读入所有串，找出中心串
		init(inputPath);

		// 记录空格的数组
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

		// MSA计算
		double cpu_time;
		double start = clock();
		int workCount = 0;
	}

	// 均分任务到每个MPI进程，每个进程计算一个for循环的量
	for (int i = my_rank; i < seqs.size(); i += num_procs){
		short **matrix = nw(centerSeq, seqs[my_rank]);
		backtrack(matrix, centerSeq, seqs[my_rank], my_rank, space, spaceForOther, maxLength);
		//printf("%d/%lu, sequence length:%lu\n", idx+1, seqs.size(), seqs[idx].size());
	}
	
	if (my_rank == 0) {
		// 输出结果
		output(space, spaceForOther, outputPath);
		delete[] space;
		delete[] spaceForOther;
		stop = MPI_Wtime();
		printf("total time: %f\n", stop - start);
	}

	MPI_Finalize();
	return 0;
}


