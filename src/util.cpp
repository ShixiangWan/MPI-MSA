#include <fstream>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
//#include <cuda.h>
#include "unistd.h"
#include "getopt.h"
#include "util.h"
#include "global.h"
using namespace std;

FastaSeqs readFastaFile(const char *path) {
	vector<string> titles;
	vector<string> seqs;
	string buff, line, title;

	ifstream file;
	file.open(path);
	assert(file);

	while (getline(file, buff)) {
		if (buff.empty() || buff[0] == '>') {
			if (!line.empty()) {
				seqs.push_back(line);
				titles.push_back(title);
			}
			if (buff[0] == '>')
				title = buff;
			line = "";
			continue;
		}
		else {
			line += buff;
		}
	}

	if (!line.empty() && !title.empty()) {
		seqs.push_back(line);
		titles.push_back(title);
	}

	file.close();
	FastaSeqs fastaSeqs = { titles, seqs };
	return fastaSeqs;
}


void writeFastaFile(const char* path, vector<string> titles, vector<string> alignedSeqs) {
	ofstream file(path);
	if (file.is_open()) {
		for (int i = 0; i<alignedSeqs.size(); i++) {
			file << titles[i] << endl;
			int lines = alignedSeqs[i].size() / 60;        // 60���ַ�һ��
			lines = alignedSeqs[i].size() % 60 == 0 ? lines : lines + 1;
			for (int k = 0; k < lines; k++)
				file << alignedSeqs[i].substr(k * 60, 60) << endl;
		}
	}

	file.close();
}

void displayUsage() {
	printf("Usage :\n");
	printf("./msa.out [options] input_path output_path\n");
	printf("Options:\n");
	printf("\t-g\t: use GPU only (default use both GPU and CPU)\n");
	printf("\t-c\t: use CPU only (default use both GPU and CPU)\n");
	printf("\t-w <int>\t: specify the workload ratio of CPU / CPU\n");
	printf("\t-b <int>\t: specify the number of blocks per grid\n");
	printf("\t-t <int>\t: specify the number of threads per block\n");
	printf("\t-n <int>\t: specify the number of GPU devices should be used\n");
}


int parseOptions(int argc, char* argv[]) {
	if (argc < 3) {
		displayUsage();
		return -1;                          // ��ִ�г���
	}

	int oc;
	while ((oc = getopt(argc, argv, "gcw:b:t:n:")) != -1) {
		switch (oc) {
		case 'g':                       // ֻʹ��GPU
			MODE = GPU_ONLY;
			break;
		case 'c':                       // ֻʹ��CPU (OpenMP)
			MODE = CPU_ONLY;
			break;
		case 'w':                       // �����������
			WORKLOAD_RATIO = atof(optarg);
			break;
		case 'b':                       // ����Blocks����
			BLOCKS = atoi(optarg);
			break;
		case 't':                       // ����Threads����
			THREADS = atoi(optarg);
			break;
		case 'n':
			GPU_NUM = atoi(optarg);     // ����ʹ�õ�GPU����
			break;
		case '?':                       // �������ѡ���ִ�г���
			displayUsage();
			return -1;
		}
	}

	return optind;
}


//bool configureKernel(int centerSeqLength, int maxLength, unsigned long sumLength) {
//
//	size_t freeMem, totalMem;
//	cudaMemGetInfo(&freeMem, &totalMem);
//
//	// ÿ������ƥ���DP��������Ҫ�Ŀռ�
//	size_t matrixSize = sizeof(short)* (centerSeqLength + 1) * (maxLength + 1);
//
//	// �õ�ÿ��Kernel����ִ�еĴ��������ɲ��������߳���BLOCKS*THREADS��
//	// ��Ӧ��ʹ�����еĿ����ڴ棬�ڴ�����һ����20%
//	freeMem = (freeMem - sizeof(char)* sumLength) / 10 * 8;        // *0.8�м��ת��double�п��ܽض�
//	int seqs = freeMem / matrixSize;
//
//	printf("freeMem: %luMB, sumLengthSize: %luMB, matrix :%luKB, seqs: %d\n", freeMem / 1024 / 1024, sumLength / 1024 / 1024, matrixSize / 1024, seqs);
//
//	// ���ж��û����õ�<BLOCKS, THREADS>�Ƿ������ڴ�����
//	// ��������㣬���Զ�����һ��<BLOCKS, THREADS>
//	if (seqs >= BLOCKS*THREADS)
//		return true;
//
//	// �������ڴ����Ƶ�ǰ���£�
//	// ����BLOCKS >= 3 �� THREADS >= 32�������GPUִ��
//	int b, t;
//	for (t = THREADS; t >= 32; t -= 32) {
//		b = seqs / t;
//		if (b >= 3 && t >= 32) {
//			BLOCKS = b;
//			THREADS = t;
//			return true;
//		}
//	}
//
//	return false;
//}
