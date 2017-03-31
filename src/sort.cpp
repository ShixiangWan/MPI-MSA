#include <stdio.h>
#include <assert.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include "util.h"
using namespace std;

bool compareLength(const string &str1, const string &str2) {
	return str1.size() < str2.size();
}

/**
* ����FASTA��ʽ�ļ��Ĺ���
* ���ܣ������ļ�Ȼ���մ��ĳ�����������
*
* �÷�
* ./sort input_path output_path
*/

//int main(int argc, char* argv[]) {
//	assert(argc >= 2);
//
//	FastaSeqs fs = readFastaFile(argv[1]);
//	vector<string> seqs = fs.seqs;
//
//	vector<int> seqsSize;
//
//	printf("Sorting ...\n");
//	sort(seqs.begin(), seqs.end(), compareLength);
//	printf("Done. Max: %lu, Min: %lu\n", seqs[0].size(), seqs[seqs.size() - 1].size());
//
//	// ����ļ�
//	char *outputPath = argv[2];
//	printf("Write %lu Sequences to %s ...\n", seqs.size(), outputPath);
//	writeFastaFile(outputPath, fs.titles, seqs);
//	printf("Done.\n");
//
//	return 0;
//
//}
