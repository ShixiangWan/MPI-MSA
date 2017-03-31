#ifndef _UTIL_H_
#define _UTIL_H_
#include <string>
#include <vector>
#include <iostream>
//#include <cuda.h>

typedef struct FastaSeqs_t {
	std::vector<std::string> titles;
	std::vector<std::string> seqs;
} FastaSeqs;

/**
* ����Fasta��ʽ���ļ�
*/
FastaSeqs readFastaFile(const char *path);

/**
* ��MSA�Ľ��д�뵽Fasta��ʽ�ļ�
*/
void writeFastaFile(const char *path, std::vector<std::string> titles, std::vector<std::string> alignedSeqs);

/**
* ����main�����Ĳ���������ȫ�ֱ���
* ����argv�в���ѡ����±꣬��input��ouput��·��
*/
int parseOptions(int argc, char* argv[]);

/**
* ��ʾ������Ϣ
*/
void displayUsage();


bool configureKernel(int centerSeqLength, int maxLength, unsigned long sumLength);

#endif
