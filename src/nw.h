#ifndef _NW_H_
#define _NW_H_
#include <vector>
#include <string>

using namespace std;

/**
* centerSeq    ���Ĵ�
* seqs         �����Ĵ�������д�
* startIdx     ��ʼִ�еĴ����o
* maxLength    ����ĳ���
*/
void cpu_msa(string centerSeq, vector<string> seqs, int startIdx, short *space, short *spaceForOther, int maxLength);
void backtrack(short **matrix, string centerSeq, string seq, int seqIdx, short *space, short *spaceForOther, int maxLength);
short** nw(string str1, string str2);

#endif
