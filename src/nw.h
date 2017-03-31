#ifndef _NW_H_
#define _NW_H_
#include <vector>
#include <string>

using namespace std;

/**
* centerSeq    中心串
* seqs         除中心串外的所有串
* startIdx     开始执行的串编号o
* maxLength    最长串的长度
*/
void cpu_msa(string centerSeq, vector<string> seqs, int startIdx, short *space, short *spaceForOther, int maxLength);
void backtrack(short **matrix, string centerSeq, string seq, int seqIdx, short *space, short *spaceForOther, int maxLength);
short** nw(string str1, string str2);

#endif
