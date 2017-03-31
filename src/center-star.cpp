#include <stdio.h>
#include <bitset>
#include <string.h>
#include "center-star.h"
using namespace std;


/**
* ��8 char = 16 bits ת����һ��������Ϊ�±�
* ��������ʶ�����ĸֱ�ӷ���-1
*/
int charsToIndex(const char *str) {
	bitset<16> bits(0x0000);
	for (int i = 0; i<8; i++) {
		switch (str[i]) {
		case 'A':       // 00
			break;
		case 'C':       // 01
			bits[i * 2 + 1] = 1;
			break;
		case 'T':       // 10
		case 'U':
			bits[i * 2] = 1;
			break;
		case 'G':       // 11
			bits[i * 2] = 1;
			bits[i * 2 + 1] = 1;
			break;
		default:        // ������ʶ�����ĸ������N,X��
			return -1;
		}
	}
	return (int)(bits.to_ulong());
}

/**
* һ�����е�ÿ���������ֻ������һ��
* ʹ��һ�������bool[65536]����¼�Ƿ��Ѿ��ӹ�һ��
*/
void setOccVector(const char *str, int *vec) {
	bool flag[65536] = { false };

	int len = strlen(str);
	int n = len / 8;
	for (int i = 0; i<n; i++) {
		int index = charsToIndex(str + i * 8);
		if (index >= 0 && !flag[index]) {
			vec[index]++;
			flag[index] = true;
		}
	}
}


/**
* ��ѯһ�����е�ÿһ�����������г��ֵĴ���
* ��������һ������Ϊ���Ĵ�
*/
int countSequences(const char *str, int *vec) {
	int len = strlen(str);
	int n = len / 8;
	int count = 0;
	for (int i = 0; i<n; i++) {
		int index = charsToIndex(str + i * 8);
		if (index >= 0)
			count += vec[index];
	}

	return count;
}


/**
* ��ÿ������Ϊp��С�Σ�ÿ��С�γ���Ϊ8�ֽ�
* 8 char = 16 bits��2^16���Ϊ65535
* 'A' = 00
* 'C' = 01
* 'T' = 10, 'U' = 10
* 'G' = 11
* ʹ��һ��int[65536]������ͳ��ÿһ�εĳ��ִ���
* �ٴ˱������д����ҳ������������ظ�����С�������Ĵ���Ϊ���Ĵ�
*/
int findCenterSequence(vector<string> sequences) {
	int vec[65536] = { 0 };

	for (int i = 0; i < sequences.size(); i++) {
		setOccVector(sequences[i].c_str(), vec);
	}

	int maxIndex = 0, maxCount = 0;
	for (int i = 0; i < sequences.size(); i++) {
		int count = countSequences(sequences[i].c_str(), vec);
		if (count > maxCount) {
			maxIndex = i;
			maxCount = count;
		}
		//printf("seq: %d, count: %d,\n", i++, count);
	}

	//printf("maxIndex: %d, maxCount:%d\n", maxIndex, maxCount);

	return maxIndex;
}

