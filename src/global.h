#ifndef _GLOBAL_H_
#define _GLOBAL_H_

/**
* ����ȫ�ֱ���
*/

#define MISMATCH -1
#define MATCH 0
#define GAP -1
#define GAP_START 0
#define GAP_EXTEND -1
#define GAP_SE -1           // GAP_START + GAP_EXTEND

#define MAX_THREADS 1024
#define MAX_BLOCKS 128

#define MIN_SCORE -32700


// ÿ��Kernel�е�Block����
extern int BLOCKS;

// ÿ��Block��Thread����
extern int THREADS;

// �������ı���, GPU/CPU
extern double WORKLOAD_RATIO;

// ʹ��GPU������, Ĭ��Ϊ0
// �ɳ����ȡGPU����������ʼֵ��Ϊ0˵���û��Լ�ָ��
extern int GPU_NUM;

// ���з�ʽ
#define GPU_ONLY 1      // ֻʹ��GPU
#define CPU_ONLY 2      // ֻʹ��CPU
#define CPU_GPU 3       // ͬʱʹ��GPU��CPU
extern int MODE;        // 1, 2, 3

#endif
