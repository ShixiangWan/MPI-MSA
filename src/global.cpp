#include "global.h"

int BLOCKS = 12;
int THREADS = 256;

double WORKLOAD_RATIO = 1;  // Ĭ��GPU/CPU�������1:1�������Ը���һ��Ĵ�

int MODE = CPU_GPU;         // Ĭ��ͬʱʹ��GPU��CPU

int GPU_NUM = 0;            // �ɳ����Զ���ȡGPU�������û�Ҳ����ͨ������ָ��
