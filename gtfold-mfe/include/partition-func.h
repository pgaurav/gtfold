#ifndef _PARTITION_DANGLE_H
#define _PARTITION_DANGLE_H


#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

extern	double ** u;
extern	double ** up;
extern	double ** upm;
extern	double ** ud;
extern	double ** u1d;

extern	double ** s1;
extern	double ** s2;
extern	double ** s3;
extern	double ** u1;

extern int part_len;
double ED3_new(int i, int j, int k);
double ED5_new(int i, int j, int k);
double f(int j, int h, int l);
void calculate_partition(int len);
void free_partition();

#ifdef __cplusplus
}
#endif

#endif
