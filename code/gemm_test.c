// Copyright 2017 Mobvoi Inc. All Rights Reserved.
// Author: congfu@mobvoi.com (Cong Fu)

/*
 * gcc -O3 gemm_test.c -o gemm_test
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#define DIM 64

typedef int int32;

long long currentTimeInMicroseconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return ((tv.tv_sec * 1000 * 1000) + tv.tv_usec);
}

// Main: perform gemm of (A, B, C) (that is, C = alpha * A * B + beta * C)
int main(int argc, const char *argv[]) {

  long long slow_cost = 0;
  float alpha = 1.1f, beta = 0.5f;
  printf("the dim is %d\n", DIM);
  float A[DIM][DIM] = {0};
  float B[DIM][DIM] = {0};
  float C[DIM][DIM] = {0};

  int32 iter = 10;

  for (int32 n = 0; n < iter; ++n) {
    printf("C operation: iteration %d\n", n);
    for (int i = 0; i < DIM; ++i) {
      for (int j = 0; j < DIM; ++j) {
        A[i][j] = -0.5f + (float)rand() / RAND_MAX;
        B[i][j] = -0.5f + (float)rand() / RAND_MAX;
        C[i][j] = -0.5f + (float)rand() / RAND_MAX;
      }
    }

    long long c_start = currentTimeInMicroseconds();
    for (int i = 0; i < DIM; ++i) {
      for (int j = 0; j < DIM; ++j) {
        float temp = 0;
        for (int k = 0; k < DIM; ++k) {
          temp += A[i][k] * B[k][j];
        }
        C[i][j] *= beta;
        C[i][j] += alpha * temp;
      }
    }
    long long c_end = currentTimeInMicroseconds();
    slow_cost += c_end - c_start;
  }

  printf("C operation takes: %f ms\n", slow_cost / 1000.0f);

  return 0;
}
