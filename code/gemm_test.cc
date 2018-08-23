// Copyright 2017 Mobvoi Inc. All Rights Reserved.
// Author: congfu@mobvoi.com (Cong Fu)

/*
 * g++ -O3 gemm_test.cc -o gemm_test -std=c++11
 *
 * mt2601 189818ms
 * arm8 android phone 24911ms
 * x86_64 linux pc 4613ms
 */

#include <cstdlib>
#include <sys/time.h>
#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

typedef int32_t int32;

typedef enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, CblasConjNoTrans=114} CBLAS_TRANSPOSE;

template <typename Real>
Real *CreateMatrix(const int32 row, const int32 col, const int32 stride,
                   const int32 seed = 777);

template <typename Real> void DestroyMatrix(Real *data);

template <typename Real>
void ShowMatrix(const std::string id, const Real *data, const int32 row,
                const int32 col, const int32 stride);

template <typename Real>
Real EuclideanDistance(const Real *X, const Real *Y, const int32 row,
                       const int32 col, const int32 stride);

template <typename Real>
Real slow_dot(const int32 dim, const Real *x, const int32 x_inc, const Real *y,
              const int32 y_inc);

template <typename Real>
void slow_gemm(const Real alpha, const CBLAS_TRANSPOSE A_trans, const Real *A,
               const int32 A_rows, const int32 A_cols, const int32 A_stride,
               const CBLAS_TRANSPOSE B_trans, const Real *B,
               const int32 B_stride, const Real beta, Real *C,
               const int32 C_rows, const int32 C_cols, const int32 C_stride);

long long currentTimeInMilliseconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return ((tv.tv_sec * 1000) + (tv.tv_usec / 1000));
}

// Main: perform gemm of (A, B, C) (that is, C = alpha * A * B + beta * C)
// A, B, C are matrices here
// alpha, beta are constants here
// We use two versions of C, X and Y, to test the speed and accuracy of
// C and OpenBlas operations.

int main(int argc, const char *argv[]) {
  int32 A_rows = 512, A_cols = 512, A_stride = 512, B_rows = 512, B_cols = 512,
        B_stride = 512, C_rows = 512, C_cols = 512, C_stride = 512;
  //duration<double> time_span;
  float alpha = 0.5, beta = 1.5;

  int32 iter = 10;

  long long slow_cost = 0;
  for (int32 n = 0; n < iter; ++n) {
    std::cout << "C operation: iteration " << n << std::endl;
    float *A = CreateMatrix<float>(A_rows, A_cols, A_stride, 'A');
    float *B = CreateMatrix<float>(B_rows, B_cols, B_stride, 'B');
    float *X = CreateMatrix<float>(C_rows, C_cols, C_stride, 'C');

    long long c_start = currentTimeInMilliseconds();
    slow_gemm<float>(alpha, CblasNoTrans, A, A_rows, A_cols, A_stride,
                     CblasNoTrans, B, B_stride, beta, X, C_rows, C_cols,
                     C_stride);
    long long c_end = currentTimeInMilliseconds();

    slow_cost += c_end - c_start;
    DestroyMatrix<float>(A);
    DestroyMatrix<float>(B);
    DestroyMatrix<float>(X);
  }
  std::cout << "C operation takes: " << slow_cost << " ms" << std::endl;
  
  return 0;
}

// Implementations of functions.

template <typename Real>
Real *CreateMatrix(const int32 row, const int32 col, const int32 stride,
                   const int32 seed) {
  Real *data = static_cast<Real *>(malloc(row * stride * sizeof(Real)));
  srand(seed);
  for (int32 r = 0; r < row; ++r) {
    for (int32 c = 0; c < col; ++c) {
      data[r * stride + c] = -1.0 + static_cast<Real>(rand()) /
                                        (static_cast<Real>(RAND_MAX / 2.0));
    }
  }
  return data;
}

template <typename Real> void DestroyMatrix(Real *data) { free(data); }

template <typename Real>
void ShowMatrix(const std::string id, const Real *data, const int32 row,
                const int32 col, const int32 stride) {
  std::cout << id << " [";
  for (int32 r = 0; r < row; ++r) {
    std::cout << "\n  ";
    for (int32 c = 0; c < col; ++c) {
      std::cout << data[r * stride + c] << " ";
    }
  }
  std::cout << "]\n";
}

template <typename Real>
Real EuclideanDistance(const Real *X, const Real *Y, const int32 row,
                       const int32 col, const int32 stride) {
  Real res = static_cast<Real>(0);
  for (int32 r = 0; r < row; ++r, X += stride, Y += stride) {
    for (int32 c = 0; c < col; ++c) {
      Real diff = X[c] - Y[c];
      res += (diff * diff);
    }
  }
  return res;
}

template <typename Real>
Real slow_dot(const int32 dim, const Real *x, const int32 x_inc, const Real *y,
              const int32 y_inc) {
  Real sum = static_cast<Real>(0);
  for (int32 d = 0; d < dim; ++d, x += x_inc, y += y_inc)
    sum += (*y) * (*x);
  return sum;
}

template <typename Real>
void slow_gemm(const Real alpha, const CBLAS_TRANSPOSE A_trans, const Real *A,
               const int32 A_rows, const int32 A_cols, const int32 A_stride,
               const CBLAS_TRANSPOSE B_trans, const Real *B,
               const int32 B_stride, const Real beta, Real *C,
               const int32 C_rows, const int32 C_cols, const int32 C_stride) {
  int32 A_row_inc = A_stride;
  int32 A_col_inc = 1;
  if (A_trans == CblasTrans) {
    std::swap(A_row_inc, A_col_inc);
  }
  int32 B_row_inc = B_stride;
  int32 B_col_inc = 1;
  if (B_trans == CblasTrans) {
    std::swap(B_row_inc, B_col_inc);
  }
  const Real *A_row_data = A;
  Real *C_row_data = C;
  int32 dim = (A_trans == CblasNoTrans) ? A_cols : A_rows;
  for (int32 r = 0; r < C_rows;
       ++r, C_row_data += C_stride, A_row_data += A_row_inc) {
    const Real *B_col_data = B;
    for (int32 c = 0; c < C_cols; ++c, B_col_data += B_col_inc) {
      C_row_data[c] *= beta;
      C_row_data[c] +=
          alpha * slow_dot(dim, A_row_data, A_col_inc, B_col_data, B_row_inc);
    }
  }
}
