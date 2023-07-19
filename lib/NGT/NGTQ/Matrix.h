//
// Copyright (C) 2021 Yahoo Japan Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <cstring>
#include <vector>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <sstream>
#include <random>

extern "C" {
  // svd
  void dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a,
	       int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
	       double* work, int* lwork, int* info);
  void sgesvd_(char* jobu, char* jobvt, int* m, int* n, float* a,
	       int* lda, float* s, float* u, int* ldu, float* vt, int* ldvt,
	       float* work, int* lwork, int* info);
  // multiplication
  void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
	      double *alpha, double *a, int *lda, double *b, int *ldb,
	      double *beta , double *c, int *ldc);
  void sgemm_(char *transa, char *transb, int *m, int *n, int *k,
	      float *alpha, float *a, int *lda, float *b, int *ldb,
	      float *beta , float *c, int *ldc);
  // {D|S}GEQRF computes a QR factorization of a M-by-N matrix A: A = Q * R.
  // {D|S}ORGQR return Q  =  H(1) H(2) . . . H(k) as returned by {D|S}GEQRF.
  void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
  void dorgqr_(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

  void sgeqrf_(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
  void sorgqr_(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
}



template <typename T> class Matrix {
public:
  Matrix(size_t r = 0, size_t c = 0, const float *v = 0): row(r), col(c), matrix(0) { construct(r, c, v); }
  Matrix(const Matrix<T> &m): row(0), col(0), matrix(0) { *this = m; }
  Matrix(const std::vector<std::vector<float>> &v): row(0), col(0), matrix(0) { construct(v); }

  ~Matrix() { delete[] matrix; }

  Matrix<T> &operator=(const Matrix<T> &m) {
    allocate(m.row, m.col);
    std::memcpy(matrix, m.matrix, row * col * sizeof(T));
    return *this;
  }

  void construct(size_t r, size_t c, const double *v) {
    allocate(r, c);
    set(v);
  }

  void construct(size_t r, size_t c, const float *v) {
    allocate(r, c);
    set(v);
  }

  void construct(const std::vector<std::vector<float>> &v) {
#if !defined(NGT_DISABLE_BLAS)
    allocate(v[0].size(), v.size());
#else
    allocate(v.size(), v[0].size());
#endif
    set(v);
  }

  void allocate(size_t r, size_t c) {
    if (matrix != 0) {
      delete[] matrix;
    }
    row = r;
    col = c;
    if (r == 0 && c == 0) {
      matrix = 0;
    } else {
      matrix = new T[r * c];
    }
  }

  bool isEmpty() { return (col == 0) && (row == 0); }

  static void
    tokenize(const std::string &str, std::vector<std::string> &token, const std::string seps) {
    std::string::size_type current = 0;
    std::string::size_type next;
    while ((next = str.find_first_of(seps, current)) != std::string::npos) {
      token.push_back(str.substr(current, next - current));
      current = next + 1;
    }
    std::string t = str.substr(current);
    token.push_back(t);
  }

  void set(const double *v) {
    if (v == 0) {
      return;
    }
    size_t l = row * col;
    for (size_t p = 0; p < l; p++) {
      matrix[p] = *v++;
    }
  }

  void set(const float *v) {
    if (v == 0) {
      return;
    }
    size_t l = row * col;
    for (size_t p = 0; p < l; p++) {
      matrix[p] = *v++;
    }
  }

  void set(const std::vector<std::vector<float>> &v) {
    T *m = matrix;
#if !defined(NGT_DISABLE_BLAS)
    assert(row == v[0].size());
    assert(col == v.size());
    for (size_t c = 0; c < col; c++) {
      for (size_t r = 0; r < row; r++) {
	*m++ = v[c][r];
      }
    }
#else
    assert(row == v.size());
    assert(col == v[0].size());
    for (size_t r = 0; r < row; r++) {
      for (size_t c = 0; c < col; c++) {
	*m++ = v[r][c];
      }
    }
#endif
  }

  void set(size_t pr, size_t pc, T v) {
#if !defined(NGT_DISABLE_BLAS)
    matrix[pc * row + pr] = v;
#else
    matrix[pr * col + pc] = v;
#endif
  }

  void put(size_t pr, size_t pc, const Matrix<T> &m) {
    for (size_t r = 0; r < m.row; r++) {
      if (pr + r < row) {
	for (size_t c = 0; c < m.col; c++) {
	  if (pc + c < col) {
	    matrix[(pr + r) * col + (pc + c)] = m.matrix[r * m.col + c];
	  }
	}
      }
    }
  }

  void horzcat(const Matrix<T> &m){
    assert(row == m.row);
    size_t nc = col + m.col;
    T *mtx = new T[row * nc];
    for (size_t r = 0; r < row; r++) {
      for (size_t c = 0; c < col; c++) {
	mtx[r * nc + c] = matrix[r * col + c];
      }
    }
    put(0, col, m);
    col = nc;
    delete[] matrix;
    matrix = mtx;
  }

  void vert(const Matrix<T> &m) {
    if (row == 0 && col == 0) {
      construct(m.row, m.col, m.matrix);
      return;
    }
    assert(col == m.col);
    size_t nr = row + m.row;
    T *mtx = new T[nr * col];
    for (size_t r = 0; r < row; r++) {
      for (size_t c = 0; c < col; c++) {
	mtx[r * col + c] = matrix[r * col + c];
      }
    }
    put(row, 0, m);
    row = nr;
    delete[] matrix;
    matrix = mtx;
  }

  void zero(size_t r, size_t c = 0) {
    allocate(r, c);
    zero();
  }

  void zero() {
    size_t l = row * col;
    for (size_t p = 0; p < l; p++) {
      matrix[p] = 0.0;
    }
  }

  void random(size_t r, size_t c = 0) {
    allocate(r, c);
    random();
  }

  void random() {
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());

    auto maxNum = engine.max();
    size_t l = row * col;
    for (size_t p = 0; p < l; p++) {
      auto randomNum = engine();
      matrix[p] = static_cast<T>(randomNum) / maxNum;
    }
  }

  static void random(std::vector<T> &a) {
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    auto maxNum = engine.max();
    for (auto &i : a) {
      auto randomNum = engine();
      i = static_cast<T>(randomNum) / maxNum;
    }
  }

  void transpose() {
    T *m = new T[col * row];
    T *msrc = matrix;
    size_t nr = col;
    size_t nc = row;
#if !defined(NGT_DISABLE_BLAS)
    for (size_t r = 0; r < nr; r++) {
      for (size_t c = 0; c < nc; c++) {
	m[c * nr + r] = *msrc++;
	//std::cerr << r * nc + c << std::endl;
      }
    }
#else
    for (size_t c = 0; c < nc; c++) {
      for (size_t r = 0; r < nr; r++) {
	m[r * nc + c] = *msrc++;
	//std::cerr << r * nc + c << std::endl;
      }
    }
#endif
    row = nr;
    col = nc;
    delete[] matrix;
    matrix = m;
  }

#if !defined(NGT_DISABLE_BLAS)
  void mul(const std::vector<std::vector<float>> &v) {
    Matrix<T> m(v);
    mulBlas(m, true);
  }


  void mul(const Matrix<T> &mtx) {
    mulBlas(mtx);
  }


  void mulBlas(const Matrix<T> &mtx, bool transpose = false) {
    char transa = 'N';
    char transb = 'N';
    int m = row;
    int n = mtx.col;
    int k = col;
    if (transpose) {
      transb = 'T';
      if (row != mtx.row) {
	std::cerr << "mul:" << row << "x" << mtx.row << std::endl;
      }
      assert(row == mtx.row);
      n = mtx.row;
      row = m;
      col = mtx.row;
    } else {
      if (col != mtx.row) {
	std::cerr << "mul:" << col << "x" << mtx.row << std::endl;
      }
      assert(col == mtx.row);
      row = m;
      col = n;
    }
    float alpha = 1.0;
    float beta = 0.0;
    T *tmpmtx = new T[m * n];
    if (transpose) {
      int ldb = mtx.row;
      gemm(&transa, &transb, &m, &n, &k, &alpha, matrix, &m, mtx.matrix, &ldb, &beta, tmpmtx, &m);
    } else {
      gemm(&transa, &transb, &m, &n, &k, &alpha, matrix, &m, mtx.matrix, &k, &beta, tmpmtx, &m);
    }
    delete[] matrix;
    matrix = tmpmtx;
  }
#else
  void mul(const Matrix<T> &mtx) {
    mulNaive(mtx);
  }
#endif


  void mulNaive(const Matrix<T> &mtx) {
#ifdef MATRIX_TRACE
    cerr << row << "x" << col << " mtx=" << mtx.row << "x" << mtx.col << std::endl;
    std::cerr << mtx << std::endl;
#endif
    if (col != mtx.row) {
      std::cerr << "mul:" << col << "x" << mtx.row << std::endl;
    }
    assert(col == mtx.row);
    size_t nr = row;
    size_t nc = mtx.col;
    T *tmpmtx = new T[nr * nc];
    for (size_t r = 0; r < nr; r++) {
      for (size_t c = 0; c < nc; c++) {
#if !defined(NGT_DISABLE_BLAS)
	T &sum = tmpmtx[c * nr + r];
	sum = 0;
	for (size_t p = 0; p < col; p++) {
	  sum += matrix[p * row + r] * mtx.matrix[c * mtx.row + p];
	}
#else
	T &sum = tmpmtx[r * nc + c];
	sum = 0;
	for (size_t p = 0; p < col; p++) {
	  sum += matrix[r * col + p] * mtx.matrix[p * mtx.col + c];
	}
#endif
      }
    }
    row = nr;
    col = nc;
    delete[] matrix;
    matrix = tmpmtx;
  }

  void diag(const Matrix<T> &m) {
    if (m.row != 1 && m.col != 1) {
      std::cerr << "Error : not vector. " << m.row << "x" << m.col << std::endl;
      return;
    }
    size_t length = m.row > m.col ? m.row : m.col;
    zero(length, length);
    for (size_t i = 0; i < length; i++) {
#if !defined(NGT_DISABLE_BLAS)
      matrix[i * row + i] = m.matrix[i];
#else
      matrix[i * col + i] = m.matrix[i];
#endif
    }
  }

  void reshape(size_t r, size_t c) {
    if (r == row && c == col) {
      return;
    }
    size_t l = r * c;
    T *m = new T[l];
    for (size_t i = 0; i < l; i++) {
      m[i] = 0.0;
    }
#if !defined(NGT_DISABLE_BLAS)
    for (size_t sc = 0; sc < col; sc++) {
      for (size_t sr = 0; sr < row; sr++) {
	m[sc * r + sr] = matrix[sc * row + sr];
      }
    }
#else
    for (size_t sr = 0; sr < row; sr++) {
      for (size_t sc = 0; sc < col; sc++) {
        m[sr * c + sc] = matrix[sr * col + sc];
      }
    }
#endif
    row = r;
    col = c;
    delete[] matrix;
    matrix = m;
  }

  static void mulSquare(std::vector<float> &a, Matrix<T> &b) {
    if (b.col != b.row) {
      std::stringstream msg;
      msg << "mulSquare : Invalid # of cols and rows. " << b.col << ":" << b.row << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    if (a.size() != b.row) {
      std::stringstream msg;
      msg << "mulSquare : Invalid # of rows and size. " << a.size() << ":" << b.row << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    std::vector<float> vec;
#if !defined(NGT_DISABLE_BLAS)
    for (size_t c = 0; c < a.size(); c++) {
      T sum = 0;
      for (size_t p = 0; p < b.col; p++) {
	sum += a[p] * b.matrix[c * b.row + p];
      }
      vec.push_back(sum);
    }
#else
    for (size_t c = 0; c < a.size(); c++) {
      T sum = 0;
      for (size_t p = 0; p < b.col; p++) {
	sum += a[p] * b.matrix[p * b.col + c];
      }
      vec.push_back(sum);
    }
#endif
    a = vec;

  }

  static void mulSquare(std::vector<std::vector<float>> &a, Matrix<T> &b) {
    assert(b.col == b.row);
    for (size_t r = 0; r < a.size(); r++) {
      mulSquare(a[r], b);
    }
  }

  static void gesvd(char* jobu, char* jobvt, int* m, int* n, double* a,
		    int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
		    double* work, int* lwork, int* info) {
    dgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
  }

  static void gesvd(char* jobu, char* jobvt, int* m, int* n, float* a,
		    int* lda, float* s, float* u, int* ldu, float* vt, int* ldvt,
		    float* work, int* lwork, int* info) {
    sgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
  }

  static void gemm(char *transa, char *transb, int *m, int *n, int *k,
		   double *alpha, double *a, int *lda, double *b, int *ldb,
		   double *beta , double *c, int *ldc) {
    dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta , c, ldc);
  }

  static void gemm(char *transa, char *transb, int *m, int *n, int *k,
		   float *alpha, float *a, int *lda, float *b, int *ldb,
		   float *beta , float *c, int *ldc) {
    sgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta , c, ldc);
  }

  static void geqrf(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info) {
    dgeqrf_(m, n, a, lda, tau, work, lwork, info);
  }

  static void geqrf(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info) {
    sgeqrf_(m, n, a, lda, tau, work, lwork, info);
  }

  static void orgqr(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info) {
    dorgqr_(m, n, k, a, lda, tau, work, lwork, info);
  }

  static void orgqr(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info) {
    sorgqr_(m, n, k, a, lda, tau, work, lwork, info);
  }

  static void svd(Matrix<T> &a, Matrix<T> &u, Matrix<T> &s, Matrix<T> &v) {
    Matrix svda(a);
#if !defined(NGT_DISABLE_BLAS)
    int m = svda.row;
    int n = svda.col;
#else
    svda.transpose();
    int n = svda.row;
    int m = svda.col;
#endif
    char jobu = 'A';
    char jobvt = 'A';
    int lda = m, ldu = m, ldvt = n, info;
    int max, min;
    if (m > n) {
      max = m;
      min = n;
    } else {
      max = n;
      min = m;
    }
    int v1 = 3 * min + max;
    int v2 = 5 * min;
    int lwork = v1 > v2 ? v1 : v2;
    T work[lwork];

    Matrix<T> sd;
    sd.allocate(m, 1);
    u.allocate(m, m);
    v.allocate(n, n);
                                                    // S       U               VT
    gesvd(&jobu, &jobvt, &m, &n, svda.matrix, &lda, sd.matrix, u.matrix, &ldu, v.matrix, &ldvt, work, &lwork, &info);
    s.diag(sd);
    s.reshape(m, n);
#if !defined(NGT_DISABLE_BLAS)
    v.transpose();
#else
    u.transpose();
#endif
  }

  void eye(size_t d) {
    zero(d, d);
    for (size_t p = 0; p < row; p++) {
      matrix[p * col + p] = 1.0;
    }
  }

  void randomRotation(size_t d) {
    random(d, d);
    qr(d);
  }

  void qr(size_t d) {
    auto di = static_cast<int>(d);
    T tau[di];
    std::vector<T> work(1);
    int lwork = -1;
    int info;
    geqrf(&di, &di, matrix, &di, &tau[0], work.data(), &lwork, &info);
    work.resize(static_cast<int>(work[0]));

    geqrf(&di, &di, matrix, &di, tau, work.data(), &lwork, &info);
    orgqr(&di, &di, &di, matrix, &di, tau, work.data(), &lwork, &info);
  }

  void printmat() {
    T mtmp;
    printf("[ ");
    for (size_t i = 0; i < row; i++) {
      printf("[ ");
      for (size_t j = 0; j < col; j++) {
        mtmp = matrix[i + j * col];
        printf("%5.2e", mtmp);
        if (j < col - 1) printf(", ");
      }
      if (i < row - 1) printf("]; ");
      else printf("] ");
    }
    printf("]");
    std::cout << std::endl;
  }

  static void save(const std::string &file, const Matrix<T> &m) {
    std::ofstream os(file);
    for (size_t r = 0; r < m.row; r++) {
      for (size_t c = 0; c < m.col; c++) {
#if !defined(NGT_DISABLE_BLAS)
	os << m.matrix[c * m.row + r];
#else
        os << m.matrix[r * m.col + c];
#endif
	if (c + 1 != m.col) {
	  os << "\t";
        }
      }
      os << std::endl;
    }
  }

  static void
    convert(std::vector<std::string> &strings, std::vector<T> &vector) {
    vector.clear();
    for (auto it = strings.begin(); it != strings.end(); ++it) {
      try {
	vector.push_back(stod(*it));
      } catch(...) {
	break;
      }
    }
  }

  static void
    extractVector(const std::string &str, std::vector<T> &vec)
  {
    std::vector<std::string> tokens;
    tokenize(str, tokens, " \t");
    convert(tokens, vec);
  }

#if !defined(NGT_DISABLE_BLAS)
  static
    void load(const std::string &file, Matrix<T> &m)
  {
    loadVectors(file, m);
    m.transpose();
  }

  static
    void loadVectors(const std::string &file, Matrix<T> &m)
#else
  static
    void load(const std::string &file, Matrix<T> &m)
#endif
  {
    std::ifstream is(file);
    if (!is) {
      std::stringstream msg;
      msg << "Matrix::load: Cannot load. " << file;
      throw std::runtime_error(msg.str().c_str());
    }
    std::string line;
    size_t row = 0, col = 0;
    std::vector<T> tmpv;
    while (getline(is, line)) {
      std::vector<T> v;
      extractVector(line, v);
#if !defined(NGT_DISABLE_BLAS)
      if (row == 0) {
	row = v.size();
      } else if (row != v.size()) {
	std::cerr << "somthing wrong." << std::endl;
	abort();
      }
      col++;
#else
      if (col == 0) {
	col = v.size();
      } else if (col != v.size()) {
	std::cerr << "somthing wrong." << std::endl;
	abort();
      }
      row++;
#endif
      for (size_t i = 0; i < v.size(); i++) {
	tmpv.push_back(v[i]);
      }
    }
    m.construct(row, col, &tmpv[0]);
  }

  friend std::ostream& operator<<(std::ostream &os, const Matrix<T> &m) {
    os << m.row << " x " << m.col << "=" << std::endl;
    os << "[";
    for (size_t r = 0; r < m.row; r++) {
      os << r << ":[";
      for (size_t c = 0; c < m.col; c++) {
#if !defined(NGT_DISABLE_BLAS)
	os << m.matrix[c * m.row + r];
#else
	os << m.matrix[r * m.col + c];
#endif
	if (c + 1 != m.col) {
	  os << "\t";
        }
      }
      os << "]";
      if (r + 1 != m.row) {
	os << std::endl;
      }
    }
    os << "]";
    return os;
  }

  size_t row;
  size_t col;
  T *matrix;
};


