#pragma once
#include <iostream>
#include <vector>

using namespace std;

// ============================================== conjugate gradient ====================================================

constexpr double tol = 1.0e-10; // 1e-10f;



inline void hostCusparseSpMV(int N, vector<int> &h_row, vector<int> &h_col, vector<double> &h_val, vector<double> &h_vec, vector<double> &h_result)
{

  for (int i = 0; i < N; i++)
  {
    h_result[i] = 0.0;

    for (int j = h_row[i]; j < h_row[i + 1]; j++)
    {
      h_result[i] += h_val[j] * h_vec[h_col[j]];
    }
  }
}

// cublasDaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
// cublasDscal(cublasHandle, N, &b, d_p, 1);
inline void hostDaxpy(int N, double alpha, double beta, vector<double> &h_a, vector<double> &h_b)
{
  for (int i = 0; i < N; i++)
  {
    h_b[i] = alpha * h_a[i] + beta * h_b[i];
  }
}

// cublasDcopy(cublasHandle, N, d_r, 1, d_p, 1)
inline void hostDcopy(int N, vector<double> &h_a, vector<double> &h_b)
{
  hostDaxpy(N, 1.0, 0.0, h_a, h_b);
}

// cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
inline void hostDdot(int N, vector<double> &h_a, vector<double> &h_b, double &result)
{
  result = 0.0;
  for (int i = 0; i < N; i++)
  {
    result += h_b[i] * h_a[i];
  }
}

int cg_host(int N, int nz, vector<int> &h_row, vector<int> &h_col, vector<double> &h_val, vector<double> &h_b, vector<double> &h_x, const int print_n)
{
  vector<double> h_r(N);
  hostDcopy(N, h_b, h_r);
  vector<double> h_p(N);
  vector<double> h_Ap(N);
  double r0, r1, dot, alpha, beta;
  int k = 1;

  r0 = 0.0;
  hostCusparseSpMV(N, h_row, h_col, h_val, h_x, h_Ap);
  hostDaxpy(N, -1.0, 1.0, h_Ap, h_r);
  hostDdot(N, h_r, h_r, r1);

  clock_t startCheckTime = clock();

  while (r1 > (tol * tol) && k <= N)
  {
    if (k > 1)
    {
      beta = r1 / r0;
      hostDaxpy(N, 1.0, beta, h_r, h_p);
    }
    else
    {
      hostDcopy(N, h_r, h_p);
    }
    hostCusparseSpMV(N, h_row, h_col, h_val, h_p, h_Ap);
    hostDdot(N, h_p, h_Ap, dot);
    alpha = r1 / dot;
    hostDaxpy(N, alpha, 1.0, h_p, h_x);
    hostDaxpy(N, -alpha, 1.0, h_Ap, h_r);
    r0 = r1;
    hostDdot(N, h_r, h_r, r1);
    if (k % print_n == 0 || k == N || r1 <= tol * tol)
      printf("iteration = %3d, residual = %e, r1/r0 = %e \n", k, sqrt(r1), r1 / r0);
    k++;
  }
  printf("\ntotal time for %s is %lf \n", "CG without preconditioning", (double)(clock() - startCheckTime) / CLOCKS_PER_SEC);

  double err = 0.0;

  for (int i = 0; i < N; i++)
  {
    auto value = fabs(h_r[i]);
    if (value > err)
      err = value;
  }
  printf("\nTest Summary:  Error amount = %.6e\n", err);
  return (k <= N && err < tol) ? 0 : 1;
}

inline void genTridiag(vector<int> &I, vector<int> &J, vector<double> &val, int N, int &nz)
{
  nz = (N - 2) * 3 + 4;

  I.resize(N + 1);
  J.resize(nz);
  val.resize(nz);

  I[0] = 0, J[0] = 0, J[1] = 1;
  val[0] = (double)rand() / RAND_MAX + 10.0f;
  val[1] = (double)rand() / RAND_MAX;
  int start;

  for (int i = 1; i < N; i++)
  {
    if (i > 1)
    {
      I[i] = I[i - 1] + 3;
    }
    else
    {
      I[1] = 2;
    }

    start = (i - 1) * 3 + 2;
    J[start] = i - 1;
    J[start + 1] = i;

    if (i < N - 1)
    {
      J[start + 2] = i + 1;
    }

    val[start] = val[start - 1];
    val[start + 1] = (double)rand() / RAND_MAX + 10.0f;

    if (i < N - 1)
    {
      val[start + 2] = (double)rand() / RAND_MAX;
    }
  }

  I[N] = nz;
}

int cgTest()
{
  int N = 1048576, nz;
  vector<int> I;
  vector<int> J;
  vector<double> val;
  vector<double> x(N);
  vector<double> rhs(N);

  /* Generate a random tridiagonal symmetric matrix in CSR format */
  genTridiag(I, J, val, N, nz);

  for (int i = 0; i < N; i++)
  {
    rhs[i] = 1.0;
    x[i] = 0.0;
  }

  const int status = cg_host(N, nz, I, J, val, rhs, x, 1);

  return status;
}

// ============================================== helper functions ====================================================
void pd(double v)
{
  if (v < 0)
    printf("%-15.3lf", v);
  else
    printf(" %-14.3lf", v);
}
void pe(double v)
{
  if (v < 0)
    printf("%-15.3e", v);
  else
    printf(" %-14.3e", v);
}
void printABCDet(int ie, int edim, double det, vector<double> a, vector<double> b, vector<double> c)
{
  printf("\n\n%s %d %s %-.3lf\n", "elem", ie, "det", det);
  printf("%-15s %-15s %-15s\n", "a", "b", "c");
  for (int i = 0; i < edim; i++)
  {
    pd(a[i]);
    pd(b[i]);
    pd(c[i]);
    printf("\n");
  }
}

void printTensorH()
{
  printf("\n\n%-14s %-14s %-14s %-14s %-14s %-14s %-14s\n", "elem", "e11", "e22", "g12", "s11", "s22", "s12");
}
void printTensor(int ie, vector<double> epsilon, vector<double> sigma)
{
  printf("%-14d", ie);
  for (int i = 0; i < epsilon.size(); i++)
  {
    pe(epsilon[i]);
  }
  for (int i = 0; i < epsilon.size(); i++)
  {
    pe(sigma[i]);
  }
  printf("\n");
}

void printStiffness(int N, vector<double> stiffness)
{
  // printf("matrix %dx%d\n", N, N);
  for (int i = 0; i < stiffness.size(); i++)
  {
    if (stiffness[i] < 0)
      printf("%-12.1lf", stiffness[i]);
    else
      printf(" %-11.1lf", stiffness[i]);

    if ((i + 1) % N == 0)
    {
      printf("\n");
    }
  }
  printf("\n");
}

void printStiffness3(int N, vector<double> stiffness)
{
  // printf("matrix %dx%d\n", N, N);
  for (int i = 0; i < stiffness.size(); i++)
  {
    if (stiffness[i] < 0)
      printf("%-12.3lf", stiffness[i]);
    else
      printf(" %-11.3lf", stiffness[i]);

    if ((i + 1) % N == 0)
    {
      printf("\n");
    }
  }
  printf("\n");
}

void pd12(double v)
{
  if (v < 0)
    printf("%-12.1lf", v);
  else
    printf(" %-11.1lf", v);
}
void printStiffnessF(int N, vector<double> stiffness, vector<double> f)
{
  printf("matrix %dx%d\n", N, N);
  int j = 0;
  for (int i = 0; i < N * N; i++)
  {
    if (stiffness[i] == 0)
      printf("%-12.4s", " 0");
    else
      pd12(stiffness[i]);

    if ((i + 1) % N == 0)
    {
      if (f.size() > 0)
      {
        printf("|");
        pd(f[j]);
        printf("\n");
      }
      else
        printf("\n");
      j++;
    }
  }
  printf("\n");
}

void printSCR(int M, int N, int nz, vector<double> val, vector<int> I, vector<int> J)
{

  printf("\nnz=%d ", nz);
  printf("I[%d] J[%d] val[%d]\n", M + 1, nz, nz);
  printf("%-10s%-10s%-10s\n", "I", "J", "val");
  for (int i = 0; i < nz; i++)
  {
    if (i < M + 1)
    {
      printf("%-10d%-10d", I[i], J[i]);
      pd(val[i]);
      printf("\n");
    }
    else
    {
      printf("%-10s%-10d", " ", J[i]);
      pd(val[i]);
      printf("\n");
    }
  }
  printf("\n\n");
}

void resultPrint(int nsize, int ndim, int udim, vector<double> node, vector<double> result, string str = "u")
{
  printf("\n%-7s", "i");
  for (int i = 0; i < ndim; i++)
  {
    printf("x%-9i", i);
  }
  for (int i = 0; i < udim; i++)
  {
    printf("%s%-14i", str.c_str(), i);
  }
  printf("\n");
  for (int p = 0; p < nsize; p++)
  {
    printf("%-6d", p);
    for (int i = 0; i < ndim; i++)
    {
      if (node[ndim * p + i] < 0.0)
        printf("%-10.3lf", node[ndim * p + i]);
      else
        printf(" %-09.3lf", node[ndim * p + i]);
    }
    for (int i = 0; i < udim; i++)
    {

      if (result[udim * p + i] < 0.0)
        printf("%-15.5e", result[udim * p + i]);
      else
        printf(" %-14.5e", result[udim * p + i]);
    }
    printf("\n");
  }
  printf("\n");
}

void axialForcesPrint(int nsize, int ndim, int udim, vector<double> node, vector<int> P, vector<double> result, string str = "F")
{
  printf("\n%-7s", "P");
  for (int i = 0; i < ndim; i++)
  {
    printf("x%-9i", i);
  }
  for (int i = 0; i < udim; i++)
  {
    printf("%s%-14i", str.c_str(), i);
  }
  printf("\n");
  for (int p = 0; p < nsize; p++)
  {
    printf("%-6d", P[p]);
    for (int i = 0; i < ndim; i++)
    {
      if (node[ndim * P[p] + i] < 0.0)
        printf("%-10.3lf", node[ndim * P[p] + i]);
      else
        printf(" %-09.3lf", node[ndim * P[p] + i]);
    }
    for (int i = 0; i < udim; i++)
    {
      if (result[udim * p + i] < 0.0)
        printf("%-15.5e", result[udim * p + i]);
      else
        printf(" %-14.5e", result[udim * p + i]);
    }
    printf("\n");
  }
  printf("\n");
}