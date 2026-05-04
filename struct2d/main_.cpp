#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <cctype>
#include <functional> // std::plus // isinf()

#include <stdlib.h>

#include "helper.h"

using namespace std;

const int schemaNO = 0;

const int edim = 2;
const int ndim = 2;
const int udim = 3;
const int kdim = udim * edim;

const vector<int> elem =
{
  0, 1,
  1, 2,
};

const int esize = elem.size() / edim;

vector<double> node =
{// [cm]
0, 0,
282.842712474619, 282.842712474619,
682.842712474619, 282.842712474619
};

const int nsize = node.size() / ndim;

const vector<double> eparam = {// E|N/cm2|, A|cm2|, J|cm4| by element
  20e+6, 30, 400,
  20e+6, 30, 400,
};
const int pdim = 3;


vector<double> load =
{ // [ux,uy,rz]
  0.0,
  0.0,
  0.0,
  0.0,
  -100000,
  25000,
  0.0,
  0.0,
  0.0,
};

vector<double> bound =
{ // [cm]
  0.0,
  0.0,
  0.0,
  INFINITY,
  INFINITY,
  INFINITY,
  0.0,
  0.0,
  0.0,
};

void getLT(vector<double> a, vector<double> b, double &L, vector<double> &T)
{
  // const double dx = ;
  // const double dy = ;
  // L = ;
  // const double ix = ;
  // const double iy = ;
  // const double iz = ;  
  // const double jx = ;
  // const double jy = ;
  // const double jz = ;    
  /*
    ix iy iz  0  0  0
    jx jy jz  0  0  0
     0  0  1  0  0  0 
     0  0  0  ix iy iz
     0  0  0  jx jy jz
     0  0  0   0  0  1 
  */
  // T[0 * kdim + 0] = ix; T[0 * kdim + 1] = iy; T[0 * kdim + 2] = iz;  
}

void setLocalStiffness(double L, double E, double A, double J, vector<double> T, vector<double> &localStiffness){
  //
  // const double kN = A * L2 / J * k1;
  //
  for (int i = 0; i < localStiffness.size(); i++)
  {
    localStiffness[i] = 0;
  }
  // N
  // localStiffness[0 * kdim + 0] = kN; localStiffness[0 * kdim + 3] =-kN;
  // localStiffness[3 * kdim + 0] =-kN; localStiffness[3 * kdim + 3] = kN;
};

void multiplyAB(const int N, vector<double> mA,vector<double> mB,vector<double> &result, bool isTranspA){// square matrix !!!

};

void multiplyAb(const int N, vector<double> mA,vector<double> b,vector<double> &result, bool isTranspA){// square matrix !!!

};

void applyKinematic(int N, vector<double> u, vector<double> &f, vector<double> &stiffness)
{
  for (int p = 0; p < N; p++)
  {
    for (int q = 0; q < N; q++)
    {
      if (isinf(u[p]))
        continue;
      if (p == q)
      {
        f[p] = stiffness[p * N + q] * u[p];
      }
      else
      {
        stiffness[p * N + q] = 0;
        f[q] -= stiffness[q * N + p] * u[p];
        stiffness[q * N + p] = 0;
      }
    }
  }
  printf("\nKinematic boundary condition is applied\n\n");
}

int kinematicTest()
{ // Л.Сегерлинд. Применение метода конечных элементов стр. 110-112
  vector<double> u = {
      150,
      INFINITY,
      INFINITY,
      INFINITY,
      40,
  };
  vector<double> f = {
      500,
      2000,
      1000,
      2000,
      900,
  };
  vector<double> stiffness = {
      55, -46,   4,   0,   0,
     -46, 140, -46,   0,   0,
       4, -46, 110, -46,   4,
       0,   0, -46, 142, -46,
       0,   0,   4, -46,  65,
  };
  const int N = u.size();
  printStiffnessF(N, stiffness, f);
  applyKinematic(N, u, f, stiffness);
  printStiffnessF(N, stiffness, f);

  double sum =
      stiffness[0 * N + 0] - 55.0 + f[0] - 8250 +
      stiffness[1 * N + 1] - 140.0 + stiffness[1 * N + 2] + 46.0 + f[1] - 8900 +
      stiffness[2 * N + 1] + 46.0 + stiffness[2 * N + 2] - 110.0 + stiffness[2 * N + 3] + 46.0 + f[2] - 240 +
      stiffness[3 * N + 2] + 46.0 + stiffness[3 * N + 3] - 142.0 + f[3] - 3840 +
      stiffness[4 * N + 4] - 65.0 + f[4] - 2600;

  return sum == 0.0 ? 0 : 1;
}

// Преобразовать вектор-матрицу жёсткости в формат CSR https://docs.nvidia.com/cuda/cusparse/index.html#csr-format
void transformMatrixToCsr(int M, int N, int &nz, vector<double> &stiffness, vector<double> &val, vector<int> &I, vector<int> &J)
{

}

void transformMatrixTest()
{
  const int M = 4;
  const int N = 5;

  vector<double> stiffness = {
      1.0, 4.0, 0.0, 0.0, 0.0,
      0.0, 2.0, 3.0, 0.0, 0.0,
      5.0, 0.0, 0.0, 7.0, 8.0,
      0.0, 0.0, 9.0, 0.0, 6.0,
  };

  vector<double> val;
  vector<int> I;
  vector<int> J;
  int nz;
  // transformMatrixToCsr(M, N, nz, stiffness, val, I, J);
  // printSCR(M, N, nz, val, I, J);
}

int main(int argc, char **argv)
{

  const int N = nsize * udim;

  vector<double> globalStiffness(N * N);
  for (int i = 0; i < N * N; i++)
  {
    globalStiffness[i] = 0;
  }

  vector<int> P{0, 0}, Q{0, 0};
  vector<double> T{
    0,0,0,0,0,0,
    0,0,0,0,0,0,
    0,0,0,0,0,0,
    0,0,0,0,0,0,
    0,0,0,0,0,0,
    0,0,0,0,0,0,
  };
  vector<double> localStiffness{
    0,0,0,0,0,0,
    0,0,0,0,0,0,
    0,0,0,0,0,0,
    0,0,0,0,0,0,
    0,0,0,0,0,0,
    0,0,0,0,0,0,
  };
  vector<double> a{0, 0};
  vector<double> b{0, 0};   
  vector<double> elemStiffness(kdim * kdim);
  printf("%s\n\n", "struct elements...");   
  for (int ie = 0; ie < esize; ie++)
  {
    const double E = 1;//eparam[...];
    const double A = 1;//eparam[...];
    const double J = 1;//eparam[...];

    P[0] = 0;//elem[...];
    P[1] = 0;//elem[...];

    Q[0] = 0; //elem[...];
    Q[1] = 0; //elem[...];

    double L = 1;
    a[0] = 0;//node[...];
    a[1] = 0;//node[...];
    b[0] = 1;//node[...];
    b[1] = 1;//node[...];    

    printf("ie=%d; ", ie);      
    getLT(a, b, L, T);
    printf("L=%lf;\n", L);
    printf("matrixTransformation\n");  
    printStiffness(kdim, T);

    printf("localStiffness [K] (in local basis)\n"); 
    // setLocalStiffness.....
    printStiffness(kdim, localStiffness);

    printf("localStiffness [T]*[K][T] (in global basis)\n"); 
    // multiplyAB... ....
    printStiffness(kdim, localStiffness);

    // apply
    // globalStiffness[...] += localStiffness[...];
  }
  printf("global stiffness K=\n");
  printStiffnessF(N, globalStiffness, load);
  // applyKinematic...
  printf("global stiffness K with kinematic\n");
  printStiffnessF(N, globalStiffness, load);

  int status = kinematicTest();

  if (status != 0)
  {
    printf("Kinematic error\n");
    int exit(status);
  }
  else
  {
    printf("Kinematic test ok\n");
  }



  vector<double> val;
  vector<int> I;
  vector<int> J;
  int nz;
  transformMatrixTest();
  // transformMatrixToCsr(N, N, nz, globalStiffness, val, I, J);
  // printSCR(N, N, nz, val, I, J);

  vector<double> result(N);
  for (int i = 0; i < N; i++)
  {
    result[i] = 0.0;
  }

  // status = cg_host(N, nz, I, J, val, load, result, 1);

  printf("kinematic result vector\n"); 
  resultPrint(nsize, ndim, udim, node, result);

  vector<double> aforces(kdim);
  vector<double> elemu(kdim);
  printf("axial forces by element\n"); 
  for (int ie = 0; ie < esize; ie++)
  {
    const double E = 1;//eparam[...];
    const double A = 1;//eparam[...];
    const double J = 1;//eparam[...];

    Q[0] = 0;//elem[...];
    Q[1] = 0;//elem[...];

    double L = 1;
    a[0] = 0;//node[...];
    a[1] = 0;//node[...];
    b[0] = 1;//node[...];
    b[1] = 1;//node[...];    

    printf("ie=%d; ", ie);      
    getLT(a, b, L, T);
    // setLocalStiffness...
    // transformation ...

    // elemu[...] = result[...];
    // multiplyAb...;
    // aforces[0] = ...;
    // aforces[2] = ...;
    // aforces[4] = ...;
    // axialForcesPrint(edim, ndim, udim, node, Q, aforces, "f");
  }

  std::cin.get();
  int exit(status);
}