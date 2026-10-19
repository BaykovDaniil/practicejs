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

#include "struct2d/helper.h"

using namespace std;

const int schemaNO = 0;

const int edim = 2;
const int ndim = 2;
const int udim = 3;
const int kdim = udim * edim;

const vector<int> elem =
{
  0, 2,
  2, 3,
  3, 4,
  4, 5,
  5, 1,
};

vector<double> node =
{// [cm]
  0.0,   0.0,
  900.0, 0.0,
  0.0,   600.0,
  300.0, 600.0,
  600.0, 600.0,
  900.0, 600.0,
};

const double E0 = 20e+6;    // [N/cm^2]
const double F0 = 100e+3;   // [N]
const double M0 = 25e+3;    // [N*cm]
const double A0 = 30.0;     // [cm^2]
const double J0 = 400.0;    // [cm^4]

const vector<double> eparam = {// E|N/cm2|, A|cm2|, J|cm4| by element
  E0, A0, J0,
  E0, A0, J0,
  E0, A0, J0,
  E0, A0, J0,
  E0, A0, J0,
};

vector<double> load =
{ // [ux, uy, rz] 
  0.0, 0.0, 0.0, // node 0
  0.0, 0.0, 0.0, // node 1
  0.0, 0.0, 0.0, // node 2
  0.0, -F0, M0,  // node 3
  0.0, -F0, M0,  // node 4
  0.0, 0.0, 0.0, // node 5
};

vector<double> bound =
{ // [ux, uy, rz]
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  INFINITY, INFINITY, INFINITY,
  INFINITY, INFINITY, INFINITY,
  INFINITY, INFINITY, INFINITY,
  INFINITY, INFINITY, INFINITY,
};

const int ddim = 2;

// q = q0  [N/cm]      transverse distributed force
// m = m0  [N]         distributed bending moment, because [N*cm/cm] = [N]
//
// Row per element:
// q0, m0
const vector<double> dload =
{
  // element 0
  -200000.0 / 600.0, 0.0,

  // element 1
  0.0, 0.0,

  0.0, 0.0,

  0.0, 0.0,

  0.0, 0.0,

};


void printElementBendingValues(
  int ie,
  double L,
  double q0,
  double m0,
  const vector<double> &localDof,
  const vector<double> &endForces,
  double normalSign
)
{
  printf("\nM, Q, v by element, ie=%d\n", ie);
  printf("%10s %18s %18s %18s\n", "x/L", "v", "M", "Q");

  const double xis[3] = {0.0, 0.5, 1.0};

  for (int ip = 0; ip < 3; ++ip)
  {
    const double xi = xis[ip];
    const double x = xi * L;
    const double xi2 = xi * xi;
    const double xi3 = xi2 * xi;

    const double N1 = 2.0 * xi3 - 3.0 * xi2 + 1.0;
    const double N2 = L * (xi3 - 2.0 * xi2 + xi);
    const double N3 = 3.0 * xi2 - 2.0 * xi3;
    const double N4 = L * (xi3 - xi2);

    const double v1   = localDof[1];
    const double phi1 = localDof[2];
    const double v2   = localDof[4];
    const double phi2 = localDof[5];

    const double v = normalSign * (N1 * v1 + N2 * phi1 + N3 * v2 + N4 * phi2);
    const double M = -endForces[2] + (endForces[1] - m0) * x + q0 * x * x / 2.0;
    const double Q = endForces[1] - 2.0 * m0 + q0 * x;

    printf("%10.3lf %18.8e %18.8e %18.8e\n", xi, v, M, Q);
  }
}


void setLocalEquivalentLoad(double L, double q0, double m0, vector<double> &localLoad)
{
  for (int i = 0; i < kdim; ++i)
  {
    localLoad[i] = 0.0;
  }

  const double L2 = L * L;

  /*
    Local element dofs:

      0 u1
      1 v1
      2 phi1
      3 u2
      4 v2
      5 phi2

      N1a = 1 - x/L
      N2a = x/L
  */

  /*
    Transverse load part:
      q(x) = q0

    Bending Hermite interpolation:

      v(x) = N1*v1 + N2*phi1 + N3*v2 + N4*phi2

    Equivalent nodal load:

      f_v = ∫ N^T q(x) dx
  */

  localLoad[1] += q0 * L / 2.0;
  localLoad[2] += q0 * L2 / 12.0;
  localLoad[4] += q0 * L / 2.0;
  localLoad[5] += -q0 * L2 / 12.0;

  localLoad[1] += -m0;
  localLoad[2] +=  0.0;
  localLoad[4] +=  m0;
  localLoad[5] +=  0.0;
}
// const vector<int> elem =
// {
//   0, 1,
//   1, 2,
// };


// vector<double> node =
// {// [cm]
// 0, 0,
// 282.842712474619, 282.842712474619,
// 682.842712474619, 282.842712474619
// };


// const vector<double> eparam = {// E|N/cm2|, A|cm2|, J|cm4| by element
//   20e+6, 30, 400,
//   20e+6, 30, 400,
// };

// vector<double> load =
// { // [ux,uy,rz]
//   0.0,
//   0.0,
//   0.0,
//   0.0,
//   -100000,
//   25000,
//   0.0,
//   0.0,
//   0.0,
// };

// vector<double> bound =
// { // [cm]
//   0.0,
//   0.0,
//   0.0,
//   INFINITY,
//   INFINITY,
//   INFINITY,
//   0.0,
//   0.0,
//   0.0,
// };
const int esize = elem.size() / edim;
const int nsize = node.size() / ndim;
const int pdim = 3;


void getLT(vector<double> a, vector<double> b, double &L, vector<double> &T)
{
  const double dx = b[0] - a[0];
  const double dy = b[1] - a[1];

  L = sqrt(dx * dx + dy * dy);

  const double ix = dx / L;
  const double iy = dy / L;
  const double iz = 0.0;

  const double jx = -iy;
  const double jy =  ix;
  const double jz = 0.0;

  for (int i = 0; i < kdim * kdim; ++i)
  {
    T[i] = 0.0;
  }

  /*
    ix iy iz  0  0  0
    jx jy jz  0  0  0
     0  0  1  0  0  0 
     0  0  0  ix iy iz
     0  0  0  jx jy jz
     0  0  0   0  0  1 
  */

  T[0 * kdim + 0] = ix; T[0 * kdim + 1] = iy; T[0 * kdim + 2] = iz;
  T[1 * kdim + 0] = jx; T[1 * kdim + 1] = jy; T[1 * kdim + 2] = jz;
  T[2 * kdim + 0] = 0.0; T[2 * kdim + 1] = 0.0; T[2 * kdim + 2] = 1.0;

  T[3 * kdim + 3] = ix; T[3 * kdim + 4] = iy; T[3 * kdim + 5] = iz;
  T[4 * kdim + 3] = jx; T[4 * kdim + 4] = jy; T[4 * kdim + 5] = jz;
  T[5 * kdim + 3] = 0.0; T[5 * kdim + 4] = 0.0; T[5 * kdim + 5] = 1.0;
}

double getVNormalSign(const vector<double> &a, const vector<double> &b)
{
  const double dx = b[0] - a[0];
  const double dy = b[1] - a[1];
  const double L = sqrt(dx * dx + dy * dy);
  const double jx = -dy / L;
  const double jy =  dx / L;
  const double normal = fabs(jx) > fabs(jy) ? jx : jy;

  return normal < 0.0 ? -1.0 : 1.0;
}

void setLocalStiffness(double L, double E, double A, double J, vector<double> T, vector<double> &localStiffness)
{
  for (int i = 0; i < localStiffness.size(); i++)
  {
    localStiffness[i] = 0.0;
  }

  const double kN  = E * A / L;

  const double k11 = 12.0 * E * J / (L * L * L);
  const double k12 =  6.0 * E * J / (L * L);
  const double k22 =  4.0 * E * J / L;
  const double k24 =  2.0 * E * J / L;

  // axial part: u1, u2
  localStiffness[0 * kdim + 0] =  kN;
  localStiffness[0 * kdim + 3] = -kN;
  localStiffness[3 * kdim + 0] = -kN;
  localStiffness[3 * kdim + 3] =  kN;

  // bending part: {v1, phi1, v2, phi2}
  localStiffness[1 * kdim + 1] =  k11;
  localStiffness[1 * kdim + 2] =  k12;
  localStiffness[1 * kdim + 4] = -k11;
  localStiffness[1 * kdim + 5] =  k12;

  localStiffness[2 * kdim + 1] =  k12;
  localStiffness[2 * kdim + 2] =  k22;
  localStiffness[2 * kdim + 4] = -k12;
  localStiffness[2 * kdim + 5] =  k24;

  localStiffness[4 * kdim + 1] = -k11;
  localStiffness[4 * kdim + 2] = -k12;
  localStiffness[4 * kdim + 4] =  k11;
  localStiffness[4 * kdim + 5] = -k12;

  localStiffness[5 * kdim + 1] =  k12;
  localStiffness[5 * kdim + 2] =  k24;
  localStiffness[5 * kdim + 4] = -k12;
  localStiffness[5 * kdim + 5] =  k22;
}

void multiplyAB(const int N, vector<double> mA, vector<double> mB, vector<double> &result, bool isTranspA)
{
  for (int i = 0; i < N * N; ++i)
  {
    result[i] = 0.0;
  }

  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      double sum = 0.0;

      for (int k = 0; k < N; ++k)
      {
        const double Aik = isTranspA ? mA[k * N + i] : mA[i * N + k];
        const double Bkj = mB[k * N + j];

        sum += Aik * Bkj;
      }

      result[i * N + j] = sum;
    }
  }
}

void multiplyAb(const int N, vector<double> mA, vector<double> b, vector<double> &result, bool isTranspA)
{
  for (int i = 0; i < N; ++i)
  {
    result[i] = 0.0;
  }

  for (int i = 0; i < N; ++i)
  {
    double sum = 0.0;

    for (int k = 0; k < N; ++k)
    {
      const double Aik = isTranspA ? mA[k * N + i] : mA[i * N + k];
      sum += Aik * b[k];
    }

    result[i] = sum;
  }
}


void applyKinematic(int N, vector<double> u, vector<double> &f, vector<double> &stiffness)
{
  for (int p = 0; p < N; p++)
  {
    if (std::isinf(u[p]))
    {
      continue;
    }

    const double up = u[p];
    const double diag = stiffness[p * N + p];

    // Move known displacement contribution to RHS of all free equations.
    for (int q = 0; q < N; q++)
    {
      if (q == p || !std::isinf(u[q]))
      {
        continue;
      }

      f[q] -= stiffness[q * N + p] * up;
    }

    // Zero row and column of constrained dof.
    for (int q = 0; q < N; q++)
    {
      stiffness[p * N + q] = 0.0;
      stiffness[q * N + p] = 0.0;
    }

    // Keep original diagonal and set RHS to diag * prescribed_value.
    stiffness[p * N + p] = diag;
    f[p] = diag * up;
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
  val.clear();
  I.clear();
  J.clear();

  I.resize(M + 1);
  nz = 0;
  I[0] = 0;

  for (int i = 0; i < M; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      const double aij = stiffness[i * N + j];

      if (aij != 0.0)
      {
        val.push_back(aij);
        J.push_back(j);
        ++nz;
      }
    }

    I[i + 1] = nz;
  }
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
  transformMatrixToCsr(M, N, nz, stiffness, val, I, J);
  printSCR(M, N, nz, val, I, J);
}

int main(int argc, char **argv)
{
  vector<double> localEquivalentLoad(kdim);
  vector<double> elemEquivalentLoad(kdim);
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
    const double E = eparam[ie * pdim + 0];
    const double A = eparam[ie * pdim + 1];
    const double J = eparam[ie * pdim + 2];

    Q[0] = elem[ie * edim + 0];
    Q[1] = elem[ie * edim + 1];

    P[0] = Q[0] * udim;
    P[1] = Q[1] * udim;

    a[0] = node[Q[0] * ndim + 0];
    a[1] = node[Q[0] * ndim + 1];

    b[0] = node[Q[1] * ndim + 0];
    b[1] = node[Q[1] * ndim + 1];

    double L = 1.0;

    printf("ie=%d; ", ie);
    getLT(a, b, L, T);
    printf("L=%lf;\n", L);

    printf("matrixTransformation\n");
    printStiffness(kdim, T);

    printf("localStiffness [K] (in local basis)\n");
    setLocalStiffness(L, E, A, J, T, localStiffness);
    printStiffness(kdim, localStiffness);

    printf("localStiffness [T]*[K][T] (in global basis)\n");

    vector<double> temp(kdim * kdim);

    // temp = Kloc * T
    multiplyAB(kdim, localStiffness, T, temp, false);

    // elemStiffness = T^T * Kloc * T
    multiplyAB(kdim, T, temp, elemStiffness, true);

    printStiffness(kdim, elemStiffness);

    setLocalEquivalentLoad(
      L,
      dload[ie * ddim + 0],
      dload[ie * ddim + 1],
      localEquivalentLoad
    );

    multiplyAb(kdim, T, localEquivalentLoad, elemEquivalentLoad, true);
    // assemble into global stiffness
    for (int ilocal = 0; ilocal < kdim; ++ilocal)
    {
      const int inodeLocal = ilocal / udim;
      const int idofLocal  = ilocal % udim;
      const int iglobal = Q[inodeLocal] * udim + idofLocal;

      load[iglobal] += elemEquivalentLoad[ilocal];

      for (int jlocal = 0; jlocal < kdim; ++jlocal)
      {
        const int jnodeLocal = jlocal / udim;
        const int jdofLocal  = jlocal % udim;
        const int jglobal = Q[jnodeLocal] * udim + jdofLocal;

        globalStiffness[iglobal * N + jglobal] += elemStiffness[ilocal * kdim + jlocal];
      }
    }
  }
  printf("global stiffness K=\n");
  printStiffnessF(N, globalStiffness, load);
  applyKinematic(N, bound, load, globalStiffness);
  printf("global stiffness K with kinematic\n");
  printStiffnessF(N, globalStiffness, load);

  int status = kinematicTest();

  if (status != 0)
  {
    printf("Kinematic error\n");
    return status;
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
  transformMatrixToCsr(N, N, nz, globalStiffness, val, I, J);
  printSCR(N, N, nz, val, I, J);

  vector<double> result(N);
  for (int i = 0; i < N; i++)
  {
    result[i] = 0.0;
  }

  status = cg_host(N, nz, I, J, val, load, result, 1);

  printf("kinematic result vector\n"); 
  resultPrint(nsize, ndim, udim, node, result);

  vector<double> aforces(kdim);
  vector<double> elemu(kdim);
  vector<double> elemuLocal(kdim);
  vector<double> rawLocalForces(kdim);

  printf("axial forces by element\n");

  for (int ie = 0; ie < esize; ie++)
  {
    const double E = eparam[ie * pdim + 0];
    const double A = eparam[ie * pdim + 1];
    const double J = eparam[ie * pdim + 2];

    Q[0] = elem[ie * edim + 0];
    Q[1] = elem[ie * edim + 1];

    a[0] = node[Q[0] * ndim + 0];
    a[1] = node[Q[0] * ndim + 1];

    b[0] = node[Q[1] * ndim + 0];
    b[1] = node[Q[1] * ndim + 1];

    double L = 1.0;

    printf("ie=%d;\n", ie);

    getLT(a, b, L, T);
    setLocalStiffness(L, E, A, J, T, localStiffness);

    // Extract global element displacement vector.
    for (int inode = 0; inode < edim; ++inode)
    {
      for (int idof = 0; idof < udim; ++idof)
      {
        const int localIndex = inode * udim + idof;
        const int globalIndex = Q[inode] * udim + idof;

        elemu[localIndex] = result[globalIndex];
      }
    }


    // global element displacement -> local element displacement
    multiplyAb(kdim, T, elemu, elemuLocal, false);

    // equivalent local load again for this element
    setLocalEquivalentLoad(
      L,
      dload[ie * ddim + 0],
      dload[ie * ddim + 1],
      localEquivalentLoad
    );

    // local nodal force vector
    multiplyAb(kdim, localStiffness, elemuLocal, rawLocalForces, false);

    // subtract equivalent distributed element load
    for (int i = 0; i < kdim; ++i)
    {
      rawLocalForces[i] -= localEquivalentLoad[i];
    }

    printElementBendingValues(
      ie,
      L,
      dload[ie * ddim + 0],
      dload[ie * ddim + 1],
      elemuLocal,
      rawLocalForces,
      getVNormalSign(a, b)
    );

    aforces[0] = -rawLocalForces[0];
    aforces[1] = rawLocalForces[1] - 2.0 * dload[ie * ddim + 1];
    aforces[2] = -rawLocalForces[2];

    aforces[3] = rawLocalForces[3];
    aforces[4] = -rawLocalForces[4] - 2.0 * dload[ie * ddim + 1];
    aforces[5] = rawLocalForces[5];

    axialForcesPrint(edim, ndim, udim, node, Q, aforces, "f");

  }

  return 0;
}
