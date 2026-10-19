import { Vector3 } from "three";
import { axialForcesPrint, cg, logm, printResult, printSCR } from "./helper.js";
// import { sqrt } from "three/tsl";

class Isoparametric2d {

  static edim = 4; // количество узлов элемента
  static ndim = 2; // степень свободы узла на плоскости 2Д
  static udim = 2; 
  static kdim = this.udim * this.edim;
  constructor() {
    this.edim = Isoparametric2d.edim;
    this.ndim = Isoparametric2d.ndim;
    this.udim = Isoparametric2d.udim;
    this.kdim = Isoparametric2d.kdim;

    this.pdim = 3;
    this.elem = [
      0, 1, 2, 3
    ];

    this.node = [// [cm]
      0, 0,
      0, 2,
      2, 2,
      2, 0
    ];

    this.elastic = [20.6, 0.3]; // { [MN/cm2], [1]}: модуль Юнга и коэфф. Пуассона
    this.thikness = 1.0; // [cm]: толщина
    this.alpha = 11.3e-06;
    this.thetta = 100;
    // this.eparam = [ // E|N/cm2|, A|cm2|, J|cm4| by element
    //   20e+6, 30, 400,
    // ];

    this.load = [ // [ux,uy]
        0.0, 0.0, 
        0.0, 0.0,
        0.0, -1000, 
        0.0, -1000
    ];

    // this.ieload = { // [py,mz]
    //     0: [10, 10],
    //     1: [10, 10]
    // };

    this.bound = [ // [cm]: Кинематические граничные условия
      0.0, 0.0,
      undefined, 0.0,
      undefined, undefined,
      0.0, undefined, 
    ];

    this.solve();
  }

  get esize() {
    return this.elem.length / this.edim;
  }

  get nsize() {
    return this.node.length / this.ndim;
  }

  getLT(a, b, L, T) {
    const dx = b[0] - a[0];
    const dy = b[1] - a[1];
    L = Math.sqrt(dx*dx + dy*dy);
    const ix = dx / L;
    const iy = dy / L;
    const iz = 0;  
    const jx = -iy;
    const jy = ix;
    const jz = 0;    
    /*
      ix iy iz  0  0  0
      jx jy jz  0  0  0
       0  0  1  0  0  0 
       0  0  0  ix iy iz
       0  0  0  jx jy jz
       0  0  0   0  0  1 
    */
    T[0 * this.kdim + 0] = ix; T[0 * this.kdim + 1] = iy; T[0 * this.kdim + 2] = iz;  
    T[1 * this.kdim + 0] = jx; T[1 * this.kdim + 1] = jy; T[1 * this.kdim + 2] = jz; 
    T[2 * this.kdim + 2] = 1;
    
    T[3 * this.kdim + 3] = ix; T[3 * this.kdim + 4] = iy; T[3 * this.kdim + 5] = iz;  
    T[4 * this.kdim + 3] = jx; T[4 * this.kdim + 4] = jy; T[4 * this.kdim + 5] = jz; 
    T[5 * this.kdim + 5] = 1;

  }

  setLocalStiffness(L, E, A, J, T, localStiffness) {
  
    const kN = E * A / Math.pow(L, 1); // ????? в записях L^2
    const a = E * J / Math.pow(L, 3);
  
    for (let i = 0; i < localStiffness.length; i++)
    {
      localStiffness[i] = 0;
    }

    localStiffness[1 * this.kdim + 1] = 12*a;
    localStiffness[1 * this.kdim + 2] = 6*L*a;
    localStiffness[1 * this.kdim + 4] = -12*a;
    localStiffness[1 * this.kdim + 5] = 6*L*a;

    localStiffness[2 * this.kdim + 1] = 6*L*a;
    localStiffness[2 * this.kdim + 2] = 4*L*L*a;
    localStiffness[2 * this.kdim + 4] = -6*L*a;
    localStiffness[2 * this.kdim + 5] = 2*L*L*a;

    localStiffness[4 * this.kdim + 1] = -12*a;
    localStiffness[4 * this.kdim + 2] = -6*L*a;
    localStiffness[4 * this.kdim + 4] = 12*a;
    localStiffness[4 * this.kdim + 5] = -6*L*a;

    localStiffness[5 * this.kdim + 1] = 6*L*a;
    localStiffness[5 * this.kdim + 2] = 2*L*L*a;
    localStiffness[5 * this.kdim + 4] = -6*L*a;
    localStiffness[5 * this.kdim + 5] = 4*L*L*a;
    
    localStiffness[0 * this.kdim + 0] = kN; localStiffness[0 * this.kdim + 3] =-kN;
    localStiffness[3 * this.kdim + 0] =-kN; localStiffness[3 * this.kdim + 3] = kN;
  }

  multiplyAB(N, mA, mB, result, isTranspA) {// square matrix !!!

    if(isTranspA === 0) {
        for(let i = 0; i < N; i++) {
            for(let j = 0; j < N; j++) {
                let sum = 0;
                for(let k = 0; k < N; k++) {
                    sum += mA[i*N + k] * mB[k*N + j];
                }
                result[i*N + j] = sum;
            }
        }
    }
    if(isTranspA === 1) {
        for(let i = 0; i < N; i++) {
            for(let j = 0; j < N; j++) {
                let sum = 0;
                for(let k = 0; k < N; k++) {
                    sum += mA[i + N*k] * mB[k*N + j];
                }
                result[i*N + j] = sum;
            }
        }
    }

  }

  multiplyAb(N, mA, b, result) {// square matrix !!!

    for(let i=0; i<N; i++) {
      let sum = 0;
      for(let j=0; j<N; j++) {
        sum += mA[i*N + j];
      }
      result[i] = sum;
    }

  }


  applyKinematic(N, u, f, stiffness) {

    for (let p = 0; p < N; p++)
    {
      if (u[p] !== undefined) 
      {
        for (let q = 0; q < N; q++) // 1
        {
          if (q === p) {
            f[p] = stiffness[p*N + p] * u[p];
          }
          else {
            stiffness[p*N + q] = 0;
          }
        }

        for (let q = 0; q < N; q++) // 2
        {
          if (q === p) {
            continue;
          }
          else {
            f[q] = f[q] - stiffness[q*N + p] * u[p];
          }
        }

        for (let q = 0; q < N; q++) // 3
        {
          if (q === p) {
            continue;
          }
          else {
            stiffness[q*N + p] = 0;
          }
        }
      }
    }
    console.log("\nKinematic boundary condition is applyed\n\n");

  }

  kinematicTest() { // Л.Сегерлинд. Применение метода конечных элементов стр. 110-112
    const u = [
      150,
      undefined,
      undefined,
      undefined,
      40,
    ];
    const f = [
      500,
      2000,
      1000,
      2000,
      900,
    ];
    const stiffness = [
      55, -46,   4,   0,   0,
      -46, 140, -46,   0,   0,
      4, -46, 110, -46,   4,
      0,   0, -46, 142, -46,
      0,   0,   4, -46,  65,
    ];
    const N = u.length;
    // this.printStiffnessF(N, stiffness, f);
    logm('matrix before', N, stiffness, f, 5, 0);
    
    this.applyKinematic(N, u, f, stiffness);
    // this.printStiffnessF(N, stiffness, f);
    logm('matrix after', N, stiffness, f, 5, 0);

    let sum =
      stiffness[0 * N + 0] - 55.0 + f[0] - 8250 +
      stiffness[1 * N + 1] - 140.0 + stiffness[1 * N + 2] + 46.0 + f[1] - 8900 +
      stiffness[2 * N + 1] + 46.0 + stiffness[2 * N + 2] - 110.0 + stiffness[2 * N + 3] + 46.0 + f[2] - 240 +
      stiffness[3 * N + 2] + 46.0 + stiffness[3 * N + 3] - 142.0 + f[3] - 3840 +
      stiffness[4 * N + 4] - 65.0 + f[4] - 2600;

    console.log(sum);

    return sum === 0.0 ? 0 : 1;
  }

  sumStiffnessTest(N, stiffness) {

    const tol = 1e-9;
    let sum = 0;

    // 

    console.log("\n\nSum stiffness is %5.9e\n\n", sum);

    return Math.abs(sum) < tol ? 0 : 1;
  }

  // Преобразовать вектор-матрицу жёсткости в формат CSR https://docs.nvidia.com/cuda/cusparse/index.html#csr-format
  transformMatrixToCsr(M, N, nz, stiffness, _val, _I, _J) {

    const val = []; // массив значений
    const I = new Array(N+1); // индексы начала элементов каждой строки
    const J = []; // номера столбцов элементов

    nz = 0;
    let icounter = 0;

    for (let i = 0; i < M; i++)
    {
      I[icounter] = nz;
      icounter++;
      for (let j = 0; j < N; j++) {

        if (stiffness[i * N + j] === 0) continue;

        val.push(stiffness[i * N + j]);
        J.push(j);
        nz++;
      }
    }
    I[icounter] = nz;

    _I.length = M + 1;
    _J.length = nz;
    _val.length = nz;

    for (let i = 0; i < nz; i++) {
      if (i < M + 1) {
        _I[i] = I[i];
      }
      _val[i] = val[i];
      _J[i] = J[i];
    }
  }

  transformMatrixTest() {
    const M = 4;
    const N = 5;

    const stiffness = [
        1.0, 4.0, 0.0, 0.0, 0.0,
        0.0, 2.0, 3.0, 0.0, 0.0,
        5.0, 0.0, 0.0, 7.0, 8.0,
        0.0, 0.0, 9.0, 0.0, 6.0
    ];

    const val = [];
    const I = [];
    const J = [];
    let nz;
    this.transformMatrixToCsr(M, N, nz, stiffness, val, I, J);
    printSCR(val, I, J, 2, 2);
    // this.printSCR_Matrix(M, N, nz, val, I, J);
    // console.log(val);
    // console.log(I);
    // console.log(J);
  }
  
  printSCR_Matrix(M, N, nz, val, I, J) {

    const matrix = new Array(N * M);
    for (let i = 0; i < N * N; i++) {
        matrix[i] = 0;
    }

    for (let i = 0; i < M; i++){
      const start = I[i];
      const end = I[i+1];

      for(let j = start; j < end; j++){
        const col = J[j];
        const value = val[j];
        matrix[i * N + col] = value;
      }
    }

    logm('Matrix from CSR', N, matrix, [], 8, 3);

  }

  resultPrint(result, N) {

    let x0 = Array(this.esize);
    let x1 = Array(this.esize);
    let x3 = Array(this.esize);

    let u0 = Array(this.esize);
    let u1 = Array(this.esize);
    let u2 = Array(this.esize);

    if( this.ndim === 2) {
      for(let i = 0; i < this.nsize; i++) {
        x0[i] = this.node[this.ndim*i + 0];
        x1[i] = this.node[this.ndim*i + 1];
        u0[i] = result[this.udim*i + 0];
        u1[i] = result[this.udim*i + 1];
        u2[i] = result[this.udim*i + 2];
      }
      const data = x0.map(
        (item, index) => ({
          'x0': item,
          'x1': x1[index],
          'u0': u0[index],
          'u1': u1[index],
          'u2': u2[index],
        })
      )
      console.table(data)
    }

    // if( ndim === 3) {
    //   for(let i = 0; i < nsize; i++) {
    //     x0[i] = node[ndim*i];
    //     x1[i] = node[ndim*i + 1];
    //     x2[i] = node[ndim*i + 2];
    //     u0[i] = result[ndim*i];
    //     u1[i] = result[ndim*i + 1];
    //     u2[i] = result[ndim*i + 2];
    //   }
    //   const data = x0.map(
    //     (item, index) => ({
    //       'x0': item,
    //       'x1': x1[index],
    //       'x2': x1[index],
    //       'u0': u0[index],
    //       'u1': u1[index],
    //       'u2': u1[index],
    //     })
    //   )
    //   console.table(data)
    // }

  }

  getForm(result, xi1, xi2) { // result -> N_p(xi1,xi2)

    result[0] = 1 / 4 * (1 - xi1)*(1 - xi2);
    result[1] = 1 / 4 * (1 + xi1)*(1 - xi2);
    result[2] = 1 / 4 * (1 + xi1)*(1 + xi2);
    result[3] = 1 / 4 * (1 - xi1)*(1 + xi2);

  }

  getFormDeriv(result, xi1, xi2) {

    result[0] = -1 / 4 * (1 - xi2);
    result[1] =  1 / 4 * (1 - xi2);
    result[2] =  1 / 4 * (1 + xi2);
    result[3] = -1 / 4 * (1 + xi2);

    result[4] = -1 / 4 * (1 - xi1);
    result[5] = -1 / 4 * (1 + xi1);
    result[6] =  1 / 4 * (1 + xi1);
    result[7] =  1 / 4 * (1 - xi1);

  }

  getJacMartix(result, xi1, xi2, x1, x2, x3, x4, y1, y2, y3, y4) {

    const deriv = new Array(8);
    this.getFormDeriv(deriv, xi1, xi2);
    const x = [x1, x2, x3, x4];
    const y = [y1, y2, y3, y4];

    // J = [dx/dxi1, dx/dxi2; dy/dxi1, dy/dxi2], построчно
    result[0] = 0;
    result[1] = 0;
    result[2] = 0;
    result[3] = 0;
    for (let i = 0; i < this.edim; i++) {
      result[0] += deriv[i] * x[i];
      result[1] += deriv[this.edim + i] * x[i];
      result[2] += deriv[i] * y[i];
      result[3] += deriv[this.edim + i] * y[i];
    }

  }

  getJacMartix_opt(result, deriv, x, y) {

    result[0] = 0;
    result[1] = 0;
    result[2] = 0;
    result[3] = 0;

    for (let p = 0; p < this.edim; p++) {
      result[0] += deriv[p] * x[p];
      result[1] += deriv[this.edim + p] * x[p];
      result[2] += deriv[p] * y[p];
      result[3] += deriv[this.edim + p] * y[p];
    }

  }

  getJacobian(jacMatrix) {
    return jacMatrix[0] * jacMatrix[3] - jacMatrix[1] * jacMatrix[2];
  }

  getJacInverse(result, jacMatrix) {

    const [a, b, c, d] = jacMatrix;
    const det = this.getJacobian(jacMatrix);

    result[0] = d / det;
    result[1] = -b / det;
    result[2] = -c / det;
    result[3] = a / det;

  }


  getPortrait(result) {

    result.length = this.nsize;

    for (let P = 0; P < this.nsize; P++) {
      result[P] = [];
    }

    for (let ie = 0; ie < this.esize; ie++) {
      for (let p = 0; p < this.edim; p++) {

        const P = this.elem[ie * this.edim + p];

        for (let q = 0; q < this.edim; q++) {

          const Q = this.elem[ie * this.edim + q];

          if (Q < P) {
            continue;
          }

          if (result[P].indexOf(Q) === -1) {
            result[P].push(Q);
          }
        }
      }
    }

    for (let P = 0; P < this.nsize; P++) {
      result[P].sort((a, b) => a - b);
    }

  }


  solve() {

    const portrait = [];
    this.getPortrait(portrait);
    console.log(portrait);

    const v = this.elastic[1];
    const E = this.elastic[0];

    const lambda = v*E / ( (1+v) * (1-2*v) ); // const lambda = v*E / (1-v*v); для ПН
    const mu = E / ( 2*(1+v));

    const N = this.nsize * this.udim;

    const globalStiffness = new Array(N * N);
    for (let i = 0; i < N * N; i++) {
        globalStiffness[i] = 0;
    }

    const load = this.load.slice();

    const gauss = [-1 / Math.sqrt(3), 1 / Math.sqrt(3)];
    const weight = [1, 1];

    const gaussDeriv = new Array(gauss.length * gauss.length);
    for (let i = 0; i < gauss.length; i++) {
      for (let j = 0; j < gauss.length; j++) {
        const g = i * gauss.length + j;

        gaussDeriv[g] = new Array(this.ndim * this.edim);
        this.getFormDeriv(gaussDeriv[g], gauss[i], gauss[j]);
      }
    }

    const P = new Array(this.edim);
    const x = new Array(this.edim);
    const y = new Array(this.edim);

    const deriv = new Array(this.ndim * this.edim);
    const dNdx = new Array(this.edim);
    const dNdy = new Array(this.edim);
    const jacMatrix = new Array(this.ndim * this.ndim);
    const jacInverse = new Array(this.ndim * this.ndim);

    const block = new Array(this.udim * this.udim);
    const elemStiffness = new Array(this.kdim * this.kdim);

    // this.printStiffness(N, globalStiffness);
    logm('globalStiffness matrix', N, globalStiffness, [], 3, 0);

    // Сборка матрицы
    for (let ie = 0; ie < this.esize; ie++) {

      for (let p = 0; p < this.edim; p++) {
        P[p] = this.elem[ie * this.edim + p];
        x[p] = this.node[P[p] * this.ndim + 0];
        y[p] = this.node[P[p] * this.ndim + 1];
      }

      for (let i = 0; i < this.kdim * this.kdim; i++) {
        elemStiffness[i] = 0;
      }

      // Интегрирование
      for (let i = 0; i < gauss.length; i++) {
        for (let j = 0; j < gauss.length; j++) {

          const deriv = gaussDeriv[i * gauss.length + j];

          this.getJacMartix_opt(jacMatrix, deriv, x, y);

          const det = this.getJacobian(jacMatrix);
          this.getJacInverse(jacInverse, jacMatrix);

          for (let p = 0; p < this.edim; p++) {
            dNdx[p] = jacInverse[0] * deriv[p]
                      + jacInverse[2] * deriv[this.edim + p];
            dNdy[p] = jacInverse[1] * deriv[p]
                      + jacInverse[3] * deriv[this.edim + p];
          }

          // const xi1 = gauss[i];
          // const xi2 = gauss[j];

          // this.getFormDeriv(deriv, xi1, xi2);
          // this.getJacMartix(jacMatrix, xi1, xi2,
          //   x[0], x[1], x[2], x[3],
          //   y[0], y[1], y[2], y[3]);

          // const det = this.getJacobian(jacMatrix);
          // this.getJacInverse(jacInverse, jacMatrix);

          // // Производные по x, y
          // for (let p = 0; p < this.edim; p++) {
          //   dNdx[p] = jacInverse[0] * deriv[p]
          //             + jacInverse[2] * deriv[this.edim + p];
          //   dNdy[p] = jacInverse[1] * deriv[p]
          //             + jacInverse[3] * deriv[this.edim + p];
          // }

          const k = this.thikness * Math.abs(det) * weight[i] * weight[j];

          for (let p = 0; p < this.edim; p++) {
            for (let q = 0; q < this.edim; q++) {

              block[0] = (lambda + 2*mu) * dNdx[p] * dNdx[q]
                          + mu * dNdy[p] * dNdy[q];
              block[1] = lambda * dNdx[p] * dNdy[q]
                          + mu * dNdy[p] * dNdx[q];
              block[2] = lambda * dNdy[p] * dNdx[q]
                          + mu * dNdx[p] * dNdy[q];
              block[3] = (lambda + 2*mu) * dNdy[p] * dNdy[q]
                          + mu * dNdx[p] * dNdx[q];

              for (let a = 0; a < this.udim; a++) {
                for (let b = 0; b < this.udim; b++) {
                  const row = p * this.udim + a;
                  const col = q * this.udim + b;
                  elemStiffness[row * this.kdim + col] += k * block[a * this.udim + b];
                }
              }
            }
          }
        }
      }

      // Сборка глобальной матрицы жёсткости
      for (let p = 0; p < this.edim; p++) {
        for (let q = 0; q < this.edim; q++) {
          for (let i = 0; i < this.udim; i++) {
            for (let j = 0; j < this.udim; j++) {
              const row = p * this.udim + i;
              const col = q * this.udim + j;
              const globalRow = P[p] * this.udim + i;
              const globalCol = P[q] * this.udim + j;

              globalStiffness[globalRow * N + globalCol] +=
                elemStiffness[row * this.kdim + col];
            }
          }
        }
      }
    }

    this.load = load;


    // const a = [0, 0, 0], b = [0, 0, 0], c = [0, 0, 0];
    // const P = [0, 0], Q = [0, 0];
    // const T = new Array(this.kdim * this.kdim);
    // const localStiffness = new Array(this.kdim * this.kdim);
    // for (let i = 0; i < this.kdim * this.kdim; i++) {
    //     T[i] = 0;
    //     localStiffness[i] = 0;
    // }

    // const a = [0, 0], b = [0, 0];
    // const elemStiffness = new Array(this.kdim * this.kdim);


    

  //   for (let ie = 0; ie < this.esize; ie++) {
  //   const E = this.eparam[3*ie + 0];//eparam[...];
  //   const A = this.eparam[3*ie + 1];//eparam[...];
  //   const J = this.eparam[3*ie + 2];//eparam[...];

  //   P[0] = this.elem[2*ie + 0];//elem[...];
  //   P[1] = this.elem[2*ie + 1];//elem[...];

  //   Q[0] = this.elem[2*ie + 0]; //elem[...];
  //   Q[1] = this.elem[2*ie + 1]; //elem[...];

  //   let L = 1;
  //   a[0] = this.node[2*P[0] + 0];//node[...];
  //   a[1] = this.node[2*P[0] + 1];//node[...];
  //   b[0] = this.node[2*P[1] + 0];//node[...];
  //   b[1] = this.node[2*P[1] + 1];//node[...];
  //   console.log(a);
  //   console.log(b);
  //   L = Math.sqrt( (a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]) );

  //   this.getLT(a, b, L, T);

  //   console.log(`elem ${ie}`);     
  //   console.log(`L = ${L}`);
  //   logm('Matrix Transformation', this.edim * this.udim, T, 5, 0);
 
  //   this.setLocalStiffness(L, E, A, J, T, localStiffness);
  //   logm('Local Stiffness [K] (in local basis)', this.kdim, localStiffness, 5, 0);

  //   const result_tmp = new Array(this.kdim * this.kdim);
  //   const result     = new Array(this.kdim * this.kdim);
  //   this.multiplyAB(this.kdim, T, localStiffness, result_tmp, 1);
  //   this.multiplyAB(this.kdim, result_tmp, T, result, 0);
  //   logm('localStiffness [T]*[K][T] (in global basis)', this.kdim, result, 5, 0);


  //   const dofPerNode = this.udim;
  //   const node1 = P[0];
  //   const node2 = P[1];


  //   const globalDofs = [
  //       node1 * dofPerNode + 0,  // u в узле 1
  //       node1 * dofPerNode + 1,  // v в узле 1
  //       node1 * dofPerNode + 2,  // z в узле 1
  //       node2 * dofPerNode + 0,  // u в узле 2
  //       node2 * dofPerNode + 1,   // v в узле 2
  //       node2 * dofPerNode + 2   // z в узле 2
  //   ];

    
    
  //   if (this.ieload && this.ieload.hasOwnProperty(ie)) {
  //     const q = this.ieload[ie][0] || 0; // поперечная нагрузка в локальной оси y
  //     const m = this.ieload[ie][1] || 0; // изгибающий момент вдоль z

    
  //     const feLocal = new Array(this.kdim).fill(0);
  //     const L2 = L * L;

  //   // Вклад от q
  //     const qL2 = q * L / 2.0;
  //     feLocal[1] += qL2;
  //     feLocal[4] += qL2;

  //     const qL2_12 = q * L2 / 12.0;
  //     feLocal[2] += qL2_12;
  //     feLocal[5] -= qL2_12;

  //   // Вклад от m
  //     const mL2 = m * L / 2.0;
  //     // feLocal[2] += mL2;
  //     // feLocal[5] += mL2;
  //     feLocal[1] -= m;
  //     feLocal[4] += m;

  //   // Переход к глобальной системе: F_global = T^T * F_local
  //     const feGlobal = new Array(this.kdim).fill(0);
  //     for (let i = 0; i < this.kdim; i++) {
  //       for (let j = 0; j < this.kdim; j++) {
  //         feGlobal[i] += T[j * this.kdim + i] * feLocal[j];
  //       }
  //     }

  //   // Добавление к глобальному вектору нагрузок
  //   for (let i = 0; i < this.kdim; i++) {
  //     const gdof = globalDofs[i];
  //     load[gdof] += feGlobal[i];
  //   }
  // }




  //   for (let i = 0; i < this.kdim; i++) {
  //       for (let j = 0; j < this.kdim; j++) {
  //           const globalRow = globalDofs[i];
  //           const globalCol = globalDofs[j];
  //           const index = globalRow * N + globalCol;
  //           globalStiffness[index] += result[i * this.kdim + j];
  //       }
  //   }
  // }

    // this.load = load;

    logm('globalStiffness matrix', N, globalStiffness, this.load, 8, 3);

    // let status = this.sumStiffnessTest(N, globalStiffness); // TODO

    // if (status !== 0) {
    //     console.log("\nSum stiffness error\n");
    //     return status;
    // } else {
    //     console.log("\nSum stiffness ok\n");
    // }

    
    // logm('globalStiffness matrix', N, globalStiffness, this.load, 8, 3);

    // status = this.kinematicTest();

    // if (status !== 0) {
    //     console.log("\nKinematic error\n");
    //     return status;
    // } else {
    //     console.log("\nKinematic test ok\n");
    // }

    // this.applyKinematic(N, this.bound, this.load, globalStiffness);

    // console.log(`Kinematic boundary condition is applied`); 
    // logm('global stiffness K with kinematic', N, globalStiffness, this.load, 8, 3);

    // // this.transformMatrixTest();

    // const val = [];
    // const I = [];
    // const J = [];
    // let nz;
    // this.transformMatrixToCsr(N, N, nz, globalStiffness, val, I, J);

    // const result = new Array(N).fill(0);
    // const { succesfull, error } = cg(N, I, J, val, this.load, 1, result);

    // this.resultPrint(result, N);

    // console.log(`axial forces by element`);

    // for (let ie = 0; ie < this.esize; ie++) {
    //   const E = this.eparam[3*ie + 0];//eparam[...];
    //   const A = this.eparam[3*ie + 1];//eparam[...];
    //   const J = this.eparam[3*ie + 2];//eparam[...];

    //   Q[0] = this.elem[2*ie + 0]; //elem[...];
    //   Q[1] = this.elem[2*ie + 1]; //elem[...];

    //   let L = 1;
    //   a[0] = this.node[2*P[0] + 0];//node[...];
    //   a[1] = this.node[2*P[0] + 1];//node[...];
    //   b[0] = this.node[2*P[1] + 0];//node[...];
    //   b[1] = this.node[2*P[1] + 1];//node[...];

    //   L = Math.sqrt( (a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]) );   

    //   console.log(`ie=${ie}`)
    //   this.getLT(a, b, L, T);

    //   // setLocalStiffness...
    //   this.setLocalStiffness(L, E, A, J, T, localStiffness);
    //   // transformation ...
    //   let elemu = Array(this.udim * this.edim);
    //   elemu[0] = result[this.udim*P[0] + 0];
    //   elemu[1] = result[this.udim*P[0] + 1];
    //   elemu[2] = result[this.udim*P[0] + 2];

    //   elemu[3] = result[this.udim*P[1] + 0];
    //   elemu[4] = result[this.udim*P[2] + 1];
    //   elemu[5] = result[this.udim*P[3] + 2];

    //   // elemu[...] = result[...];
    //   // multiplyAb...;
    //   const af     = new Array(this.kdim);
    //   this.multiplyAb(this.kdim, localStiffness, elemu, af);
    //   console.log(af);
    //   // aforces[0] = ...;
    //   // aforces[2] = ...;
    //   // aforces[4] = ...;
    //   // axialForcesPrint(edim, ndim, udim, node, Q, aforces, "f");
    // }


    // console.log(`\naxial forces by element`);

// for (let ie = 0; ie < this.esize; ie++) {
//     const E = this.eparam[3 * ie + 0];
//     const A = this.eparam[3 * ie + 1];
//     const J = this.eparam[3 * ie + 2];

//     const node1 = this.elem[2 * ie + 0]; // P[0]
//     const node2 = this.elem[2 * ie + 1]; // P[1]

//     // Координаты узлов
//     const x11 = this.node[2 * node1 + 0];
//     const y1 = this.node[2 * node1 + 1];
//     const x2 = this.node[2 * node2 + 0];
//     const y2 = this.node[2 * node2 + 1];

//     const L = Math.hypot(x2 - x11, y2 - y1);

//     // Матрица трансформации T (6x6 для плоского стержня с 3 степенями свободы на узел)
//     const T = new Array(this.kdim * this.kdim).fill(0);
//     this.getLT([x11, y1], [x2, y2], L, T);

//     // Глобальные перемещения узлов элемента (из решения result)
//     const dofPerNode = this.udim; // должно быть 3
//     const globalDofs = [
//         node1 * dofPerNode + 0, // u1
//         node1 * dofPerNode + 1, // v1
//         node1 * dofPerNode + 2, // theta1
//         node2 * dofPerNode + 0, // u2
//         node2 * dofPerNode + 1, // v2
//         node2 * dofPerNode + 2  // theta2
//     ];

//     const globalU = new Array(this.kdim);
//     for (let i = 0; i < this.kdim; i++) {
//         globalU[i] = result[globalDofs[i]];
//     }

//     // Перемещения в локальной системе координат: u_local = T * u_global
//     const localU = new Array(this.kdim);
//     for (let i = 0; i < this.kdim; i++) {
//         localU[i] = 0;
//         for (let j = 0; j < this.kdim; j++) {
//             localU[i] += T[i * this.kdim + j] * globalU[j];
//         }
//     }

//     // Локальная матрица жёсткости элемента (6x6) без трансформации
//     const localK = new Array(this.kdim * this.kdim).fill(0);
//     this.setLocalStiffness(L, E, A, J, T, localK);

//     // Локальные узловые силы: F_local = K_local * u_local
//     const localForces = new Array(this.kdim);
//     for (let i = 0; i < this.kdim; i++) {
//         localForces[i] = 0;
//         for (let j = 0; j < this.kdim; j++) {
//             localForces[i] += localK[i * this.kdim + j] * localU[j];
//         }
//     }

//     // Осевое усилие (N) — это горизонтальная сила в первом узле в локальной системе
//     // Обычно: localForces[0] — сила по оси X в узле 1 (осевая)
//     const axialForce = localForces[0];

//     const axialForceNode2 = localForces[3];

//     // console.log(`Element ${ie}: L = ${L.toFixed(4)}, Axial Force N = ${axialForce.toFixed(4)}`);

//     const forcesNode1 = [localForces[0], localForces[1], localForces[2]]; // N₁, Q₁, M₁
//     const forcesNode2 = [localForces[3], localForces[4], localForces[5]]; // N₂, Q₂, M₂

//     let x0 = Array(this.esize);
//     let x1 = Array(this.esize);
//     x0[0] = x11;
//     x0[1] = x2;
//     x1[0] = y1;
//     x1[1] = y2;
//     let f0 = Array(this.esize);
//     let f1 = Array(this.esize);
//     let f2 = Array(this.esize);
//     f0[0] = -forcesNode1[0];
//     f0[1] = forcesNode2[0];

//     f1[0] = forcesNode1[1];
//     f1[1] = -forcesNode2[1];

//     f2[0] = -forcesNode1[2];
//     f2[1] = forcesNode2[2];

//     const data = x0.map(
//         (item, P) => ({
//           'x0': item,
//           'x1': x1[P],
//           'f0': f0[P],
//           'f1': f1[P],
//           'f2': f2[P],
//         })
//       )
//       console.table(data)
    
//     // Если нужно сохранить в массив
//     // this.axialForces[ie] = axialForce;
// }


        

    return status;
  }
  
}

export { Isoparametric2d }
