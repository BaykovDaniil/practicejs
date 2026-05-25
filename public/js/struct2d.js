import { Vector3 } from "three";
import { axialForcesPrint, cg, logm, printResult, printSCR } from "./helper.js";
// import { sqrt } from "three/tsl";

class Struct2d {

  static edim = 2; // количество узлов стержневого элемента
  static ndim = 2; // степень свободы узла на плоскости 2Д
  static udim = 3;
  static kdim = this.udim * this.edim;
  constructor(options = {}) {
    this.edim = Struct2d.edim;
    this.ndim = Struct2d.ndim;
    this.udim = Struct2d.udim;
    this.kdim = Struct2d.kdim;

    this.pdim = 3;
    this.elem = [
      0, 2,
      2, 3,
      3, 4,
      4, 5,
      5, 1,
    ];

    this.node = [// [cm]
      0, 0,
      900, 0,
      0, 600,
      300, 600,
      600, 600,
      900, 600,
    ];

    this.elastic = [20.6, 0.3]; // { [MN/cm2], [1]}: модуль Юнга и коэфф. Пуассона
    this.thikness = 1.0; // [cm]: толщина
    this.alpha = 11.3e-06;
    this.thetta = 0;
    if (options.thetta !== undefined) {
      this.thetta = options.thetta;
    }
    this.eparam = [ // E|N/cm2|, A|cm2|, J|cm4| by element
      20e+6, 30, 400,
      20e+6, 30, 400,
      20e+6, 30, 400,
      20e+6, 30, 400,
      20e+6, 30, 400,
    ];

    this.load = [ // [ux,uy,rz]
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, -100000, 25000,
        0.0, -100000, 25000,
        0.0, 0.0, 0.0,
    ];

    this.ieload = { // [py,mz]
        0: {'q_0': 200000/600, 'm_0': 0, 'x':[0, 0.5, 1]}
        // 4: {'q_0': 200000/600, 'm_0': 0, 'x':[0, 0.5, 1]}
    };

    // { // [py,mz]
    //     1: {'q_0': 10, 'm_0': 10, 'x':[0, 0.5, 1]},
    //     2: {'q_0': 10, 'm_0': 10, 'x':[0, 0.5, 1]},
    //     3: {'q_0': 10, 'm_0': 10, 'x':[0, 0.5, 1]}
    // };

    this.bound = [ // [cm]: Кинематические граничные условия
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      undefined, undefined, undefined,
      undefined, undefined, undefined,
      undefined, undefined, undefined,
      undefined, undefined, undefined,
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

  getElementLoad(ie) {
    const load = this.ieload && this.ieload.hasOwnProperty(ie) ? this.ieload[ie] : {};

    if (Array.isArray(load)) {
      return {
        q: -(load[0] || 0),
        m: load[1] || 0,
        x: load[2] || [0, 0.5, 1],
      };
    }

    return {
      q: -(load['q_0'] || 0),
      m: load['m_0'] || 0,
      x: load['x'] || [0, 0.5, 1],
    };
  }

  getEquivalentLocalLoad(L, q, m) {
    const feLocal = new Array(this.kdim).fill(0);
    const L2 = L * L;

    feLocal[1] += q * L / 2.0;
    feLocal[4] += q * L / 2.0;

    feLocal[2] += q * L2 / 12.0;
    feLocal[5] -= q * L2 / 12.0;

    feLocal[1] -= m;
    feLocal[4] += m;

    return feLocal;
  }

  getEquivalentThermalLoad(E, A) {
    const feLocal = new Array(this.kdim).fill(0);
    const dN = E * A * this.alpha * this.thetta;

    feLocal[0] = -dN;
    feLocal[3] = dN;

    return feLocal;
  }

  addLocalLoadToGlobal(T, feLocal, globalDofs, load) {
    const feGlobal = new Array(this.kdim).fill(0);
    for (let i = 0; i < this.kdim; i++) {
      for (let j = 0; j < this.kdim; j++) {
        feGlobal[i] += T[j * this.kdim + i] * feLocal[j];
      }
    }

    for (let i = 0; i < this.kdim; i++) {
      const gdof = globalDofs[i];
      load[gdof] += feGlobal[i];
    }
  }

  calcM(x, endForces, q, m) {
    return -endForces[2] + (endForces[1] - m) * x + q * x*x / 2.0;
  }

  calcQ(x, endForces, q, m) {
    return endForces[1] - 2*m + q * x;
  }

  calcV(x, L, localU, normalSign) {
    const xi = x / L;
    const xi2 = xi * xi;
    const xi3 = xi2 * xi;

    const N1 = 2*xi3 - 3*xi2 + 1;
    const N2 = L * (xi3 - 2*xi2 + xi);
    const N3 = 3*xi2 - 2*xi3;
    const N4 = L * (xi3 - xi2);

    const v = N1*localU[1] + N2*localU[2] + N3*localU[4] + N4*localU[5];
    return normalSign * v;
  }

  getVNormalSign(x1, y1, x2, y2) {
    const L = Math.hypot(x2 - x1, y2 - y1);
    const jx = -(y2 - y1) / L;
    const jy = (x2 - x1) / L;
    const normal = Math.abs(jx) > Math.abs(jy) ? jx : jy;

    return normal < 0 ? -1 : 1;
  }

  solve() {

    // const v = this.elastic[1];
    // const E = this.elastic[0];

    // const lambda = v*E / ( (1+v) * (1-2*v) );
    // const mu = E / ( 2*(1+v));

    const N = this.nsize * this.udim;

    const globalStiffness = new Array(N * N);
    for (let i = 0; i < N * N; i++) {
        globalStiffness[i] = 0;
    }

    const load = this.load.slice();

    // this.printStiffness(N, globalStiffness);
    logm('globalStiffness matrix', N, globalStiffness, [], 3, 0);


    // const a = [0, 0, 0], b = [0, 0, 0], c = [0, 0, 0];
    const P = [0, 0], Q = [0, 0];
    const T = new Array(this.kdim * this.kdim);
    const localStiffness = new Array(this.kdim * this.kdim);
    for (let i = 0; i < this.kdim * this.kdim; i++) {
        T[i] = 0;
        localStiffness[i] = 0;
    }

    const a = [0, 0], b = [0, 0];
    const elemStiffness = new Array(this.kdim * this.kdim);
    

    for (let ie = 0; ie < this.esize; ie++) {
      const E = this.eparam[3*ie + 0];//eparam[...];
      const A = this.eparam[3*ie + 1];//eparam[...];
      const J = this.eparam[3*ie + 2];//eparam[...];

      P[0] = this.elem[2*ie + 0];//elem[...];
      P[1] = this.elem[2*ie + 1];//elem[...];

      Q[0] = this.elem[2*ie + 0]; //elem[...];
      Q[1] = this.elem[2*ie + 1]; //elem[...]; 

      let L = 1;
      a[0] = this.node[2*P[0] + 0];//node[...];
      a[1] = this.node[2*P[0] + 1];//node[...];
      b[0] = this.node[2*P[1] + 0];//node[...];
      b[1] = this.node[2*P[1] + 1];//node[...];
      console.log(a);
      console.log(b);
      L = Math.sqrt( (a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]) );

      this.getLT(a, b, L, T);

      console.log(`elem ${ie}`);     
      console.log(`L = ${L}`);
      logm('Matrix Transformation', this.edim * this.udim, T, 5, 0);
  
      this.setLocalStiffness(L, E, A, J, T, localStiffness);
      logm('Local Stiffness [K] (in local basis)', this.kdim, localStiffness, 5, 0);

      const result_tmp = new Array(this.kdim * this.kdim);
      const result     = new Array(this.kdim * this.kdim);
      this.multiplyAB(this.kdim, T, localStiffness, result_tmp, 1);
      this.multiplyAB(this.kdim, result_tmp, T, result, 0);
      logm('localStiffness [T]*[K][T] (in global basis)', this.kdim, result, 5, 0);


      const dofPerNode = this.udim;
      const node1 = P[0];
      const node2 = P[1];


      const globalDofs = [
          node1 * dofPerNode + 0,  // u в узле 1
          node1 * dofPerNode + 1,  // v в узле 1
          node1 * dofPerNode + 2,  // z в узле 1
          node2 * dofPerNode + 0,  // u в узле 2
          node2 * dofPerNode + 1,   // v в узле 2
          node2 * dofPerNode + 2   // z в узле 2
      ];

    
    
      if (this.ieload && this.ieload.hasOwnProperty(ie)) {
        const { q, m } = this.getElementLoad(ie);
        const feLocal = this.getEquivalentLocalLoad(L, q, m);
        this.addLocalLoadToGlobal(T, feLocal, globalDofs, load);
      }

      const feThermal = this.getEquivalentThermalLoad(E, A);
      this.addLocalLoadToGlobal(T, feThermal, globalDofs, load);




      for (let i = 0; i < this.kdim; i++) {
          for (let j = 0; j < this.kdim; j++) {
              const globalRow = globalDofs[i];
              const globalCol = globalDofs[j];
              const index = globalRow * N + globalCol;
              globalStiffness[index] += result[i * this.kdim + j];
          }
      }
    }

    this.load = load;

    logm('globalStiffness matrix', N, globalStiffness, this.load, 8, 3);

    let status = this.sumStiffnessTest(N, globalStiffness);

    if (status !== 0) {
        console.log("\nSum stiffness error\n");
        return status;
    } else {
        console.log("\nSum stiffness ok\n");
    }

    
    logm('globalStiffness matrix', N, globalStiffness, this.load, 8, 3);

    status = this.kinematicTest();

    if (status !== 0) {
        console.log("\nKinematic error\n");
        return status;
    } else {
        console.log("\nKinematic test ok\n");
    }

    this.applyKinematic(N, this.bound, this.load, globalStiffness);

    console.log(`Kinematic boundary condition is applied`); 
    logm('global stiffness K with kinematic', N, globalStiffness, this.load, 8, 3);

    // this.transformMatrixTest();

    const val = [];
    const I = [];
    const J = [];
    let nz;
    this.transformMatrixToCsr(N, N, nz, globalStiffness, val, I, J);

    const result = new Array(N).fill(0);
    const { succesfull, error } = cg(N, I, J, val, this.load, 1, result);

    this.resultPrint(result, N);

    console.log(`axial forces by element`);

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
    const locals = {};

    for (let ie = 0; ie < this.esize; ie++) {
      console.log(`ie=${ie}`);
      const E = this.eparam[3 * ie + 0];
      const A = this.eparam[3 * ie + 1];
      const J = this.eparam[3 * ie + 2];

      const node1 = this.elem[2 * ie + 0]; // P[0]
      const node2 = this.elem[2 * ie + 1]; // P[1]

      // Координаты узлов
      const x11 = this.node[2 * node1 + 0];
      const y1 = this.node[2 * node1 + 1];
      const x2 = this.node[2 * node2 + 0];
      const y2 = this.node[2 * node2 + 1];

      const L = Math.hypot(x2 - x11, y2 - y1);

      // Матрица трансформации T
      const T = new Array(this.kdim * this.kdim).fill(0);
      this.getLT([x11, y1], [x2, y2], L, T);

      // Глобальные перемещения узлов элемента
      const dofPerNode = this.udim; // должно быть 3
      const globalDofs = [
          node1 * dofPerNode + 0, // u1
          node1 * dofPerNode + 1, // v1
          node1 * dofPerNode + 2, // theta1
          node2 * dofPerNode + 0, // u2
          node2 * dofPerNode + 1, // v2
          node2 * dofPerNode + 2  // theta2
      ];

      const globalU = new Array(this.kdim);
      for (let i = 0; i < this.kdim; i++) {
          globalU[i] = result[globalDofs[i]];
      }

      // u_local = T * u_global
      const localU = new Array(this.kdim);
      for (let i = 0; i < this.kdim; i++) {
          localU[i] = 0;
          for (let j = 0; j < this.kdim; j++) {
              localU[i] += T[i * this.kdim + j] * globalU[j];
          }
      }

      console.log(localU);

      // Локальная матрица жёсткости(6x6) без трансформации
      const localK = new Array(this.kdim * this.kdim).fill(0);
      this.setLocalStiffness(L, E, A, J, T, localK);

      // F_local = K_local * u_local
      const localForces = new Array(this.kdim);
      for (let i = 0; i < this.kdim; i++) {
          localForces[i] = 0;
          for (let j = 0; j < this.kdim; j++) {
              localForces[i] += localK[i * this.kdim + j] * localU[j];
          }
      }

      const { q, m } = this.getElementLoad(ie);
      const feLocal = this.getEquivalentLocalLoad(L, q, m);
      const feThermal = this.getEquivalentThermalLoad(E, A);
      const endForces = localForces.map((force, index) => force - feLocal[index] - feThermal[index]);
      locals[ie] = { localU, localForces, feLocal, feThermal, endForces, q, m, L };

      // Осевое усилие (N) — это горизонтальная сила в первом узле в локальной системе
      const axialForce = localForces[0];

      const axialForceNode2 = localForces[3];

      // console.log(`Element ${ie}: L = ${L.toFixed(4)}, Axial Force N = ${axialForce.toFixed(4)}`);

      const forcesNode1 = [endForces[0], endForces[1], endForces[2]]; // N₁, Q₁, M₁
      const forcesNode2 = [endForces[3], endForces[4], endForces[5]]; // N₂, Q₂, M₂

      let x0 = Array(this.esize);
      let x1 = Array(this.esize);
      x0[0] = x11;
      x0[1] = x2;
      x1[0] = y1;
      x1[1] = y2;
      let f0 = Array(this.esize);
      let f1 = Array(this.esize);
      let f2 = Array(this.esize);
      f0[0] = -forcesNode1[0];
      // / this.eparam[1];
      f0[1] = forcesNode2[0];
      // / this.eparam[1];

      f1[0] = forcesNode1[1] - 2*m;
      f1[1] = -forcesNode2[1] - 2*m;

      f2[0] = -forcesNode1[2];
      f2[1] = forcesNode2[2];

      const data = x0.map(
          (item, P) => ({
            'x0': item,
            'x1': x1[P],
            'f0': f0[P],
            'f1': f1[P],
            'f2': f2[P],
          })
        )
        console.table(data)
      
      // Если нужно сохранить в массив
      // this.axialForces[ie] = axialForce;
    }


    // M и Q
    console.log('M and Q at several points');
    for (let ie = 0; ie < this.esize; ie++) {

      console.log(`ie=${ie}`);
      
      const elementState = locals[ie];
      const { q, m, L, endForces } = elementState;
      const x = this.getElementLoad(ie).x;
      const M = new Array(x.length);
      const shearForce = new Array(x.length);
      const deflection = new Array(x.length);

      const i1 = this.elem[2*ie + 0];
      const i2 = this.elem[2*ie + 1];

      const x11 = this.node[this.ndim*i1 + 0];
      const x12 = this.node[this.ndim*i1 + 1];
      const x21 = this.node[this.ndim*i2 + 0];
      const x22 = this.node[this.ndim*i2 + 1];
      const VNormalSign = this.getVNormalSign(x11, x12, x21, x22);
      // u0[i] = result[this.udim*i + 0];
      // const u11 = result[this.udim*i1 + 1];
      // const u12 = result[this.udim*i1 + 2];
      // const u21 = result[this.udim*i2 + 1];
      // const u22 = result[this.udim*i2 + 2];

      const E = this.eparam[3*ie + 0];
      const J = this.eparam[3*ie + 2];

      console.log(L, E, J, q, m);

      console.log(elementState.localU);
      for(let i=0; i<x.length; i++) {
        const xi = x[i] * L;
        M[i] = this.calcM(xi, endForces, q, m);
        shearForce[i] = this.calcQ(xi, endForces, q, m);
        deflection[i] = this.calcV(xi, L, elementState.localU, VNormalSign);
        
      }
      
      const data = x.map(
        (item, index) => ({
          'x/L': item,
          'M': M[index],
          'Q': shearForce[index],
          'v': deflection[index],
        })
      )
      console.table(data)

      
      
    }

    return status;
  }
  
}

export { Struct2d }
