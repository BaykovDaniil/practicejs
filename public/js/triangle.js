import { Vector3 } from "three";
import { axialForcesPrint, cg, logm, printResult, printSCR } from "./helper.js";
// import { sqrt } from "three/tsl";

class Triangle {

  static edim = 3; // количество узлов стержневого элемента
  static ndim = 2; // степень свободы узла на плоскости 2Д
  constructor() {
    this.edim = Triangle.edim;
    this.ndim = Triangle.ndim;
    this.elem = [
      0, 1, 3,
      1, 4, 3,
      0, 3, 2,
    ];

    this.node = [// [cm]
      -4.0, 0.0,
      -2.0, 0.0,
      -5.0, 1.7,
      -3.0, 1.7,
      -1.0, 1.7,
    ];

    this.elastic = [20.6, 0.3]; // { [MN/cm2], [1]}: модуль Юнга и коэфф. Пуассона
    this.thikness = 1.0; // [cm]: толщина
    this.alpha = 11.3e-06;
    this.thetta = 100;

    this.load = [ // [MN]: узловые силы
      0.0, 0.0,
      0.0, 0.0,
      0.0, -3.0e-3,
      0.0, -3.0e-3,
      0.0, 0.0,
    ];

    this.bound = [ // [cm]: Кинематические граничные условия
      0.0, 0.0,
      0.0, 0.0,
      undefined, undefined,
      undefined, undefined,
      undefined, undefined
    ];

    this.solve();
  }

  get esize() {
    return this.elem.length / this.edim;
  }

  get nsize() {
    return this.node.length / this.ndim;
  }

  gradientA(P, result) {
    // const thiss = {node:[1,2,3], ndim: 9};
    // const {node, ndim} = this;
    result[0] = this.node[P[1] * this.ndim + 0] * this.node[P[2] * this.ndim + 1] 
                -this.node[P[2] * this.ndim + 0] * this.node[P[1] * this.ndim + 1];

    result[1] = this.node[P[2] * this.ndim + 0] * this.node[P[0] * this.ndim + 1] 
                -this.node[P[0] * this.ndim + 0] * this.node[P[2] * this.ndim + 1];

    result[2] = this.node[P[0] * this.ndim + 0] * this.node[P[1] * this.ndim + 1] 
                -this.node[P[1] * this.ndim + 0] * this.node[P[0] * this.ndim + 1];
    //   result[1] = ;
    //   result[2] = ;

  }

  gradientB(P, result) {

    result[0] = this.node[P[1] * this.ndim + 1]
                -this.node[P[2] * this.ndim + 1];

    result[1] = this.node[P[2] * this.ndim + 1]
                -this.node[P[0] * this.ndim + 1];

    result[2] = this.node[P[0] * this.ndim + 1] 
                -this.node[P[1] * this.ndim + 1];

  }

  gradientC(P, result) {

    result[0] = this.node[P[2] * this.ndim + 0]
                -this.node[P[1] * this.ndim + 0];

    result[1] = this.node[P[0] * this.ndim + 0]
                -this.node[P[2] * this.ndim + 0];

    result[2] = this.node[P[1] * this.ndim + 0] 
                -this.node[P[0] * this.ndim + 0];

  }


  getDet(a) {

    return a[0] + a[1] + a[2];

  }

  blockStressStiffnessPS(det, bp, bq, cp, cq, block) {

    const v = this.elastic[1];
    const E = this.elastic[0];

    const a = det / 2 * E / (det*det*(1-v*v));

    block[0 * this.ndim + 0] = a * (bp*bq + (1-v)/2 * cp*cq);
    block[0 * this.ndim + 1] = a * (bp*cq*v + (1-v)/2 * cp*bq);

    block[1 * this.ndim + 0] = a * (cp*bq*v + (1-v)/2 * bp*cq);
    block[1 * this.ndim + 1] = a * (cp*cq + (1-v)/2 * bp*bq);

  }

  blockStressStiffnessPE(det, bp, bq, cp, cq, block) {

    const v = this.elastic[1];
    const E = this.elastic[0];

    const lambda = v*E / ( (1+v) * (1-2*v) );
    const mu = E / ( 2*(1+v));

    const a = det* 1 / (det*det) / 2;

    block[0 * this.ndim + 0] = a * ((lambda + 2*mu)*bp*bq + mu*cp*cq);
    block[0 * this.ndim + 1] = a * (lambda*bp*cq + mu*cp*bq);

    block[1 * this.ndim + 0] = a * (lambda*cp*bq + mu*bp*cq);
    block[1 * this.ndim + 1] = a * ((lambda + 2*mu)*cp*cq + mu*bp*bq);

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

  resultPrint(nsize, ndim, node, result) {

    let x0 = Array(nsize);
    let x1 = Array(nsize);
    let x2 = Array(nsize);

    let u0 = Array(nsize);
    let u1 = Array(nsize);
    let u2 = Array(nsize);

    if( ndim === 2) {
      for(let i = 0; i < nsize; i++) {
        x0[i] = node[ndim*i];
        x1[i] = node[ndim*i + 1];
        u0[i] = result[ndim*i];
        u1[i] = result[ndim*i + 1];
      }
      const data = x0.map(
        (item, index) => ({
          'x0': item,
          'x1': x1[index],
          'u0': u0[index],
          'u1': u1[index],
        })
      )
      console.table(data)
    }

    if( ndim === 3) {
      for(let i = 0; i < nsize; i++) {
        x0[i] = node[ndim*i];
        x1[i] = node[ndim*i + 1];
        x2[i] = node[ndim*i + 2];
        u0[i] = result[ndim*i];
        u1[i] = result[ndim*i + 1];
        u2[i] = result[ndim*i + 2];
      }
      const data = x0.map(
        (item, index) => ({
          'x0': item,
          'x1': x1[index],
          'x2': x1[index],
          'u0': u0[index],
          'u1': u1[index],
          'u2': u1[index],
        })
      )
      console.table(data)
    }

  }

  // CSR_dot_vector(val, I, J, vect) { // TODO

  //   const result = new Array(vect.length);

  //   // for

  // }

  // vector_norm_l2(vector) {

  //   let sum = 0;

  //   for(let i = 0; i < vector.length; i++) {
  //     sum += vector[i] ** 2;
  //   }

  //   return sqrt(sum);
  // }

  // scalar_product(v1, v2) {

  //   let result = 0;

  //   for(let i = 0; i < v1.length; i++) {
  //     result += v1[i] * v2[i];
  //   }

  //   return result;

  // }

  // cg_host(N, nz, I, J, val, load, result, t){

  //   const max_iter = 100;
  //   const status = 1;

  //   const x_old = new Array(N);
  //   const x_new = new Array(N);
  //   for(let i = 0; i < N; i++) {
  //     x_old[i] = 0;
  //   }
    
  //   const r_old = [...load];
  //   const r_new = new Array(N);
  //   const z_old = [...r_old];
  //   const z_new = new Array(N);
  //   const Ax = [...x_old];

  //   for(let k = 0; k < max_iter; k++) {

  //     const Az_old = this.CSR_dot_vector(val, I, J, z_old);
  //     alpha = this.scalar_product(r_old, r_old) / this.scalar_product(Az_old, z_old);

  //     for(let i = 0; i < N; i++) {
  //       x_new[i] = x_old + alpha * z_old[i];
  //     }

  //     for(let i = 0; i < N; i++) {
  //       r_new[i] = r_old - alpha * Az_old[i];
  //     }

  //     betta = this.scalar_product(r_new, r_new) / this.scalar_product(r_old, r_old);

  //     for(let i = 0; i < N; i++) {
  //       z_new[i] = r_new[i] + betta * z_old[i];
  //     }

  //     if(this.vector_norm_l2(r_new) / this.vector_norm_l2(load) < t) {
  //       status = 0;
  //       break;
  //     }


  //   }
    
  //   for(let i = 0; i < N; i++) {
  //     result[i] = x_new[i];
  //   }

  //   return status;

  // }

  // printStiffness(N, globalStiffness){
  //   for (let i = 0; i < N * N; i++) {
  //     console.log(globalStiffness[i]);
  //   }
  // }

  solve() {

    const v = this.elastic[1];
    const E = this.elastic[0];

    const lambda = v*E / ( (1+v) * (1-2*v) );
    const mu = E / ( 2*(1+v));

    const N = this.nsize * this.ndim;

    const globalStiffness = new Array(N * N);
    for (let i = 0; i < N * N; i++) {
        globalStiffness[i] = 0;
    }

    // this.printStiffness(N, globalStiffness);
    logm('globalStiffness matrix', N, globalStiffness, [], 3, 0);


    const a = [0, 0, 0], b = [0, 0, 0], c = [0, 0, 0];
    const P = [0, 0, 0], Q = [0, 0, 0];
    const block = new Array(this.ndim * this.ndim);

    this.load_temperature = new Array(this.nsize * this.ndim);
    for (let i = 0; i < this.nsize * this.ndim; i++) {
        this.load_temperature[i] = this.load[i];
    }
    console.log(this.load_temperature);

    for (let ie = 0; ie < this.esize; ie++) {

        P[0] = this.elem[ie * this.edim + 0];
        P[1] = this.elem[ie * this.edim + 1];
        P[2] = this.elem[ie * this.edim + 2];

        Q[0] = this.elem[ie * this.edim + 0];
        Q[1] = this.elem[ie * this.edim + 1];
        Q[2] = this.elem[ie * this.edim + 2];

        this.gradientA(P, a);
        this.gradientB(P, b);
        this.gradientC(P, c);
        let det = this.getDet(a);

        console.log(`elem ${ie}`);
        console.log(`det ${det}`);
        const data_abc = a.map(
          (item, index) => ({
            'a': item,
            'b': b[index],
            'c': c[index],
          })
        )
        console.table(data_abc)

        // const elemStiffness = new Array(this.ndim * this.edim * this.ndim * this.edim);

        for (let p = 0; p < this.edim; p++) { // кэ

          const k = this.alpha * this.thikness * this.thetta * (3*lambda + 2*mu) / 2

          this.load_temperature[Q[p]*this.ndim + 0] += b[p] * k;
          this.load_temperature[Q[p]*this.ndim + 1] += c[p] * k;

            for (let q = 0; q < this.edim; q++) {

              // координатный формат

              // симметричность
              // посмтотреть какие элементы
              // P b Q не соприкасаются, то нули

                this.blockStressStiffnessPE(det, b[p], b[q], c[p], c[q], block);

                console.log(`\nPQ[${P[p]}${Q[q]}] pq[${p}${q}]\n`);


                // print stiffness
                console.log(block[0 * this.ndim + 0], block[0 * this.ndim + 1]);
                console.log(block[1 * this.ndim + 0], block[1 * this.ndim + 1]);

                globalStiffness[(P[p] + 0)*N*this.ndim + Q[q]*this.ndim + 0] += block[0];
                globalStiffness[(P[p] + 0)*N*this.ndim + Q[q]*this.ndim + 1] += block[1];
                globalStiffness[(P[p] + 0)*N*this.ndim+N + Q[q]*this.ndim + 0] += block[2];
                globalStiffness[(P[p] + 0)*N*this.ndim+N + Q[q]*this.ndim + 1] += block[3];

                for (let i = 0; i < this.ndim; i++) {
                    for (let j = 0; j < this.ndim; j++) {
                        // globalStiffness[...] += block[...];
                        // globalStiffness[P[p]*this.nsize Q[q]*this.nsize] += block[i * this.ndim + j];
                    }
                }
            }
        }

        // this.printStiffness(this.ndim * this.edim, elemStiffness);
    }

    logm('globalStiffness matrix', N, globalStiffness, this.load, 8, 3);

    let status = this.sumStiffnessTest(N, globalStiffness); // TODO

    if (status !== 0) {
        console.log("\nSum stiffness error\n");
        return status;
    } else {
        console.log("\nSum stiffness ok\n");
    }

    // this.printStiffnessF(N, globalStiffness, this.load);
    let globalStiffness_temperature = new Array(N * N);
    for (let i = 0; i < N * N; i++) {
        globalStiffness_temperature[i] = globalStiffness[i];
    }
    this.applyKinematic(N, this.bound, this.load, globalStiffness);
    this.applyKinematic(N, this.bound, this.load_temperature, globalStiffness_temperature);
    logm('globalStiffness matrix', N, globalStiffness, this.load, 8, 3);

    status = this.kinematicTest();

    if (status !== 0) {
        console.log("\nKinematic error\n");
        return status;
    } else {
        console.log("\nKinematic test ok\n");
    }

    this.transformMatrixTest();

    const val = [];
    const I = [];
    const J = [];
    let nz;

    const val_temperature = [];
    const I_temperature = [];
    const J_temperature = [];
    let nz_temperature;

    this.transformMatrixToCsr(N, N, nz, globalStiffness, val, I, J);
    this.transformMatrixToCsr(N, N, nz_temperature, globalStiffness_temperature, val_temperature, I_temperature, J_temperature);
    this.printSCR_Matrix(N, N, nz, val, I, J);

    const result = new Array(N).fill(0);
    const result_temperature = new Array(N).fill(0);

    const { succesfull, error } = cg(N, I, J, val, this.load, 1, result);
    const { succesfull_temperature, error_temperature } = cg(N, I_temperature, J_temperature, val_temperature, this.load_temperature, 1, result_temperature);

    // status = this.cg_host(N, nz, I, J, val, this.load, result, 1);
    console.log("Без учета теплового расширения");
    printResult(this.ndim, this.node, result);
    // this.resultPrint(this.nsize, this.ndim, this.node, result);

    console.log("С учетом теплового расширения");
    printResult(this.ndim, this.node, result_temperature);
    // this.resultPrint(this.nsize, this.ndim, this.node, result_temperature);
    // logm("Ux, Uy", 2, result, []);

    // Cochi equation
    // Constutive equation
    const epsilon = [0, 0, 0];
    const sigma = [0, 0, 0];

    // const v = this.elastic[1];
    // const E = this.elastic[0];
    const la = E*v / ( (1+v)*(1-2*v) );
    // const mu = 
    // const ka = 

    // this.printTensorH();
    for (let ie = 0; ie < this.esize; ie++) {
        P[0] = this.elem[ie * this.edim + 0];
        P[1] = this.elem[ie * this.edim + 1];
        P[2] = this.elem[ie * this.edim + 2];

        this.gradientA(P, a);
        this.gradientB(P, b);
        this.gradientC(P, c);
        let det = this.getDet(a);
        const epsilon = [0, 0, 0];
        const sigma = [0, 0, 0];
        for (let p = 0; p < this.edim; p++) {
            epsilon[0] += b[p] / det * result_temperature[P[p]*this.ndim + 0]; // p-1
            epsilon[1] += c[p] / det * result_temperature[P[p]*this.ndim + 1]; // p-2
            epsilon[2] += (c[p] / det * result_temperature[P[p]*this.ndim + 0] + b[p] / det * result_temperature[P[p]*this.ndim + 1]); // p-1  p-2
            // plane strain
            // sigma[0] += lambda * (epsilon[0]+epsilon[1]) + 2*mu*epsilon[0];
            // sigma[1] += lambda * (epsilon[0]+epsilon[1]) + 2*mu*epsilon[1];
            // sigma[2] += mu*epsilon[2];

            // plane stress
            // sigma[0] += 
            // sigma[1] += 
            // sigma[2] += 
        }
        sigma[0] += la * (epsilon[0]+ epsilon[1]) + 2*mu*(epsilon[0]) - (3*la + 2*mu)* this.alpha*this.thetta;
        sigma[1] += la * (epsilon[0] + epsilon[1]) + 2*mu*(epsilon[1]) - (3*la + 2*mu)* this.alpha*this.thetta;
        sigma[2] += mu*epsilon[2];
        // console.log(epsilon);
        console.log(sigma);
        // this.printTensor(ie, epsilon, sigma);
      }

    return status;
  }
  
}

export { Triangle }