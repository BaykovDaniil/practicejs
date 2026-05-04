import { axialForcesPrint, cg, logm, printResult, printSCR } from "./helper.js";



class Truss2d {
  static edim = 2; // количество узлов стержневого элемента
  static ndim = 2; // степень свободы узла на плоскости 2Д
  constructor() {
    this.edim = Truss2d.edim;
    this.ndim = Truss2d.ndim;
    this.elem = [
      2, 0,
      0, 1,
      1, 2
    ];

    this.node = [// [cm]
      0, 0,
      500, 0,
      900, 692.82032
    ];

    this.material = [// свойство стержня: area, E, alpha, theta 
      // area - площадь сечения [cm2], 
      // E - модуль Юнга [N/cm2], 
      // alpha - коэффициент теплового расширения alpha 
      // theta - перепад температуры
      10, 20e+6, 11e-6, 5, // ie = 0
      20, 20e+6, 11e-6, 5, // ie = 1
      20, 20e+6, 11e-6, 5  // ie = 2
    ];
    this.materialdim = 4;

    this.load = [// [N]
      0.0, 0.0,
      0.0, 0.0,
      0.0, -100.0e+3
    ];

    this.bound = [
      0.0, 0.0,
      undefined, 0.0,
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
  /**
  @param {Array<Double>} a - first node coordinate [x0, y0]
  @param {Array<Double>} b - second node coordinate [x1, y1]
  @param {Array<Double>} r - result element basis [e0, e1]
  @param {Array<Double>} g - result element metric [g00, g01, g10, g11]
  @returns {Double} length of elem
  */
  getLRG(a, b, r, g) {// do calculate L, r, g parameters
    const L = Math.sqrt( (b[0] - a[0])*(b[0] - a[0]) + (b[1] - a[1])*(b[1] - a[1]) );// длина конечного элемента(по координатам определяем) 
    r[0] = (b[0]-a[0]) / L; r[1] = (b[1]-a[1]) / L; // направляющие координаты
    g[0] = r[0] * r[0]; 
    g[1] = r[0] * r[1]; 
    g[2] = r[1] * r[0]; 
    g[3] = r[1] * r[1]; 
    return L;
  }

  /**
  @param {Integer} N - rows or cols of a square matrix
  @param {Array<Double>} u - kinematic vector
  @param {Array<Double>} f - result of load vector
  @param {Array<Double>} stiffness - result of stiffness
  */
  applyKinematic(N, u, f, stiffness) {// f, stiffness are modificated params
    for (let p = 0; p < N; p++) {
      for (let q = 0; q < N; q++) {
        if (u[p] === undefined)
          continue;
        if (p === q) {
          // f[p] = ;
        }
        else {
          // stiffness[p * N + q] = ;
          // f[q] -= ;
          // stiffness[q * N + p] =  ;
        }
      }
    }
    console.log("Kinematic boundary condition is applyed");
  }
  /**
  @param {Double} kx - coefficient stiffness EA/L
  @param {Array<Double>} gradG - graient of element матрица
  @param {Array<Double>} block - result of a stiffness block матрица 2 на 2 жесткости
  */
  blockStiffness(kx, gradG, block) {
    const edim = this.edim;
    for (let i = 0; i < block.length; i++) {
      block[i] = kx * gradG[i];
    }

    // for (let p = 0; p < edim; p++) {
    //   for (let q = 0; q < edim; q++) {
    //     // block[...] += ...;
    //   }
    // }
  }

  kinematicTest() { // Л.Сегерлинд. Применение метода конечных элементов стр. 110-112
    const u = [150, undefined, undefined, undefined, 40];
    const f = [500, 2000, 1000, 2000, 900];
    const stiffness = [
      55, -46, 4, 0, 0,
      -46, 140, -46, 0, 0,
      4, -46, 110, -46, 4,
      0, 0, -46, 142, -46,
      0, 0, 4, -46, 65,
    ];
    const N = u.length;
    logm('Before kinematic boundary condition', N, stiffness, f, 5, 0);
    this.applyKinematic(N, u, f, stiffness);
    logm('', N, stiffness, f, 5, 0);

    const sum =
      stiffness[0 * N + 0] - 55.0 + f[0] - 8250 +
      stiffness[1 * N + 1] - 140.0 + stiffness[1 * N + 2] + 46.0 + f[1] - 8900 +
      stiffness[2 * N + 1] + 46.0 + stiffness[2 * N + 2] - 110.0 + stiffness[2 * N + 3] + 46.0 + f[2] - 240 +
      stiffness[3 * N + 2] + 46.0 + stiffness[3 * N + 3] - 142.0 + f[3] - 3840 +
      stiffness[4 * N + 4] - 65.0 + f[4] - 2600;

    return sum == 0.0 ? 0 : 1;
  }

  // Преобразовать вектор-матрицу жёсткости в формат CSR https://docs.nvidia.com/cuda/cusparse/index.html#csr-format
  transformMatrixToCsr(M, N, stiffness) {
    const I = new Array(M + 1).fill(0);
    const J = [], val = [];
    let nz = 0, icounter = 0;

    for (let i = 0; i < M; i++) {
      // ...
      for (let j = 0; j < N; j++) {

        // ...
      }
    }
    // ...
    return [val, I, J]
  }

  transformMatrixTest() {
    const M = 4;
    const N = 5;

    const stiffness = [
      1.0, 4.0, 0.0, 0.0, 0.0,
      0.0, 2.0, 3.0, 0.0, 0.0,
      5.0, 0.0, 0.0, 7.0, 8.0,
      0.0, 0.0, 9.0, 0.0, 6.0,
    ];
    logm(`transformMatrixTest m = ${4}x${5}`, N, stiffness, [], 3, 0);
    const [val, I, J] = this.transformMatrixToCsr(M, N, stiffness);
    printSCR(val, I, J, 3, 0);
  }

  solve() {
    const N = this.nsize * this.ndim;
    const globalStiffness = new Array(N * N).fill(0);
    logm('globalStiffness matrix', N, globalStiffness, [], 3, 0);
    const P = [0, 0], Q = [0, 0], r = [0, 0], g = [0, 0, 0, 0],
      gradV = [-1, 1], gradG = [1, -1, -1, 1], a = [0, 0], b = [0, 0];
    const node = this.node, elem = this.elem, edim = this.edim, ndim = this.ndim,
      esize = this.esize;
    const elemStiffness = new Array(this.edim * this.edim).fill(0);
    const tmp = new Array(this.edim * this.edim).fill(0);
    const elemStiffness2d = new Array(edim * edim * ndim * ndim).fill(0);
    for (let ie = 0; ie < esize; ie++) {

      const A  = this.material[ie * this.materialdim + 0];//this.material[...];
      const E  = this.material[ie * this.materialdim + 1];//this.material[...];
      const al = this.material[ie * this.materialdim + 2];//this.material[...];
      const th = this.material[ie * this.materialdim + 3];//this.material[...];

      P[0] = this.elem[ie * this.edim + 0];
      P[1] = this.elem[ie * this.edim + 1];
      Q[0] = this.elem[ie * this.edim + 0];
      Q[1] = this.elem[ie * this.edim + 1];

      a[0] = this.node[P[0] * this.ndim + 0];
      a[1] = this.node[P[0] * this.ndim + 1];
      b[0] = this.node[P[1] * this.ndim + 0];
      b[1] = this.node[P[1] * this.ndim + 1];

      let L    = this.getLRG(a, b, r, g);
      const kx = E*A/L;
      this.blockStiffness(kx, gradG, elemStiffness);
      console.log(`ie=${ie}`);
      console.log(`L=${L} r=${r}`);
      logm('gij=', ndim, g);
      logm('Kpq=', edim, elemStiffness, [], 12, 1);

      tmp[0] = g[0] * elemStiffness[0];
      tmp[1] = g[1] * elemStiffness[1];
      tmp[2] = g[2] * elemStiffness[2];
      tmp[3] = g[3] * elemStiffness[3];

      elemStiffness2d.fill(0);

      elemStiffness2d[0]  = tmp[0];
      elemStiffness2d[1]  = -tmp[1];
      elemStiffness2d[4]  = -tmp[2];
      elemStiffness2d[5]  = tmp[3];

      elemStiffness2d[2]  = -tmp[0];
      elemStiffness2d[3]  = tmp[1];
      elemStiffness2d[6]  = tmp[2];
      elemStiffness2d[7]  = -tmp[3];

      elemStiffness2d[8]  = -tmp[0];
      elemStiffness2d[9]  = tmp[1];
      elemStiffness2d[12] = tmp[2];
      elemStiffness2d[13] = -tmp[3];

      elemStiffness2d[10] = tmp[0];
      elemStiffness2d[11] = -tmp[1];
      elemStiffness2d[14] = -tmp[2];
      elemStiffness2d[15] = tmp[3];

      for (let p = 0; p < edim; p++) {
        // add termal force
        // ...;
        // ...;
        

        for (let q = 0; q < edim; q++) {

          // globalStiffness[(P[p] + 0)*N*this.ndim + Q[q]*this.ndim + 0] += elemStiffness2d[0];
          // globalStiffness[(P[p] + 0)*N*this.ndim + Q[q]*this.ndim + 1] += elemStiffness2d[1];
          // globalStiffness[(P[p] + 0)*N*this.ndim+N + Q[q]*this.ndim + 0] += elemStiffness2d[2];
          // globalStiffness[(P[p] + 0)*N*this.ndim+N + Q[q]*this.ndim + 1] += elemStiffness2d[3];

          globalStiffness[(P[p] + 0)*N*this.ndim + Q[q]*this.ndim + 0] += tmp[0];
          globalStiffness[(P[p] + 0)*N*this.ndim + Q[q]*this.ndim + 1] += Math.pow(-1,(P[p] + 0) + Q[q] + 1)*tmp[1];
          globalStiffness[(P[p] + 0)*N*this.ndim+N + Q[q]*this.ndim + 0] += Math.pow(-1,(P[p] + 1) + Q[q] + 0)*tmp[2];
          globalStiffness[(P[p] + 0)*N*this.ndim+N + Q[q]*this.ndim + 1] += tmp[3];
          for (let i = 0; i < ndim; i++) {
            for (let j = 0; j < ndim; j++) {
              // elemStiffness2d[...] += ...;
              // elemStiffness2d[p*this.ndim + i] += g[] * elemStiffness[i*2 +j];
              // globalStiffness[(P[p] + 0)*N*this.ndim + Q[q]*this.ndim + j] += elemStiffness2d[i+j];
            }
          }
        }
      }
      logm('K^ij_pq=', edim * ndim, elemStiffness2d, [], 12, 1);
    }
    logm('Before kinematic boundary condition', N, globalStiffness, this.load, 12, 1);
    this.applyKinematic(N, this.bound, this.load, globalStiffness);
    logm('', N, globalStiffness, this.load, 12, 1);
    this.kinematicTest();
    this.transformMatrixTest();

    const [val, I, J] = this.transformMatrixToCsr(N, N, globalStiffness);
    printSCR(val, I, J, 12, 1);

    const result = new Array(N).fill(0);
    const { succesfull, error } = cg(N, I, J, val, this.load, 1, result);
    
    // console.log(`Расчет выполнен ${succesfull ? `успешно` : `с ошибками`}! `, `Наибольшая невязка: ${error}`);
    printResult(ndim, this.node, result, 14, 5);

    const aforces = new Array(esize).fill(0);

    for (let ie = 0; ie < esize; ie++) {
      const A = 1.0;//this.material[...];
      const E = 1.0;//this.material[...];
      const al = 1.0;//this.material[...];
      const th = 1.0;//this.material[...];

      // P[0] = elem[...];
      // P[1] = ...;

      // Q[0] = ...;
      // Q[1] = ...;

      // a[0] = node[...];
      // a[1] = ...;
      // b[0] = ...;
      // b[1] = ...;

      const L = this.getLRG(a, b, r, g);

      // const dN = ...;
      // aforces[ie] -= ...;
      for (let q = 0; q < edim; q++) {
        for (let j = 0; j < ndim; j++) {
          // aforces[ie] += ...;
        }
      }
    }

    axialForcesPrint(aforces);

  }
}
export { Truss2d }