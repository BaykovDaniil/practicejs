import { Vector3 } from "three";
import { axialForcesPrint, cg, logm, printResult, printSCR } from "./helper.js";
// import { sqrt } from "three/tsl";
// COO -> CSR на следующую пару сделать
// на следующей парепосчитаем перемещения, напряжения, деформации с тепловыми расширениями и морфирование, поверхностные силы(численно интегрировать, [ie, грань, P])

class Isoparametric2d {

  static edim = 4; // количество узлов элемента
  static ndim = 2; // степень свободы узла на плоскости 2Д
  static udim = 2;
  static kdim = this.udim * this.edim;
  constructor(flag = "ПД") {
    this.edim = Isoparametric2d.edim;
    this.ndim = Isoparametric2d.ndim;
    this.udim = Isoparametric2d.udim;
    this.kdim = Isoparametric2d.kdim;

    this.pdim = 3;
    this.elem = [
      0, 1, 2, 3
    ];

    this.node = [// [cm]
      -0.5, -0.5,
      0.5, -0.5,
      0.5, 0.5,
      -0.5, 0.5
    ];

    this.elastic = [20.6, 0.3]; // { [MN/cm2], [1]}: модуль Юнга и коэфф. Пуассона
    this.flag = flag; // "ПД" — плоская деформация, "ПН" — плоское напряжение
    this.thikness = 1.0; // [cm]: толщина
    this.alpha = 11.3e-06;
    this.thetta = 100;
    // this.eparam = [ // E|N/cm2|, A|cm2|, J|cm4| by element
    //   20e+6, 30, 400,
    // ];

    this.load = [ // [ux,uy]
        0.0, 0.0,
        0.0, 0.0,
        0.0, -1e-3,
        0.0, -1e-3
    ];

    // this.ieload = { // [py,mz]
    //     0: [10, 10],
    //     1: [10, 10]
    // };

    // Только нулевые!!!!
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

  printSCR_Matrix(M, N, nz, val, I, J) {

    const matrix = new Array(N * M);
    for (let i = 0; i < M * N; i++) {
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

    for(let P = 0; P < this.nsize; P++) {
      result[P] = [];
    }

    for(let ie = 0; ie < this.esize; ie++) {
      for(let p = 0; p < this.edim; p++) {

        const P = this.elem[ie * this.edim + p];

        for(let q = 0; q < this.edim; q++) {

          const Q = this.elem[ie * this.edim + q];

          if(Q < P) {
            continue;
          }


          // вынести
          let i = 0;
          while(i < result[P].length && result[P][i] < Q) {
            i++;
          }

          if(result[P][i] !== Q) {
            result[P].splice(i, 0, Q);
          }

          // if(result[P].indexOf(Q) === -1) {
          //   // result[P].push(Q);
          //   if(Q==P) {
          //     result[P].push(Q);
          //   }
          //   for(const i = 0; i < result[P].length; i++) {
          //     if(result[P][i] < Q) {
          //       result[P].splice(i, 0, Q);
          //     }
          //   }
          // }
        }
      }
    }

    // for (let P = 0; P < this.nsize; P++) {
    //   result[P].sort((a, b) => a - b);
    // }

  }

  getCsr(portrait, I, J, Val, ) {

    const N = this.nsize * this.udim;

    // let N_free = 0;
    // flag.length = N;
    // for(let i = 0; i < N; i++) {
    //   flag[i] = this.bound[i] === undefined ? N_free++ : -1;
    // }

    I.length = N + 1;
    // I.length = N_free + 1;
    I.fill(0);

    // количество коэффициентов в каждой строке
    for(let P = 0; P < this.nsize; P++) {
      for(let k = 0; k < portrait[P].length; k++) {
        const Q = portrait[P][k];

        for(let i = 0; i < this.udim; i++) {
          I[P*this.udim + i + 1] += this.udim;

          if(P !== Q) {
            I[Q*this.udim + i + 1] += this.udim;
          }
        }
      }
    }

    // начала строк
    for(let i = 0; i < N; i++) {
      I[i + 1] += I[i];
    }

    J.length = I[N];
    Val.length = I[N];
    Val.fill(0);

    const next = I.slice(0, N);

    // столбцы
    for(let P = 0; P < this.nsize; P++) {
      for(let k = 0; k < portrait[P].length; k++) {
        const Q = portrait[P][k];

        for(let i = 0; i < this.udim; i++) {
          const row = P*this.udim + i;

          for(let j = 0; j < this.udim; j++) {
            const col = Q*this.udim + j;

            J[next[row]++] = col;

            if(P !== Q) {
              J[next[col]++] = row;
            }
          }
        }
      }
    }

  }

  getCsrIndex(row, col, I, J) {

    let left = I[row];
    let right = I[row + 1] - 1;

    while(left <= right) {
      const mid = Math.floor((left + right) / 2);

      if(J[mid] === col) {
        return mid;
      }

      if(J[mid] < col) {
        left = mid + 1;
      }
      else{
        right = mid - 1;
      }
    }

    throw new RangeError("Позиция отсутствует в CSR");

  }

  applyKinematicCsr(N, u, f, I, J, Val) {

    let nz = 0;
    let start = I[0];

    for(let p = 0; p < N; p++) {
      const end = I[p + 1]; // старый конец строки
      I[p] = nz;            // новое начало строки

      for(let k = start; k < end; k++) {
        const q = J[k];
        const value = Val[k];

        if(u[p] !== undefined) {
          if (q !== p)
            continue;
          f[p] = value * u[p];
        }
        else if(u[q] !== undefined) {
          f[p] -= value * u[q];
          continue;
        }

        if (value === 0)
          continue;

        J[nz] = q;
        Val[nz] = value;
        nz++;
      }

      start = end;
    }

    I[N] = nz;
    J.length = nz;
    Val.length = nz;

  }


//   applyKinematicCsr(N, u, f, I, J, Val) {
//
//     for(let p = 0; p < N; p++) {
//       if(u[p] !== undefined) {
//         for(let k = I[p]; k < I[p + 1]; k++) {
//           if(J[k] === p) {
//             f[p] = Val[k] * u[p];
//           }
//           else {
//             Val[k] = 0;
//           }
//         }
//       }
//       else{
//         for(let k = I[p]; k < I[p + 1]; k++) {
//           const q = J[k];
//           if(u[q] !== undefined) {
//             f[p] -= Val[k] * u[q];
//             Val[k] = 0;
//           }
//         }
//       }
//     }
//
//   }

  printSCR_Matrix_new(M, N, nz, Val, I, J) {

    const row = new Array(N);

    console.log('Matrix from CSR');

    for(let i = 0; i < M; i++) {
      row.fill(0);

      for(let k = I[i]; k < I[i + 1]; k++) {
        row[J[k]] = Val[k];
      }

      let line = '';

      for(let j = 0; j < N; j++) {
        const value = row[j] === 0 ? '-' : row[j].toExponential(1);
        line += value.padStart(8);
      }

      console.log(line);
    }

  }


  solve() {

    const portrait = [];
    this.getPortrait(portrait);
    console.log(portrait);
    const I = [];
    const J = [];
    const Val = [];

    this.getCsr(portrait, I, J, Val);

    const v = this.elastic[1];
    const E = this.elastic[0];

    let lambda;
    if (this.flag === "ПД") {
      lambda = v*E / ((1+v) * (1-2*v));
    } else {
      lambda = v*E / (1-v*v);
    }
    const mu = E / ( 2*(1+v));

    const N = this.nsize * this.udim;

    const load = this.load.slice();

    const gauss = [-1 / Math.sqrt(3), 1 / Math.sqrt(3)];
    const weight = [1, 1];

    const gaussDeriv = new Array(gauss.length * gauss.length);
    for(let i = 0; i < gauss.length; i++) {
      for(let j = 0; j < gauss.length; j++) {
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

              // globalStiffness[globalRow * N + globalCol] +=
              //   elemStiffness[row * this.kdim + col];
              const index = this.getCsrIndex(globalRow, globalCol, I, J);
              Val[index] += elemStiffness[row * this.kdim + col];
            }
          }
        }
      }
    }

    this.load = load;
    console.log('Matrix before applying kinematic')
    this.printSCR_Matrix_new(N, N, Val.length, Val, I, J);

    const bound = this.bound.slice();

    this.applyKinematicCsr(N, bound, load, I, J, Val);
    console.log('Matrix after applying kinematic')
    this.printSCR_Matrix_new(N, N, Val.length, Val, I, J);

    const result = new Array(N).fill(0);

    const { succesfull, error } =
        cg(N, I, J, Val, load, 1, result);

    if(!succesfull || !result.every(Number.isFinite)) {
      throw new Error(`cg: решение не получено, error = ${error}`);
    };

    printResult(this.ndim, this.node, result, 14, 6);

    // в центре конечного элемента
    const centerResults = [];
    this.getFormDeriv(deriv, 0, 0);
    for(let ie = 0; ie < this.esize; ie++) {
      for(let p = 0; p < this.edim; p++) {
        P[p] = this.elem[ie * this.edim + p];
        x[p] = this.node[P[p] * this.ndim];
        y[p] = this.node[P[p] * this.ndim + 1];
      }
      this.getJacMartix_opt(jacMatrix, deriv, x, y);
      const det = this.getJacobian(jacMatrix);
      this.getJacInverse(jacInverse, jacMatrix);

      let eps_11 = 0;
      let eps_22 = 0;
      let eps_12 = 0;

      for (let p = 0; p < this.edim; p++) {
        const dNx = jacInverse[0] * deriv[p]
                  + jacInverse[2] * deriv[this.edim + p];

        const dNy = jacInverse[1] * deriv[p]
                  + jacInverse[3] * deriv[this.edim + p];

        const ux = result[P[p] * this.udim];
        const uy = result[P[p] * this.udim + 1];

        eps_11 += dNx * ux;
        eps_22 += dNy * uy;
        eps_12 += 0.5 * (dNy * ux + dNx * uy);
      }

      centerResults.push({
        ie,
        x: x.reduce((sum, value) => sum + value, 0) / this.edim,
        y: y.reduce((sum, value) => sum + value, 0) / this.edim,

        eps_11,
        eps_22,
        eps_12,

        sigma_11: (lambda + 2 * mu) * eps_11 + lambda * eps_22,
        sigma_22: lambda * eps_11 + (lambda + 2 * mu) * eps_22,
        sigma_12: 2 * mu * eps_12,
      });
    }

    const e = value => value.toExponential(3).padStart(14);
    console.log(
      "  ie             x             y        eps_11        eps_22        eps_12      sigma_11      sigma_22      sigma_12"
    );
    for(const r of centerResults) {
      console.log(
        String(r.ie).padStart(4) +
        e(r.x) + e(r.y) +
        e(r.eps_11) + e(r.eps_22) + e(r.eps_12) +
        e(r.sigma_11) + e(r.sigma_22) + e(r.sigma_12)
      );
    }

    return 0;
  }

}

export { Isoparametric2d }
