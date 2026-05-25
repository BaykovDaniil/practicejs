const totalLength = 12, FractionDigits = 3;

const nume = (num, fd = FractionDigits) =>
  Number(num).toExponential(fd).replace(/e([+-])(\d)$/, 'e$10$2');
const numf = (num, fd = FractionDigits) => Number(num).toFixed(fd);

const setline = (num, tl = totalLength) => {
  let value = String(num), sline = "";
  if ((tl - value.length) < 0) {
    sline += ' ';
  } else for (let jj = 0; jj < (tl - value.length); jj++) {
    sline += ' ';
  }
  sline += value;
  return sline;
}

export const logm = (title, size, array, vector = [], tl, fd) => {
  title && title.length > 0 && console.log(title);
  let sline = '', jj = 0;
  for (let ii = 0; ii < array.length; ii++) {
    sline += setline(numf(array[ii], fd), tl);
    if ((ii + 1) % size === 0) {
      sline += vector.length > 0 ? '|' + setline(numf(vector[jj], fd), tl) : '';
      console.log(sline);
      sline = '';
      jj++;
    }
  }
}

export const printSCR = (val, I, J, tl, fd) => {
  console.log(`nz=${J.length} I[${I.length}] J[${J.length}] val[${val.length}]`);
  console.log(setline(`I`, 3, 0) + setline(`J`, 3, 0) + setline(`val`, Math.max(tl, 6), fd));
  for (let i = 0; i < J.length; i++) {
    const rowStart = i < I.length ? I[i] : '';
    console.log(setline(rowStart, 3, 0) + setline(J[i], 3, 0) + setline(numf(val[i], fd), tl));
  }
}


export const printResult = (ndim, node, result, tl, fd) => {
  const header = setline(`i`, 3, 0) + setline(`x0`, tl, 1) + setline(`x1`, tl, 1)
    + setline(`u0`, tl, fd) + setline(`u1`, tl, fd);
  console.log(header);

  for (let p = 0; p < result.length / 2; p++) {
    const u = [result[ndim * p + 0], result[ndim * p + 1]];
    const x = [node[ndim * p + 0], node[ndim * p + 1]];
    const line = setline(p, 3, 0) + setline(numf(x[0], 3), tl) + setline(numf(x[1], 3), tl)
     + setline(nume(u[0], fd), tl) + setline(nume(u[1], fd), tl);
    console.log(line);
  }
}

export const axialForcesPrint = (force) =>
{
  const esize = force.length;
  console.log(setline(`ie`, 3, 0) + setline(`axial force`, 12));
  for (let ie = 0; ie < esize; ie++)
  {
    console.log(setline(ie, 3, 0) + setline(numf(force[ie], 3), 12));
  }
}

// ----------------- conjugent gradient method --------------------------------

//cublasDaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
//cublasDscal(cublasHandle, N, &b, d_p, 1);
const hostDaxpy = (N, alpha, beta, a, b) => {
  for (let i = 0; i < N; i++) {
    b[i] = alpha * a[i] + beta * b[i];
  }
}

//cublasDcopy(cublasHandle, N, d_r, 1, d_p, 1)
const hostDcopy = (N, a, b) => {
  hostDaxpy(N, 1.0, 0.0, a, b);
}

const hostCusparseSpMV = (N, row, col, val, vec, result) => {

  for (let i = 0; i < N; i++) {
    result[i] = 0.0;

    for (let j = row[i]; j < row[i + 1]; j++) {
      result[i] += val[j] * vec[col[j]];
    }
  }
}

//cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
const hostDdot = (N, a, b) => {
  let result = 0.0;
  for (let i = 0; i < N; i++) {
    result += b[i] * a[i];
  }
  return result;
}

const free = (a, aBuffer) => {
  a = undefined;
  aBuffer = undefined;
}

export const cg = (N, row, col, val, b, pN, x) => {
  const rBuffer = new ArrayBuffer(N * 8);
  const r = new Float64Array(rBuffer);
  hostDcopy(N, b, r);

  const pBuffer = new ArrayBuffer(N * 8);
  const p = new Float64Array(pBuffer);
  const ApBuffer = new ArrayBuffer(N * 8);
  const Ap = new Float64Array(ApBuffer);

  let r0, r1, dot, alpha, beta;
  let k = 1;

  r0 = 0.0;
  hostCusparseSpMV(N, row, col, val, x, Ap);
  hostDaxpy(N, -1.0, 1.0, Ap, r);
  r1 = hostDdot(N, r, r);

  const tol = window.tol || 1e-10;//1.0e-4;//1e-10;
  while (r1 > (tol * tol) && k <= N) {
    if (k > 1) {
      beta = r1 / r0;
      hostDaxpy(N, 1.0, beta, r, p);
    }
    else {
      hostDcopy(N, r, p);
    }
    hostCusparseSpMV(N, row, col, val, p, Ap);
    dot = hostDdot(N, p, Ap);
    alpha = r1 / dot;
    hostDaxpy(N, alpha, 1.0, p, x);
    hostDaxpy(N, -alpha, 1.0, Ap, r);
    r0 = r1;
    r1 = hostDdot(N, r, r);
    if (k % pN === 0 || k === N || r1 <= tol * tol) {
      console.log(`iteration = ${setline(k, 5)}, residual = ${setline(nume(Math.sqrt(r1), 6), 12)}, r1/r0 = ${setline(nume(r1 / r0, 6), 12)}`);
    }
    k++;
  }
  let err = 0.0;
  // test
  // hostCusparseSpMV(N, row, col, val, x, r);
  // hostDaxpy(N, 1.0, -1.0, b, r);
  // for (let i = 0; i < N; i++) {
  //   let value = Math.abs(r[i]);
  //   if (value > err)
  //     err = value;
  // }
  //

  for (let i = 0; i < N; i++)
    {
      let value = Math.abs(r[i]);
      if (value > err)
        err = value;
    }


  free(r, rBuffer);
  free(p, pBuffer);
  free(Ap, ApBuffer);
  return (k <= N && err < tol) ? { succesfull: true, error: err } : { succesfull: false, error: err };
}

// ------- Base class for fem problem ---
