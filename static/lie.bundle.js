// ../HomogeneousTools/Lie.js/src/DynkinType.js
var _E_POS_ROOTS = { 6: 36, 7: 63, 8: 120 };
var SimpleDynkinType = class {
  /**
   * @param {string} series  One of 'A','B','C','D','E','F','G'.
   * @param {number} n       Rank.
   */
  constructor(series, n) {
    series = series.toUpperCase();
    n = Number(n);
    switch (series) {
      case "A":
        if (n < 1) throw new RangeError(`A_n requires n >= 1, got ${n}`);
        break;
      case "B":
        if (n < 2) throw new RangeError(`B_n requires n >= 2, got ${n}`);
        break;
      case "C":
        if (n < 2) throw new RangeError(`C_n requires n >= 2, got ${n}`);
        break;
      case "D":
        if (n < 4) throw new RangeError(`D_n requires n >= 4, got ${n}`);
        break;
      case "E":
        if (![6, 7, 8].includes(n)) throw new RangeError(`E_n requires n in {6,7,8}, got ${n}`);
        break;
      case "F":
        if (n !== 4) throw new RangeError(`F_n requires n = 4, got ${n}`);
        break;
      case "G":
        if (n !== 2) throw new RangeError(`G_n requires n = 2, got ${n}`);
        break;
      default:
        throw new RangeError(`Unknown Dynkin series '${series}'`);
    }
    this.series = series;
    this.n = n;
    this._key = `${series}${n}`;
  }
  /** Rank of this simple type. */
  get rank() {
    return this.n;
  }
  /** Number of positive roots. */
  get nPositiveRoots() {
    switch (this.series) {
      case "A":
        return this.n * (this.n + 1) / 2;
      case "B":
        return this.n * this.n;
      case "C":
        return this.n * this.n;
      case "D":
        return this.n * (this.n - 1);
      case "E":
        return _E_POS_ROOTS[this.n];
      case "F":
        return 24;
      case "G":
        return 6;
    }
  }
  /** Dimension of the Lie group/algebra: rank + 2 * nPositiveRoots. */
  get dimension() {
    return this.n + 2 * this.nPositiveRoots;
  }
  /** Number of simple factors (always 1 for simple type). */
  get nComponents() {
    return 1;
  }
  /** String key for memoization. */
  key() {
    return this._key;
  }
  toString() {
    return this._key;
  }
};
var ProductDynkinType = class {
  /**
   * @param {SimpleDynkinType[]} factors  Array of simple factors.
   */
  constructor(factors) {
    if (!Array.isArray(factors) || factors.length === 0) {
      throw new Error("ProductDynkinType requires a non-empty array of SimpleDynkinType");
    }
    this.factors = factors;
    const offsets = new Array(factors.length);
    let off = 0;
    for (let i = 0; i < factors.length; i++) {
      offsets[i] = off;
      off += factors[i].rank;
    }
    this.offsets = offsets;
    this._key = factors.map((f) => f.key()).join("x");
  }
  /** Total rank. */
  get rank() {
    let r = 0;
    for (const f of this.factors) r += f.rank;
    return r;
  }
  /** Total number of positive roots. */
  get nPositiveRoots() {
    let s = 0;
    for (const f of this.factors) s += f.nPositiveRoots;
    return s;
  }
  /** Dimension of the semisimple Lie group/algebra. */
  get dimension() {
    return this.rank + 2 * this.nPositiveRoots;
  }
  /** Number of simple factors. */
  get nComponents() {
    return this.factors.length;
  }
  /** String key for memoization. */
  key() {
    return this._key;
  }
  toString() {
    return this.factors.map((f) => f.toString()).join(" \xD7 ");
  }
};
function parseDynkinType(str) {
  const parts = str.split(/\s*[x×]\s*/i).filter(Boolean);
  const factors = parts.map((p) => {
    const m = p.match(/^([A-Ga-g])(\d+)$/);
    if (!m) throw new SyntaxError(`Cannot parse Dynkin type '${p}'`);
    return new SimpleDynkinType(m[1], parseInt(m[2], 10));
  });
  return factors.length === 1 ? factors[0] : new ProductDynkinType(factors);
}

// ../HomogeneousTools/Lie.js/src/Rational.js
function gcd(a, b) {
  a = Math.abs(a);
  b = Math.abs(b);
  while (b) {
    [a, b] = [b, a % b];
  }
  return a;
}
function lcm(a, b) {
  if (a === 0 && b === 0) return 0;
  return Math.abs(a / gcd(a, b)) * Math.abs(b);
}
function bigGcd(a, b) {
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;
  while (b) {
    [a, b] = [b, a % b];
  }
  return a;
}
var Rational = class _Rational {
  /**
   * @param {number} p  Numerator.
   * @param {number} q  Denominator (must be nonzero).
   */
  constructor(p, q = 1) {
    if (q === 0) throw new Error("Rational: denominator is zero");
    if (q < 0) {
      p = -p;
      q = -q;
    }
    const g = gcd(Math.abs(p), q);
    this.p = p / g;
    this.q = q / g;
  }
  /** @returns {Rational} */
  add(r) {
    return new _Rational(this.p * r.q + r.p * this.q, this.q * r.q);
  }
  /** @returns {Rational} */
  sub(r) {
    return new _Rational(this.p * r.q - r.p * this.q, this.q * r.q);
  }
  /** @returns {Rational} */
  mul(r) {
    return new _Rational(this.p * r.p, this.q * r.q);
  }
  /** @returns {Rational} */
  div(r) {
    return new _Rational(this.p * r.q, this.q * r.p);
  }
  /** @returns {Rational} */
  neg() {
    return new _Rational(-this.p, this.q);
  }
  /** @returns {number} */
  toNumber() {
    return this.p / this.q;
  }
  /** @returns {boolean} */
  isZero() {
    return this.p === 0;
  }
  /** @returns {boolean} */
  eq(r) {
    return this.p === r.p && this.q === r.q;
  }
  toString() {
    return this.q === 1 ? `${this.p}` : `${this.p}/${this.q}`;
  }
};
var rat = (p, q = 1) => new Rational(p, q);
var ZERO = Object.freeze(new Rational(0));
var ONE = Object.freeze(new Rational(1));

// ../HomogeneousTools/Lie.js/src/CartanMatrix.js
var _cartanCache = /* @__PURE__ */ new Map();
var _symCache = /* @__PURE__ */ new Map();
var _bilinearCache = /* @__PURE__ */ new Map();
var _inverseCache = /* @__PURE__ */ new Map();
var _omegaBilinearCache = /* @__PURE__ */ new Map();
function _simpleCartanMatrix(dt) {
  const n = dt.rank;
  const C = new Int32Array(n * n);
  for (let i = 0; i < n; i++) C[i * n + i] = 2;
  switch (dt.series) {
    case "A":
      for (let i = 0; i < n - 1; i++) {
        C[i * n + (i + 1)] = -1;
        C[(i + 1) * n + i] = -1;
      }
      break;
    case "B":
      for (let i = 0; i < n - 1; i++) {
        C[i * n + (i + 1)] = -1;
        C[(i + 1) * n + i] = -1;
      }
      C[(n - 1) * n + (n - 2)] = -2;
      break;
    case "C":
      for (let i = 0; i < n - 1; i++) {
        C[i * n + (i + 1)] = -1;
        C[(i + 1) * n + i] = -1;
      }
      C[(n - 2) * n + (n - 1)] = -2;
      break;
    case "D":
      for (let i = 0; i < n - 2; i++) {
        if (i < n - 3) {
          C[i * n + (i + 1)] = -1;
          C[(i + 1) * n + i] = -1;
        }
      }
      for (let i = 0; i < n - 3; i++) {
        C[i * n + (i + 1)] = -1;
        C[(i + 1) * n + i] = -1;
      }
      C[(n - 3) * n + (n - 2)] = -1;
      C[(n - 2) * n + (n - 3)] = -1;
      C[(n - 3) * n + (n - 1)] = -1;
      C[(n - 1) * n + (n - 3)] = -1;
      break;
    case "E": {
      const edges = [[0, 2]];
      for (let i = 2; i < n - 1; i++) edges.push([i, i + 1]);
      edges.push([1, 3]);
      for (const [a, b] of edges) {
        C[a * n + b] = -1;
        C[b * n + a] = -1;
      }
      break;
    }
    case "F":
      C[0 * n + 1] = -1;
      C[1 * n + 0] = -1;
      C[1 * n + 2] = -1;
      C[2 * n + 1] = -2;
      C[2 * n + 3] = -1;
      C[3 * n + 2] = -1;
      break;
    case "G":
      C[0 * n + 1] = -3;
      C[1 * n + 0] = -1;
      break;
  }
  return C;
}
function cartanMatrix(dt) {
  const k = dt.key();
  if (_cartanCache.has(k)) return _cartanCache.get(k);
  let C;
  if (dt.factors !== void 0) {
    const r = dt.rank;
    C = new Int32Array(r * r);
    for (const fac of dt.factors) {
      const off = dt.offsets[dt.factors.indexOf(fac)];
      const block = _simpleCartanMatrix(fac);
      const fr = fac.rank;
      for (let i = 0; i < fr; i++) {
        for (let j = 0; j < fr; j++) {
          C[(off + i) * r + (off + j)] = block[i * fr + j];
        }
      }
    }
  } else {
    C = _simpleCartanMatrix(dt);
  }
  _cartanCache.set(k, C);
  return C;
}
function cartanSymmetrizer(dt) {
  const k = dt.key();
  if (_symCache.has(k)) return _symCache.get(k);
  if (dt.factors !== void 0) {
    const parts = dt.factors.map((f) => cartanSymmetrizer(f));
    const r = dt.rank;
    const d2 = new Int32Array(r);
    let off = 0;
    for (let idx = 0; idx < dt.factors.length; idx++) {
      const p = parts[idx];
      for (let i = 0; i < p.length; i++) d2[off + i] = p[i];
      off += dt.factors[idx].rank;
    }
    _symCache.set(k, d2);
    return d2;
  }
  const n = dt.rank;
  const C = cartanMatrix(dt);
  const dRat = new Array(n).fill(null);
  dRat[0] = rat(1);
  const visited = new Uint8Array(n);
  visited[0] = 1;
  const queue = [0];
  while (queue.length > 0) {
    const i = queue.shift();
    for (let j = 0; j < n; j++) {
      if (!visited[j] && C[i * n + j] !== 0) {
        const cij = C[i * n + j];
        const cji = C[j * n + i];
        dRat[j] = dRat[i].mul(rat(Math.abs(cij), Math.abs(cji)));
        visited[j] = 1;
        queue.push(j);
      }
    }
  }
  let lcmDen = 1;
  for (const r of dRat) lcmDen = lcm(lcmDen, r.q);
  const d = new Int32Array(n);
  for (let i = 0; i < n; i++) d[i] = dRat[i].p * (lcmDen / dRat[i].q);
  let g = d[0];
  for (let i = 1; i < n; i++) g = gcd(g, d[i]);
  for (let i = 0; i < n; i++) d[i] /= g;
  _symCache.set(k, d);
  return d;
}
function cartanBilinearForm(dt) {
  const k = dt.key();
  if (_bilinearCache.has(k)) return _bilinearCache.get(k);
  const n = dt.rank;
  const C = cartanMatrix(dt);
  const d = cartanSymmetrizer(dt);
  const B = new Int32Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      B[i * n + j] = d[i] * C[i * n + j];
    }
  }
  _bilinearCache.set(k, B);
  return B;
}
function cartanMatrixInverse(dt) {
  const k = dt.key();
  if (_inverseCache.has(k)) return _inverseCache.get(k);
  const n = dt.rank;
  const C = cartanMatrix(dt);
  const M = new Array(n);
  for (let i = 0; i < n; i++) {
    M[i] = new Array(2 * n);
    for (let j = 0; j < n; j++) M[i][j] = rat(C[i * n + j]);
    for (let j = 0; j < n; j++) M[i][n + j] = i === j ? rat(1) : rat(0);
  }
  for (let col = 0; col < n; col++) {
    let pivot = -1;
    for (let row = col; row < n; row++) {
      if (!M[row][col].isZero()) {
        pivot = row;
        break;
      }
    }
    if (pivot < 0) throw new Error("Cartan matrix is singular");
    if (pivot !== col) [M[col], M[pivot]] = [M[pivot], M[col]];
    const scale = M[col][col];
    for (let j = 0; j < 2 * n; j++) M[col][j] = M[col][j].div(scale);
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const f = M[row][col];
      if (f.isZero()) continue;
      for (let j = 0; j < 2 * n; j++) {
        M[row][j] = M[row][j].sub(f.mul(M[col][j]));
      }
    }
  }
  const inv = new Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      inv[i * n + j] = M[i][n + j];
    }
  }
  _inverseCache.set(k, inv);
  return inv;
}
function omegaBilinearFormScaled(dt) {
  const k = dt.key();
  if (_omegaBilinearCache.has(k)) return _omegaBilinearCache.get(k);
  const n = dt.rank;
  const Cinv = cartanMatrixInverse(dt);
  const d = cartanSymmetrizer(dt);
  const C = cartanMatrix(dt);
  const BCinv = new Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let s = rat(0);
      for (let k2 = 0; k2 < n; k2++) {
        s = s.add(rat(d[i] * C[i * n + k2]).mul(Cinv[k2 * n + j]));
      }
      BCinv[i * n + j] = s;
    }
  }
  const Bomega = new Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let s = rat(0);
      for (let k2 = 0; k2 < n; k2++) {
        s = s.add(Cinv[k2 * n + i].mul(BCinv[k2 * n + j]));
      }
      Bomega[i * n + j] = s;
    }
  }
  let scale = 1;
  for (let i = 0; i < n * n; i++) scale = lcm(scale, Bomega[i].q);
  const matrix = new Int32Array(n * n);
  for (let i = 0; i < n * n; i++) matrix[i] = Bomega[i].p * (scale / Bomega[i].q);
  const result = { scale, matrix };
  _omegaBilinearCache.set(k, result);
  return result;
}
function cartanDeterminant(dt) {
  if (dt.factors !== void 0) {
    let d = 1;
    for (const f of dt.factors) d *= cartanDeterminant(f);
    return d;
  }
  switch (dt.series) {
    case "A":
      return dt.n + 1;
    case "B":
      return 2;
    case "C":
      return 2;
    case "D":
      return 4;
    case "E":
      if (dt.n === 6) return 3;
      if (dt.n === 7) return 2;
      return 1;
    case "F":
      return 1;
    case "G":
      return 1;
  }
}

// ../HomogeneousTools/Lie.js/src/RootSystem.js
var _rsCache = /* @__PURE__ */ new Map();
var NEG = 4294967295;
function _computePositiveRootsAndReflections(C, n) {
  const posRoots = [];
  const posCoroots = [];
  for (let i = 0; i < n; i++) {
    const r = new Int32Array(n);
    r[i] = 1;
    posRoots.push(r);
    const cr = new Int32Array(n);
    cr[i] = 1;
    posCoroots.push(cr);
  }
  const rootIndex = /* @__PURE__ */ new Map();
  for (let i = 0; i < n; i++) rootIndex.set(posRoots[i].join(","), i);
  const reflData = /* @__PURE__ */ new Map();
  for (let s = 0; s < n; s++) reflData.set(`${s},${s}`, NEG);
  let idx = 0;
  while (idx < posRoots.length) {
    const rootI = posRoots[idx];
    const corootI = posCoroots[idx];
    for (let s = 0; s < n; s++) {
      const key = `${s},${idx}`;
      if (reflData.has(key)) continue;
      let pairing = 0;
      for (let j = 0; j < n; j++) pairing += C[s * n + j] * rootI[j];
      let copairing = 0;
      for (let j = 0; j < n; j++) copairing += corootI[j] * C[j * n + s];
      if (pairing * copairing >= 4) {
        reflData.set(key, NEG);
        continue;
      }
      const newRoot = new Int32Array(n);
      for (let j = 0; j < n; j++) newRoot[j] = rootI[j] - (j === s ? pairing : 0);
      const newCoroot = new Int32Array(n);
      for (let j = 0; j < n; j++) newCoroot[j] = corootI[j] - (j === s ? copairing : 0);
      const allNonNeg = newRoot.every((v) => v >= 0);
      const allNonPos = newRoot.every((v) => v <= 0);
      if (allNonNeg) {
        const nrKey = newRoot.join(",");
        let nrIdx = rootIndex.get(nrKey);
        if (nrIdx === void 0) {
          nrIdx = posRoots.length;
          posRoots.push(newRoot);
          posCoroots.push(newCoroot);
          rootIndex.set(nrKey, nrIdx);
        }
        reflData.set(key, nrIdx);
        reflData.set(`${s},${nrIdx}`, idx);
      } else {
        reflData.set(key, NEG);
      }
    }
    idx++;
  }
  const nPos = posRoots.length;
  const heights = posRoots.map((r) => {
    let s = 0;
    for (let i = 0; i < n; i++) s += r[i];
    return s;
  });
  const perm = Array.from({ length: nPos }, (_, i) => i);
  perm.sort((a, b) => heights[a] - heights[b] || a - b);
  const invPerm = new Uint32Array(nPos);
  for (let i = 0; i < nPos; i++) invPerm[perm[i]] = i;
  const sortedRoots = perm.map((i) => posRoots[i]);
  const sortedCoroots = perm.map((i) => posCoroots[i]);
  const refl = new Array(n);
  for (let s = 0; s < n; s++) {
    refl[s] = new Uint32Array(nPos);
    for (let i = 0; i < nPos; i++) {
      const oldI = perm[i];
      const rkey = `${s},${oldI}`;
      const oldV = reflData.get(rkey);
      refl[s][i] = oldV === void 0 || oldV === NEG ? NEG : invPerm[oldV];
    }
  }
  return { positiveRoots: sortedRoots, positiveCoroots: sortedCoroots, refl };
}
function rootSystem(dt) {
  const k = dt.key();
  if (_rsCache.has(k)) return _rsCache.get(k);
  const n = dt.rank;
  const C = cartanMatrix(dt);
  const { positiveRoots, positiveCoroots, refl } = _computePositiveRootsAndReflections(C, n);
  const nPos = positiveRoots.length;
  let hrIdx = 0;
  let hrHeight = 0;
  for (let i = 0; i < nPos; i++) {
    let h = 0;
    for (let j = 0; j < n; j++) h += positiveRoots[i][j];
    if (h > hrHeight) {
      hrHeight = h;
      hrIdx = i;
    }
  }
  let hcrIdx = 0;
  let hcrHeight = 0;
  for (let i = 0; i < nPos; i++) {
    let h = 0;
    for (let j = 0; j < n; j++) h += positiveCoroots[i][j];
    if (h > hcrHeight) {
      hcrHeight = h;
      hcrIdx = i;
    }
  }
  const rs = {
    rank: n,
    nPositiveRoots: nPos,
    positiveRoots,
    positiveCoroots,
    refl,
    highestRootIdx: hrIdx,
    highestCorootIdx: hcrIdx,
    _dt: dt
  };
  _rsCache.set(k, rs);
  return rs;
}
function highestRoot(rs) {
  return rs.positiveRoots[rs.highestRootIdx];
}
function highestShortRoot(rs) {
  return rs.positiveRoots[rs.highestCorootIdx];
}
function rootHeight(v) {
  let h = 0;
  for (let i = 0; i < v.length; i++) h += v[i];
  return h;
}
function coxeterCoefficients(dt) {
  if (dt.factors !== void 0) {
    const r = dt.rank;
    const result = new Int32Array(r);
    let off = 0;
    for (const f of dt.factors) {
      const cc = coxeterCoefficients(f);
      for (let i = 0; i < f.rank; i++) result[off + i] = cc[i];
      off += f.rank;
    }
    return result;
  }
  const rs = rootSystem(dt);
  return Int32Array.from(highestRoot(rs));
}
function dualCoxeterCoefficients(dt) {
  if (dt.factors !== void 0) {
    const r = dt.rank;
    const result = new Int32Array(r);
    let off = 0;
    for (const f of dt.factors) {
      const dc = dualCoxeterCoefficients(f);
      for (let i = 0; i < f.rank; i++) result[off + i] = dc[i];
      off += f.rank;
    }
    return result;
  }
  const rs = rootSystem(dt);
  return Int32Array.from(rs.positiveCoroots[rs.highestRootIdx]);
}
function coxeterNumber(dt) {
  const cc = coxeterCoefficients(dt);
  let s = 1;
  for (let i = 0; i < cc.length; i++) s += cc[i];
  return s;
}
function dualCoxeterNumber(dt) {
  const dc = dualCoxeterCoefficients(dt);
  let s = 1;
  for (let i = 0; i < dc.length; i++) s += dc[i];
  return s;
}
var _degrees = {
  A: (n) => Array.from({ length: n }, (_, i) => i + 2),
  B: (n) => Array.from({ length: n }, (_, i) => 2 * (i + 1)),
  C: (n) => Array.from({ length: n }, (_, i) => 2 * (i + 1)),
  D: (n) => {
    const d = Array.from({ length: n - 1 }, (_, i) => 2 * (i + 1));
    d.push(n);
    return d;
  },
  E6: [2, 5, 6, 8, 9, 12],
  E7: [2, 6, 8, 10, 12, 14, 18],
  E8: [2, 8, 12, 14, 18, 20, 24, 30],
  F: [2, 6, 8, 12],
  G: [2, 6]
};
function degreesFundamentalInvariants(dt) {
  if (dt.factors !== void 0) {
    return dt.factors.flatMap((f) => degreesFundamentalInvariants(f));
  }
  const { series, n } = dt;
  if (series === "E") return _degrees[`E${n}`].slice();
  if (series === "F" || series === "G") return _degrees[series].slice();
  return _degrees[series](n);
}

// ../HomogeneousTools/Lie.js/src/WeightLattice.js
function fundamentalWeight(rank, i) {
  const w = new Int32Array(rank);
  w[i] = 1;
  return w;
}
function weylVector(rank) {
  const w = new Int32Array(rank);
  for (let i = 0; i < rank; i++) w[i] = 1;
  return w;
}
function isDominant(w) {
  for (let i = 0; i < w.length; i++) if (w[i] < 0) return false;
  return true;
}
function reflect(w, s, C, rank) {
  const out = new Int32Array(rank);
  const pairing = w[s];
  for (let j = 0; j < rank; j++) {
    out[j] = w[j] - pairing * C[j * rank + s];
  }
  return out;
}
function reflectInPlace(w, s, C, rank) {
  const pairing = w[s];
  for (let j = 0; j < rank; j++) {
    w[j] -= pairing * C[j * rank + s];
  }
}
function conjugateDominantWeight(w, C, rank) {
  const v = Int32Array.from(w);
  let s = 0;
  while (s < rank) {
    if (v[s] < 0) {
      const p = v[s];
      for (let j = 0; j < rank; j++) v[j] -= p * C[j * rank + s];
      s = 0;
    } else {
      s++;
    }
  }
  return v;
}
function conjugateDominantWeightWithLength(w, C, rank) {
  const v = Int32Array.from(w);
  let len = 0;
  let s = 0;
  while (s < rank) {
    if (v[s] < 0) {
      const p = v[s];
      for (let j = 0; j < rank; j++) v[j] -= p * C[j * rank + s];
      len++;
      s = 0;
    } else {
      s++;
    }
  }
  return { weight: v, length: len };
}
function weightAdd(a, b) {
  const n = a.length;
  const c = new Int32Array(n);
  for (let i = 0; i < n; i++) c[i] = a[i] + b[i];
  return c;
}
function weightSub(a, b) {
  const n = a.length;
  const c = new Int32Array(n);
  for (let i = 0; i < n; i++) c[i] = a[i] - b[i];
  return c;
}
function weightScale(w, k) {
  const n = w.length;
  const c = new Int32Array(n);
  for (let i = 0; i < n; i++) c[i] = k * w[i];
  return c;
}
function weightNeg(w) {
  const n = w.length;
  const c = new Int32Array(n);
  for (let i = 0; i < n; i++) c[i] = -w[i];
  return c;
}
function weightKey(w) {
  return w.join(",");
}
function rootWeightPairing(root, weight, d) {
  let s = 0;
  for (let i = 0; i < root.length; i++) s += root[i] * d[i] * weight[i];
  return s;
}
function rootToWeight(root, C, rank) {
  const w = new Int32Array(rank);
  for (let i = 0; i < rank; i++) {
    let s = 0;
    for (let j = 0; j < rank; j++) s += C[i * rank + j] * root[j];
    w[i] = s;
  }
  return w;
}

// ../HomogeneousTools/Lie.js/src/WeylGroup.js
var _denomCache = /* @__PURE__ */ new Map();
var _scaledRootsCache = /* @__PURE__ */ new Map();
var _longestCache = /* @__PURE__ */ new Map();
function _weylDenomAndScaledRoots(dt) {
  const k = dt.key();
  if (_denomCache.has(k)) return { denom: _denomCache.get(k), scaledRoots: _scaledRootsCache.get(k) };
  const rs = rootSystem(dt);
  const n = rs.rank;
  const d = cartanSymmetrizer(dt);
  const scaledRoots = [];
  let denom = 1n;
  for (let ri = 0; ri < rs.nPositiveRoots; ri++) {
    const alpha = rs.positiveRoots[ri];
    const sa = new Int32Array(n);
    let ip = 0;
    for (let i = 0; i < n; i++) {
      sa[i] = d[i] * alpha[i];
      ip += sa[i];
    }
    scaledRoots.push(sa);
    denom *= BigInt(ip);
  }
  _denomCache.set(k, denom);
  _scaledRootsCache.set(k, scaledRoots);
  return { denom, scaledRoots };
}
function weylDimension(dt, hw) {
  const { denom, scaledRoots } = _weylDenomAndScaledRoots(dt);
  const n = dt.rank;
  const rho = weylVector(n);
  const lrho = new Int32Array(n);
  for (let i = 0; i < n; i++) lrho[i] = hw[i] + rho[i];
  let numer = 1n;
  for (const sa of scaledRoots) {
    let ip = 0;
    for (let i = 0; i < n; i++) ip += lrho[i] * sa[i];
    numer *= BigInt(ip);
  }
  return numer / denom;
}
function borelWeylBott(dt, lambda) {
  const n = dt.rank;
  const C = cartanMatrix(dt);
  const mu = new Int32Array(n);
  for (let i = 0; i < n; i++) mu[i] = lambda[i] + 1;
  const { weight: muDom, length: deg } = conjugateDominantWeightWithLength(mu, C, n);
  for (let i = 0; i < n; i++) {
    if (muDom[i] === 0) return null;
  }
  const hw = new Int32Array(n);
  for (let i = 0; i < n; i++) hw[i] = muDom[i] - 1;
  return { degree: deg, weight: hw };
}
function weylOrder(dt) {
  if (dt.factors !== void 0) {
    let o = 1n;
    for (const f of dt.factors) o *= weylOrder(f);
    return o;
  }
  const { series, n } = dt;
  const facN = _factorial(n);
  switch (series) {
    case "A":
      return _factorial(n + 1);
    case "B":
      return facN * 2n ** BigInt(n);
    case "C":
      return facN * 2n ** BigInt(n);
    case "D":
      return facN * 2n ** BigInt(n - 1);
    case "E":
      if (n === 6) return 51840n;
      if (n === 7) return 2903040n;
      return 696729600n;
    case "F":
      return 1152n;
    case "G":
      return 12n;
  }
}
function _factorial(n) {
  let f = 1n;
  for (let i = 2; i <= n; i++) f *= BigInt(i);
  return f;
}
function longestElement(dt) {
  const k = dt.key();
  if (_longestCache.has(k)) return _longestCache.get(k);
  const n = dt.rank;
  const C = cartanMatrix(dt);
  const wt = new Int32Array(n);
  for (let i = 0; i < n; i++) wt[i] = 1;
  const word = [];
  outer: while (true) {
    for (let s = 0; s < n; s++) {
      if (wt[s] > 0) {
        word.push(s);
        const p = wt[s];
        for (let j = 0; j < n; j++) wt[j] -= p * C[j * n + s];
        continue outer;
      }
    }
    break;
  }
  _longestCache.set(k, word);
  return word;
}
function dominantWeights(dt, hw) {
  const rs = rootSystem(dt);
  const n = rs.rank;
  const C = cartanMatrix(dt);
  const posRootsW = rs.positiveRoots.map((alpha) => rootToWeight(alpha, C, n));
  const found = /* @__PURE__ */ new Set();
  const hwKey = weightKey(hw);
  found.add(hwKey);
  let todo = [hw];
  const result = [Int32Array.from(hw)];
  while (todo.length > 0) {
    const newTodo = [];
    for (const w of todo) {
      for (const alphaW of posRootsW) {
        const wSub = weightSub(w, alphaW);
        if (!isDominant(wSub)) continue;
        const k = weightKey(wSub);
        if (found.has(k)) continue;
        found.add(k);
        result.push(wSub);
        newTodo.push(wSub);
      }
    }
    todo = newTodo;
  }
  const Cinv = cartanMatrixInverse(dt);
  const levelOf = (w) => {
    let lvl = 0;
    for (let i = 0; i < n; i++) {
      let s = 0;
      for (let j = 0; j < n; j++) s += Cinv[i * n + j].toNumber() * (hw[j] - w[j]);
      lvl += s;
    }
    return lvl;
  };
  result.sort((a, b) => levelOf(a) - levelOf(b));
  return result;
}
function isSingular(dt, w) {
  const C = cartanMatrix(dt);
  const n = dt.rank;
  const { weight: dom } = conjugateDominantWeightWithLength(w, C, n);
  for (let i = 0; i < n; i++) if (dom[i] === 0) return true;
  return false;
}

// ../HomogeneousTools/Lie.js/src/Weylloop.js
function _w2e_A(w, e, n) {
  let s = 0;
  e[n] = 0;
  for (let i = n - 1; i >= 0; i--) {
    s += w[i];
    e[i] = s;
  }
}
function _e2w_A(e, w, n) {
  for (let i = 0; i < n; i++) w[i] = e[i] - e[i + 1];
}
function _w2e_B(w, e, n) {
  let s = w[n - 1];
  e[n - 1] = s;
  for (let i = n - 2; i >= 0; i--) {
    s += 2 * w[i];
    e[i] = s;
  }
}
function _e2w_B(e, w, n) {
  for (let i = 0; i < n - 1; i++) w[i] = e[i] - e[i + 1] >> 1;
  w[n - 1] = e[n - 1];
}
function _w2e_C(w, e, n) {
  let s = w[n - 1];
  e[n - 1] = s;
  for (let i = n - 2; i >= 0; i--) {
    s += w[i];
    e[i] = s;
  }
}
function _e2w_C(e, w, n) {
  for (let i = 0; i < n - 1; i++) w[i] = e[i] - e[i + 1];
  w[n - 1] = e[n - 1];
}
function _w2e_D(w, e, n) {
  let s = w[n - 1] - w[n - 2];
  e[n - 1] = s;
  for (let i = n - 2; i >= 0; i--) {
    s += 2 * w[i];
    e[i] = s;
  }
}
function _e2w_D(e, w, n) {
  for (let i = 0; i < n - 1; i++) w[i] = e[i] - e[i + 1] >> 1;
  w[n - 1] = w[n - 2] + e[n - 1];
}
function _w2e_E6(w, e) {
  let s = w[1] - w[2];
  e[4] = s;
  const sumsum0 = 4 * w[0] - s;
  s += 2 * w[2];
  e[3] = s;
  let sumsum = sumsum0 + s;
  s += 2 * w[3];
  e[2] = s;
  sumsum += s;
  s += 2 * w[4];
  e[1] = s;
  sumsum += s;
  s += 2 * w[5];
  e[0] = s;
  sumsum += s;
  e[5] = -sumsum;
}
function _e2w_E6(e, w) {
  let s = e[5];
  s += e[0];
  w[5] = e[0] - e[1] >> 1;
  s += e[1];
  w[4] = e[1] - e[2] >> 1;
  s += e[2];
  w[3] = e[2] - e[3] >> 1;
  s += e[3];
  w[2] = e[3] - e[4] >> 1;
  w[1] = w[2] + e[4];
  w[0] = e[4] - s >> 2;
}
function _w2e_E7(w, e) {
  e[7] = 0;
  let s = 0;
  for (let i = 6; i >= 2; i--) {
    s += w[i];
    e[i] = s;
  }
  e[1] = s + w[0];
  e[0] = e[6] + e[5] + e[4] - e[3] - e[2] - e[1] - 2 * w[1];
}
function _e2w_E7(e, w) {
  w[0] = e[1] - e[2];
  for (let i = 2; i <= 6; i++) w[i] = e[i] - e[i + 1];
  w[1] = e[7] + e[6] + e[5] + e[4] - e[3] - e[2] - e[1] - e[0] >> 1;
}
function _w2e_E8(w, e) {
  let s = w[1] - w[2];
  e[7] = s;
  const sumsum0 = 4 * w[0] - s;
  s += 2 * w[2];
  e[6] = s;
  let sumsum = sumsum0 + s;
  s += 2 * w[3];
  e[5] = s;
  sumsum += s;
  s += 2 * w[4];
  e[4] = s;
  sumsum += s;
  s += 2 * w[5];
  e[3] = s;
  sumsum += s;
  s += 2 * w[6];
  e[2] = s;
  sumsum += s;
  s += 2 * w[7];
  e[1] = s;
  sumsum += s;
  e[0] = -sumsum;
}
function _e2w_E8(e, w) {
  let s = e[0];
  s += e[1];
  w[7] = e[1] - e[2] >> 1;
  s += e[2];
  w[6] = e[2] - e[3] >> 1;
  s += e[3];
  w[5] = e[3] - e[4] >> 1;
  s += e[4];
  w[4] = e[4] - e[5] >> 1;
  s += e[5];
  w[3] = e[5] - e[6] >> 1;
  s += e[6];
  w[2] = e[6] - e[7] >> 1;
  w[1] = w[2] + e[7];
  w[0] = e[7] - s >> 2;
}
function _w2e_F4(w, e) {
  e[3] = w[2];
  e[2] = e[3] + 2 * w[1];
  e[1] = e[2] + 2 * w[0];
  e[0] = -2 * w[3] - e[1] - e[2] - e[3];
}
function _e2w_F4(e, w) {
  w[2] = e[3];
  w[1] = e[2] - e[3] >> 1;
  w[0] = e[1] - e[2] >> 1;
  w[3] = -(e[0] + e[1] + e[2] + e[3]) >> 1;
}
function _w2e_G2(w, e) {
  e[2] = 0;
  e[1] = w[1];
  e[0] = -w[0] - w[1];
}
function _e2w_G2(e, w) {
  w[1] = e[1] - e[2];
  w[0] = e[2] - w[1] - e[0];
}
var _typeInfo = {
  A: (n) => ({ subtype: "A", epsDim: n + 1, permSize: n + 1, w2e: (w, e) => _w2e_A(w, e, n), e2w: (e, w) => _e2w_A(e, w, n) }),
  B: (n) => ({ subtype: "B", epsDim: n, permSize: n, w2e: (w, e) => _w2e_B(w, e, n), e2w: (e, w) => _e2w_B(e, w, n) }),
  C: (n) => ({ subtype: "B", epsDim: n, permSize: n, w2e: (w, e) => _w2e_C(w, e, n), e2w: (e, w) => _e2w_C(e, w, n) }),
  D: (n) => ({ subtype: "D", epsDim: n, permSize: n, w2e: (w, e) => _w2e_D(w, e, n), e2w: (e, w) => _e2w_D(e, w, n) }),
  E6: () => ({ subtype: "D", epsDim: 6, permSize: 5, w2e: _w2e_E6, e2w: _e2w_E6 }),
  E7: () => ({ subtype: "A", epsDim: 8, permSize: 8, w2e: _w2e_E7, e2w: _e2w_E7 }),
  E8: () => ({ subtype: "D", epsDim: 8, permSize: 8, w2e: _w2e_E8, e2w: _e2w_E8 }),
  F4: () => ({ subtype: "B", epsDim: 4, permSize: 4, w2e: _w2e_F4, e2w: _e2w_F4 }),
  G2: () => ({ subtype: "A", epsDim: 3, permSize: 3, w2e: _w2e_G2, e2w: _e2w_G2 })
};
function _getTypeInfo(dt) {
  const { series, n } = dt;
  if (series === "E") return _typeInfo[`E${n}`]();
  if (series === "F") return _typeInfo.F4();
  if (series === "G") return _typeInfo.G2();
  return _typeInfo[series](n);
}
function _nextperm(w, N) {
  if (N <= 1) return false;
  let i = N - 2;
  while (i >= 0 && w[i] >= w[i + 1]) i--;
  if (i < 0) return false;
  let j = N - 1;
  while (w[i] >= w[j]) j--;
  [w[i], w[j]] = [w[j], w[i]];
  let lo = i + 1, hi = N - 1;
  while (lo < hi) {
    [w[lo], w[hi]] = [w[hi], w[lo]];
    lo++;
    hi--;
  }
  return true;
}
function _normalform(w, N, subtype) {
  if (subtype === "A") {
    const sub = w.subarray ? w.subarray(0, N) : w.slice(0, N);
    _sortSubarray(w, N);
    if (w[0] !== 0) {
      const c = w[0];
      for (let i = 0; i < N; i++) w[i] -= c;
    }
  } else {
    let parity = 0;
    for (let i = 0; i < N; i++) {
      if (w[i] < 0) {
        w[i] = -w[i];
        parity++;
      }
    }
    _sortSubarray(w, N);
    if (subtype === "D" && parity & 1 && w[0] !== 0) {
      w[0] = -w[0];
    }
  }
}
function _sortSubarray(arr, N) {
  for (let i = 1; i < N; i++) {
    const v = arr[i];
    let j = i - 1;
    while (j >= 0 && arr[j] > v) {
      arr[j + 1] = arr[j];
      j--;
    }
    arr[j + 1] = v;
  }
}
function _weylMatrix(C, rank, word) {
  const M = new Int32Array(rank * rank);
  for (let i = 0; i < rank; i++) M[i * rank + i] = 1;
  for (const k of word) {
    for (let i = 0; i < rank; i++) {
      const p = M[i * rank + k];
      if (p === 0) continue;
      for (let j = 0; j < rank; j++) {
        M[i * rank + j] -= p * C[j * rank + k];
      }
    }
  }
  return M;
}
function _matMul(A, B, n) {
  const R = new Int32Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let s = 0;
      for (let k = 0; k < n; k++) s += A[i * n + k] * B[k * n + j];
      R[i * n + j] = s;
    }
  }
  return R;
}
function _identity(n) {
  const I = new Int32Array(n * n);
  for (let i = 0; i < n; i++) I[i * n + i] = 1;
  return I;
}
var _cosetCache = /* @__PURE__ */ new Map();
function _matKey(M) {
  return M.join(",");
}
function _buildCosetReps(C, rank, coxWord, coxOrder, xWords) {
  const cox = _weylMatrix(C, rank, coxWord);
  const X = xWords.map((w) => _weylMatrix(C, rank, w));
  const all = [];
  let cur = _identity(rank);
  for (let k = 0; k < coxOrder; k++) {
    for (const x of X) {
      all.push(_matMul(cur, x, rank));
    }
    if (k < coxOrder - 1) cur = _matMul(cur, cox, rank);
  }
  const seen = /* @__PURE__ */ new Set();
  const unique = [];
  for (const m of all) {
    const k = _matKey(m);
    if (!seen.has(k)) {
      seen.add(k);
      unique.push(m);
    }
  }
  return unique;
}
function _getCosetReps(dt) {
  const k = dt.key();
  if (_cosetCache.has(k)) return _cosetCache.get(k);
  const C = cartanMatrix(dt);
  const n = dt.rank;
  let reps;
  if ("ABCD".includes(dt.series)) {
    reps = [_identity(n)];
  } else if (dt.series === "E" && dt.n === 6) {
    reps = _buildCosetReps(
      C,
      6,
      [5, 1, 4, 3, 2, 0],
      12,
      [[], [1, 3, 2, 0], [0, 3, 4, 1, 3, 2, 0]]
    );
  } else if (dt.series === "E" && dt.n === 7) {
    const x1w = [6, 5, 4, 3, 2, 1, 0, 4, 2, 3, 1];
    const C7 = C;
    const x1 = _weylMatrix(C7, 7, x1w);
    const x2 = _matMul(x1, x1, 7);
    const x5 = _matMul(_matMul(x2, x2, 7), x1, 7);
    const cox = _weylMatrix(C7, 7, [0, 1, 2, 3, 4, 5, 6]);
    const all = [];
    let cur = _identity(7);
    for (let k2 = 0; k2 < 18; k2++) {
      for (const x of [_identity(7), x1, x2, x5]) {
        all.push(_matMul(cur, x, 7));
      }
      if (k2 < 17) cur = _matMul(cur, cox, 7);
    }
    const seen = /* @__PURE__ */ new Set();
    reps = [];
    for (const m of all) {
      const mk = _matKey(m);
      if (!seen.has(mk)) {
        seen.add(mk);
        reps.push(m);
      }
    }
  } else if (dt.series === "E" && dt.n === 8) {
    reps = _buildCosetReps(
      C,
      8,
      [6, 4, 2, 1, 7, 5, 3, 0],
      15,
      [
        [],
        [6, 3, 2, 1, 3, 4, 5, 0, 2, 3, 4, 1, 3, 2, 0]
      ]
    );
    const x1 = _weylMatrix(C, 8, [6, 3, 2, 1, 3, 4, 5, 0, 2, 3, 4, 1, 3, 2, 0]);
    const x1_2 = _matMul(x1, x1, 8);
    const x1_3 = _matMul(x1_2, x1, 8);
    const x1_4 = _matMul(x1_2, x1_2, 8);
    const x2 = _weylMatrix(C, 8, [7, 5, 6, 4, 5, 2, 3, 4, 1, 3, 2, 0]);
    const x2_2 = _matMul(x2, x2, 8);
    const x2_3 = _matMul(x2_2, x2, 8);
    const x2_4 = _matMul(x2_3, x2, 8);
    const x2_5 = _matMul(x2_4, x2, 8);
    const cox = _weylMatrix(C, 8, [6, 4, 2, 1, 7, 5, 3, 0]);
    const Xmats = [_identity(8), x1, x1_3, x1_4, x2, x2_2, x2_3, x2_4, x2_5];
    const all = [];
    let cur = _identity(8);
    for (let k2 = 0; k2 < 15; k2++) {
      for (const x of Xmats) all.push(_matMul(cur, x, 8));
      if (k2 < 14) cur = _matMul(cur, cox, 8);
    }
    const seen = /* @__PURE__ */ new Set();
    reps = [];
    for (const m of all) {
      const mk = _matKey(m);
      if (!seen.has(mk)) {
        seen.add(mk);
        reps.push(m);
      }
    }
  } else if (dt.series === "F") {
    reps = _buildCosetReps(C, 4, [0, 1, 2, 3], 3, [[]]);
  } else if (dt.series === "G") {
    reps = _buildCosetReps(C, 2, [0, 1], 2, [[]]);
  }
  _cosetCache.set(k, reps);
  return reps;
}
function weylloop(dt, v, action) {
  if (dt.factors !== void 0) {
    _weylloopProduct(dt, v, action);
    return;
  }
  const info = _getTypeInfo(dt);
  const { subtype, epsDim, permSize, w2e, e2w } = info;
  const rank = dt.rank;
  const cosetReps = _getCosetReps(dt);
  const tmpW = new Int32Array(rank);
  const altW = new Int32Array(rank);
  const suborbitEps = [];
  for (const rep of cosetReps) {
    for (let j = 0; j < rank; j++) {
      let s = 0;
      for (let i = 0; i < rank; i++) s += v[i] * rep[i * rank + j];
      altW[j] = s;
    }
    const eTmp = new Int32Array(epsDim);
    w2e(altW, eTmp);
    _normalform(eTmp, permSize, subtype);
    suborbitEps.push(eTmp);
  }
  suborbitEps.sort((a, b) => {
    for (let i = 0; i < a.length; i++) {
      if (a[i] !== b[i]) return a[i] - b[i];
    }
    return 0;
  });
  const uniqueEps = [suborbitEps[0]];
  for (let i = 1; i < suborbitEps.length; i++) {
    let dup = true;
    for (let j = 0; j < suborbitEps[i].length; j++) {
      if (suborbitEps[i][j] !== suborbitEps[i - 1][j]) {
        dup = false;
        break;
      }
    }
    if (!dup) uniqueEps.push(suborbitEps[i]);
  }
  const inx = new Int32Array(permSize);
  for (const eRep of uniqueEps) {
    const w = Int32Array.from(eRep);
    if (subtype === "A") {
      do {
        e2w(w, tmpW);
        action(tmpW);
      } while (_nextperm(w, permSize));
    } else {
      const alternate = subtype === "D" && w[0] !== 0;
      do {
        let nNonzero = 0;
        for (let i = 0; i < permSize; i++) {
          if (w[i] !== 0) inx[nNonzero++] = i;
        }
        let signcount = 0;
        while (true) {
          e2w(w, tmpW);
          action(tmpW);
          let ii = 0;
          let bits = signcount;
          signcount++;
          if (alternate) {
            w[0] = -w[0];
            ii++;
          }
          while (ii < nNonzero && bits & 1) {
            bits >>= 1;
            ii++;
          }
          if (ii >= nNonzero) {
            ii--;
            if (ii >= 0) w[inx[ii]] = -w[inx[ii]];
            break;
          }
          w[inx[ii]] = -w[inx[ii]];
        }
        const minus = alternate && w[0] < 0;
        if (minus) w[0] = -w[0];
        const hasNext = _nextperm(w, permSize);
        if (minus) w[0] = -w[0];
        if (!hasNext) break;
      } while (true);
    }
  }
}
function _weylloopProduct(dt, v, action) {
  const factors = dt.factors;
  const nFactors = factors.length;
  const offsets = dt.offsets;
  const totalRank = dt.rank;
  const factorOrbits = [];
  for (let fi = 0; fi < nFactors; fi++) {
    const fdt = factors[fi];
    const fr = fdt.rank;
    const off = offsets[fi];
    const vf = new Int32Array(fr);
    for (let i = 0; i < fr; i++) vf[i] = v[off + i];
    const orbits = [];
    weylloop(fdt, vf, (w) => orbits.push(Int32Array.from(w)));
    factorOrbits.push(orbits);
  }
  const tmpFull = new Int32Array(totalRank);
  const _recurse = (depth) => {
    if (depth >= nFactors) {
      action(tmpFull);
      return;
    }
    const off = offsets[depth];
    const fr = factors[depth].rank;
    for (const orb of factorOrbits[depth]) {
      for (let i = 0; i < fr; i++) tmpFull[off + i] = orb[i];
      _recurse(depth + 1);
    }
  };
  _recurse(0);
}

// ../HomogeneousTools/Lie.js/src/Characters.js
var _domCharCache = /* @__PURE__ */ new Map();
var _tensorCache = /* @__PURE__ */ new Map();
var _symCache2 = /* @__PURE__ */ new Map();
var _extCache = /* @__PURE__ */ new Map();
function dominantCharacter(dt, hw) {
  const ck = dt.key() + ":" + weightKey(hw);
  if (_domCharCache.has(ck)) return _domCharCache.get(ck);
  const n = dt.rank;
  const rs = rootSystem(dt);
  const C = cartanMatrix(dt);
  const d = cartanSymmetrizer(dt);
  const B = cartanBilinearForm(dt);
  const { scale: S, matrix: BomegaS } = omegaBilinearFormScaled(dt);
  const nPos = rs.nPositiveRoots;
  const alphaW = rs.positiveRoots.map((alpha) => rootToWeight(alpha, C, n));
  const alphaDotAlpha = new Int32Array(nPos);
  for (let k = 0; k < nPos; k++) {
    const v = rs.positiveRoots[k];
    let s = 0;
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        s += v[i] * B[i * n + j] * v[j];
    alphaDotAlpha[k] = s;
  }
  const dAlpha = rs.positiveRoots.map((v) => {
    const da = new Int32Array(n);
    for (let i = 0; i < n; i++) da[i] = d[i] * v[i];
    return da;
  });
  const rho = weylVector(n);
  const domWeights = dominantWeights(dt, hw);
  const nDom = domWeights.length;
  const domIndex = /* @__PURE__ */ new Map();
  for (let i = 0; i < nDom; i++) domIndex.set(weightKey(domWeights[i]), i);
  const domMults = new Int32Array(nDom);
  domMults[0] = 1;
  const hwRho = weightAdd(hw, rho);
  let firstTermS = 0;
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++)
      firstTermS += hwRho[i] * BomegaS[i * n + j] * hwRho[j];
  const rootMults = new Int32Array(nPos);
  for (let idx = 1; idx < nDom; idx++) {
    const mu = domWeights[idx];
    let sigma = 0;
    rootMults.fill(1);
    let hasZero = false;
    for (let j = 0; j < n; j++) if (mu[j] === 0) {
      hasZero = true;
      break;
    }
    if (hasZero) {
      const lastRoot = rs.positiveRoots[nPos - 1];
      let maxH = 0;
      for (let i = 0; i < n; i++) maxH += lastRoot[i];
      for (let lev = 1; lev <= maxH; lev++) {
        for (let j = 0; j < n; j++) {
          if (mu[j] !== 0) continue;
          for (let k = 0; k < nPos; k++) {
            if (rootMults[k] === 0) continue;
            let h = 0;
            for (let i = 0; i < n; i++) h += rs.positiveRoots[k][i];
            if (h !== lev) continue;
            const e = -alphaW[k][j];
            if (e > 0) {
              const target = rs.refl[j][k];
              if (target !== 4294967295) {
                rootMults[target] += rootMults[k];
                rootMults[k] = 0;
              }
            }
          }
        }
      }
    }
    for (let k = 0; k < nPos; k++) {
      if (rootMults[k] === 0) continue;
      let muDotAlpha = 0;
      for (let i = 0; i < n; i++) muDotAlpha += mu[i] * dAlpha[k][i];
      let nu = weightAdd(mu, alphaW[k]);
      let step = 1;
      while (true) {
        const v = Int32Array.from(nu);
        let s = 0;
        while (s < n) {
          if (v[s] < 0) {
            const p = v[s];
            for (let j = 0; j < n; j++) v[j] -= p * C[j * n + s];
            s = 0;
          } else {
            s++;
          }
        }
        const di = domIndex.get(weightKey(v));
        if (di === void 0) break;
        if (domMults[di] !== 0) {
          const ip = muDotAlpha + step * alphaDotAlpha[k];
          sigma += rootMults[k] * domMults[di] * ip;
        }
        nu = weightAdd(nu, alphaW[k]);
        step++;
      }
    }
    if (sigma !== 0) {
      const muRho = weightAdd(mu, rho);
      let secondTermS = 0;
      for (let i = 0; i < n; i++)
        for (let j = 0; j < n; j++)
          secondTermS += muRho[i] * BomegaS[i * n + j] * muRho[j];
      const denomS = firstTermS - secondTermS;
      const numerator = 2 * S * sigma;
      domMults[idx] = numerator / denomS;
    }
  }
  const result = /* @__PURE__ */ new Map();
  for (let i = 0; i < nDom; i++) {
    if (domMults[i] !== 0) {
      result.set(weightKey(domWeights[i]), { weight: domWeights[i], mult: domMults[i] });
    }
  }
  _domCharCache.set(ck, result);
  return result;
}
function _brauerKlimykDominant(dt, domChar, mu) {
  const n = dt.rank;
  const C = cartanMatrix(dt);
  const result = /* @__PURE__ */ new Map();
  for (const [, { weight: lambdaDom, mult: m }] of domChar) {
    if (m === 0) continue;
    weylloop(dt, lambdaDom, (orbitWt) => {
      const dr = new Int32Array(n);
      for (let j = 0; j < n; j++) dr[j] = mu[j] + orbitWt[j];
      let sign = 1;
      while (true) {
        let done = true;
        for (let s = 0; s < n; s++) {
          const c = dr[s];
          if (c === -1) {
            sign = 0;
            done = true;
            break;
          }
          if (c <= -2) {
            sign = -sign;
            const p = c + 1;
            for (let j = 0; j < n; j++) dr[j] -= p * C[j * n + s];
            done = false;
            break;
          }
        }
        if (done) break;
      }
      if (sign !== 0) {
        const k = weightKey(dr);
        result.set(k, (result.get(k) || 0) + sign * m);
      }
    });
  }
  for (const [k, v] of result) {
    if (v === 0) result.delete(k);
  }
  return result;
}
function tensorProduct(dt, lambda, mu) {
  const k1 = dt.key() + ":" + weightKey(lambda) + ":" + weightKey(mu);
  if (_tensorCache.has(k1)) return _tensorCache.get(k1);
  const k2 = dt.key() + ":" + weightKey(mu) + ":" + weightKey(lambda);
  if (_tensorCache.has(k2)) return _tensorCache.get(k2);
  const dL = weylDimension(dt, lambda);
  const dM = weylDimension(dt, mu);
  let result;
  if (dL <= dM) {
    result = _brauerKlimykDominant(dt, dominantCharacter(dt, lambda), mu);
  } else {
    result = _brauerKlimykDominant(dt, dominantCharacter(dt, mu), lambda);
  }
  _tensorCache.set(k1, result);
  return result;
}
function tensorProductChars(dt, charA, charB) {
  const result = /* @__PURE__ */ new Map();
  for (const [kA, mA] of charA) {
    const lambdaArr = kA.split(",").map(Number);
    const lambda = Int32Array.from(lambdaArr);
    for (const [kB, mB] of charB) {
      const muArr = kB.split(",").map(Number);
      const mu = Int32Array.from(muArr);
      const tp = tensorProduct(dt, lambda, mu);
      const coeff = mA * mB;
      for (const [wk, m] of tp) {
        result.set(wk, (result.get(wk) || 0) + coeff * m);
      }
    }
  }
  for (const [k, v] of result) if (v === 0) result.delete(k);
  return result;
}
function dualWeight(dt, hw) {
  const n = dt.rank;
  const C = cartanMatrix(dt);
  const w0 = longestElement(dt);
  const v = Int32Array.from(hw);
  for (const s of w0) {
    const p = v[s];
    for (let j = 0; j < n; j++) v[j] -= p * C[j * n + s];
  }
  for (let j = 0; j < n; j++) v[j] = -v[j];
  return v;
}
function dualChar(dt, char) {
  let isDomChar = false;
  for (const [, v] of char) {
    isDomChar = typeof v === "object" && v !== null && "mult" in v;
    break;
  }
  const result = /* @__PURE__ */ new Map();
  for (const [k, v] of char) {
    const m = isDomChar ? v.mult : v;
    const w = Int32Array.from(k.split(",").map(Number));
    const wd = dualWeight(dt, w);
    const wk = weightKey(wd);
    if (isDomChar) {
      const existing = result.get(wk);
      result.set(wk, { weight: wd, mult: (existing ? existing.mult : 0) + m });
    } else {
      result.set(wk, (result.get(wk) || 0) + m);
    }
  }
  for (const [k, v] of result) {
    const val = isDomChar ? v.mult : v;
    if (val === 0) result.delete(k);
  }
  return result;
}
function _adamsOperatorDom(dt, hw, k) {
  const domChar = dominantCharacter(dt, hw);
  const result = /* @__PURE__ */ new Map();
  for (const [, { weight, mult }] of domChar) {
    const scaled = weightScale(weight, k);
    const sk = weightKey(scaled);
    result.set(sk, { weight: scaled, mult });
  }
  return result;
}
function symmetricPower(dt, hw, k) {
  if (k < 0) return /* @__PURE__ */ new Map();
  if (k === 0) {
    const z = new Int32Array(dt.rank);
    return /* @__PURE__ */ new Map([[weightKey(z), 1]]);
  }
  if (k === 1) return /* @__PURE__ */ new Map([[weightKey(hw), 1]]);
  const ck = dt.key() + ":sym:" + weightKey(hw) + ":" + k;
  if (_symCache2.has(ck)) return _symCache2.get(ck);
  const result = /* @__PURE__ */ new Map();
  for (let r = 1; r <= k; r++) {
    const adams = _adamsOperatorDom(dt, hw, r);
    const prev = symmetricPower(dt, hw, k - r);
    for (const [pk, pm] of prev) {
      const mu = Int32Array.from(pk.split(",").map(Number));
      const bk = _brauerKlimykDominant(dt, adams, mu);
      for (const [wk, m] of bk) {
        result.set(wk, (result.get(wk) || 0) + pm * m);
      }
    }
  }
  for (const [wk, m] of result) {
    result.set(wk, m / k);
  }
  for (const [k2, v] of result) if (v === 0) result.delete(k2);
  _symCache2.set(ck, result);
  return result;
}
function exteriorPower(dt, hw, k) {
  if (k < 0) return /* @__PURE__ */ new Map();
  if (k === 0) {
    const z = new Int32Array(dt.rank);
    return /* @__PURE__ */ new Map([[weightKey(z), 1]]);
  }
  if (k === 1) return /* @__PURE__ */ new Map([[weightKey(hw), 1]]);
  const dim = Number(weylDimension(dt, hw));
  if (k > dim) return /* @__PURE__ */ new Map();
  const ck = dt.key() + ":ext:" + weightKey(hw) + ":" + k;
  if (_extCache.has(ck)) return _extCache.get(ck);
  const result = /* @__PURE__ */ new Map();
  for (let r = 1; r <= k; r++) {
    const sign = r % 2 === 0 ? -1 : 1;
    const adams = _adamsOperatorDom(dt, hw, r);
    const prev = exteriorPower(dt, hw, k - r);
    for (const [pk, pm] of prev) {
      const mu = Int32Array.from(pk.split(",").map(Number));
      const bk = _brauerKlimykDominant(dt, adams, mu);
      for (const [wk, m] of bk) {
        result.set(wk, (result.get(wk) || 0) + sign * pm * m);
      }
    }
  }
  for (const [wk, m] of result) {
    result.set(wk, m / k);
  }
  for (const [k2, v] of result) if (v === 0) result.delete(k2);
  _extCache.set(ck, result);
  return result;
}
function charAdd(a, b) {
  const r = new Map(a);
  for (const [k, m] of b) r.set(k, (r.get(k) || 0) + m);
  for (const [k, v] of r) if (v === 0) r.delete(k);
  return r;
}
function charSub(a, b) {
  const r = new Map(a);
  for (const [k, m] of b) r.set(k, (r.get(k) || 0) - m);
  for (const [k, v] of r) if (v === 0) r.delete(k);
  return r;
}
function charScale(c, s) {
  if (s === 0) return /* @__PURE__ */ new Map();
  const r = /* @__PURE__ */ new Map();
  for (const [k, m] of c) r.set(k, s * m);
  return r;
}
function charDimension(dt, char) {
  let total = 0n;
  let isDomChar = false;
  for (const [, v] of char) {
    isDomChar = typeof v === "object" && v !== null && "mult" in v;
    break;
  }
  if (isDomChar) {
    for (const [, v] of char) {
      let orbitSize = 0;
      weylloop(dt, v.weight, () => {
        orbitSize++;
      });
      total += BigInt(v.mult) * BigInt(orbitSize);
    }
  } else {
    for (const [k, m] of char) {
      const w = Int32Array.from(k.split(",").map(Number));
      total += BigInt(m) * weylDimension(dt, w);
    }
  }
  return total;
}
function irrepChar(hw) {
  return /* @__PURE__ */ new Map([[weightKey(hw), 1]]);
}
function adamsOperator(dt, hw, k) {
  const domChar = dominantCharacter(dt, hw);
  const n = dt.rank;
  const C = cartanMatrix(dt);
  const totalCounts = /* @__PURE__ */ new Map();
  for (const [, { weight, mult }] of domChar) {
    weylloop(dt, weight, (orbitWt) => {
      const scaledOrb = weightScale(orbitWt, k);
      const dom = conjugateDominantWeight(scaledOrb, C, n);
      const dk = weightKey(dom);
      if (!totalCounts.has(dk)) {
        totalCounts.set(dk, { weight: dom, count: 0 });
      }
      totalCounts.get(dk).count += mult;
    });
  }
  const result = /* @__PURE__ */ new Map();
  for (const [key, { weight, count }] of totalCounts) {
    let orbitSize = 0;
    weylloop(dt, weight, () => {
      orbitSize++;
    });
    const mult = count / orbitSize;
    if (mult !== 0) {
      result.set(key, { weight, mult });
    }
  }
  return result;
}
function isEffective(char) {
  for (const [, m] of char) if (m < 0) return false;
  return true;
}
function isIrreducible(char) {
  if (char.size !== 1) return false;
  for (const [, m] of char) if (m !== 1) return false;
  return true;
}
function highestWeight(dt, char) {
  let best = null;
  let bestLevel = -Infinity;
  for (const [k] of char) {
    const w = Int32Array.from(k.split(",").map(Number));
    let level = 0;
    for (let i = 0; i < w.length; i++) level += w[i];
    if (level > bestLevel) {
      bestLevel = level;
      best = w;
    }
  }
  return best;
}
function weylOrbit(dt, w) {
  const orbit = [];
  weylloop(dt, w, (orbitWt) => {
    orbit.push(Int32Array.from(orbitWt));
  });
  return orbit;
}
function applyW0(dt, w) {
  const n = dt.rank;
  const C = cartanMatrix(dt);
  const w0 = longestElement(dt);
  const v = Int32Array.from(w);
  for (const s of w0) {
    const p = v[s];
    for (let j = 0; j < n; j++) v[j] -= p * C[j * n + s];
  }
  return v;
}
function lrTensorProduct(dt, lambda, mu) {
  if (dt.series !== "A") {
    throw new Error("lrTensorProduct is only implemented for type A");
  }
  return tensorProduct(dt, lambda, mu);
}
function plethysm(dt, partition, hw) {
  if (partition.length === 0) {
    const z = new Int32Array(dt.rank);
    return /* @__PURE__ */ new Map([[weightKey(z), 1]]);
  }
  if (partition.length === 1 && partition[0] === 1) {
    return /* @__PURE__ */ new Map([[weightKey(hw), 1]]);
  }
  if (partition.length === 1) {
    return symmetricPower(dt, hw, partition[0]);
  }
  if (partition.every((p) => p === 1)) {
    return exteriorPower(dt, hw, partition.length);
  }
  const r = partition.length;
  const maxK = Math.max(...partition.map((p, i) => p - i + r - 1));
  const syms = /* @__PURE__ */ new Map();
  for (let k = 0; k <= maxK; k++) {
    syms.set(k, symmetricPower(dt, hw, k));
  }
  if (r === 1) return syms.get(partition[0]) || /* @__PURE__ */ new Map();
  if (r === 2) {
    const a = partition[0];
    const b = partition[1];
    const S2 = (k) => k >= 0 ? syms.get(k) || symmetricPower(dt, hw, k) : /* @__PURE__ */ new Map();
    const AD = tensorProductChars(dt, S2(a), S2(b));
    const BC = tensorProductChars(dt, S2(a + 1), S2(b - 1));
    return charSub(AD, BC);
  }
  function _subPlethysm(subPartition) {
    return plethysm(dt, subPartition, hw);
  }
  const M = [];
  for (let i = 0; i < r; i++) {
    M.push([]);
    for (let j = 0; j < r; j++) {
      const k = partition[i] - i + j;
      M[i].push(k >= 0 ? syms.get(k) || symmetricPower(dt, hw, k) : null);
    }
  }
  function det(indices2) {
    if (indices2.length === 1) {
      const k = partition[indices2[0]] - indices2[0] + 0;
      return k >= 0 ? syms.get(k) || symmetricPower(dt, hw, k) : /* @__PURE__ */ new Map();
    }
    let result2 = /* @__PURE__ */ new Map();
    for (let jj = 0; jj < indices2.length; jj++) {
      const i = indices2[0];
      const k = partition[i] - i + jj;
      if (k < 0) continue;
      const entry = syms.get(k) || symmetricPower(dt, hw, k);
      if (entry.size === 0) continue;
      const subIndices = indices2.slice(1);
      const subMat = subIndices;
    }
    return result2;
  }
  const charTensor = tensorProductChars(
    dt,
    symmetricPower(dt, hw, partition[0]),
    plethysm(dt, partition.slice(1), hw)
  );
  if (r === 3) {
    const S2 = (k) => k >= 0 ? syms.get(k) || symmetricPower(dt, hw, k) : /* @__PURE__ */ new Map();
    const m00 = S2(partition[0]), m01 = S2(partition[0] + 1), m02 = S2(partition[0] + 2);
    const m10 = S2(partition[1] - 1), m11 = S2(partition[1]), m12 = S2(partition[1] + 1);
    const m20 = S2(partition[2] - 2), m21 = S2(partition[2] - 1), m22 = S2(partition[2]);
    const t = (a, b) => a.size > 0 && b.size > 0 ? tensorProductChars(dt, a, b) : /* @__PURE__ */ new Map();
    const minor00 = charSub(t(m11, m22), t(m12, m21));
    const minor01 = charSub(t(m10, m22), t(m12, m20));
    const minor02 = charSub(t(m10, m21), t(m11, m20));
    const term0 = t(m00, minor00);
    const term1 = t(m01, minor01);
    const term2 = t(m02, minor02);
    return charSub(charAdd(term0, term2), term1);
  }
  function perms(arr) {
    if (arr.length <= 1) return [arr];
    const result2 = [];
    for (let i = 0; i < arr.length; i++) {
      const rest = [...arr.slice(0, i), ...arr.slice(i + 1)];
      for (const p of perms(rest)) result2.push([arr[i], ...p]);
    }
    return result2;
  }
  function permSign(perm) {
    let inv = 0;
    for (let i = 0; i < perm.length; i++)
      for (let j = i + 1; j < perm.length; j++)
        if (perm[i] > perm[j]) inv++;
    return inv % 2 === 0 ? 1 : -1;
  }
  const indices = Array.from({ length: r }, (_, i) => i);
  const allPerms = perms(indices);
  const S = (k) => k >= 0 ? syms.get(k) || symmetricPower(dt, hw, k) : /* @__PURE__ */ new Map();
  let result = /* @__PURE__ */ new Map();
  for (const perm of allPerms) {
    const sign = permSign(perm);
    let prod = irrepChar(new Int32Array(dt.rank));
    prod = /* @__PURE__ */ new Map([[weightKey(new Int32Array(dt.rank)), 1]]);
    let zero = false;
    for (let i = 0; i < r; i++) {
      const k = partition[i] - i + perm[i];
      const factor = S(k);
      if (factor.size === 0) {
        zero = true;
        break;
      }
      prod = tensorProductChars(dt, prod, factor);
      if (prod.size === 0) {
        zero = true;
        break;
      }
    }
    if (zero) continue;
    if (sign === 1) {
      result = charAdd(result, prod);
    } else {
      result = charSub(result, prod);
    }
  }
  for (const [k, v] of result) if (v === 0) result.delete(k);
  return result;
}

// ../HomogeneousTools/Lie.js/src/MarkedDynkinType.js
var _leviTypeCache = /* @__PURE__ */ new Map();
var _leviPermCache = /* @__PURE__ */ new Map();
var _decompMatrixCache = /* @__PURE__ */ new Map();
var _decompInvCache = /* @__PURE__ */ new Map();
var _anticanCache = /* @__PURE__ */ new Map();
var _tangentWtCache = /* @__PURE__ */ new Map();
var _centralScaleCache = /* @__PURE__ */ new Map();
var MarkedDynkinType = class {
  /**
   * @param {SimpleDynkinType | ProductDynkinType} dt  Ambient Dynkin type.
   * @param {number[]} marked  Sorted 1-based indices of marked (crossed-out) nodes.
   */
  constructor(dt, marked) {
    if (typeof dt === "string") dt = parseDynkinType(dt);
    const sorted = [...marked].sort((a, b) => a - b);
    const R = dt.rank;
    for (const m of sorted) {
      if (m < 1 || m > R) throw new RangeError(`Marked node ${m} out of range [1, ${R}]`);
    }
    if (new Set(sorted).size !== sorted.length) throw new Error("Marked nodes must be distinct");
    this.dt = dt;
    this.marked = sorted;
    this.unmarked = [];
    for (let i = 1; i <= R; i++) {
      if (!sorted.includes(i)) this.unmarked.push(i);
    }
    this._key = dt.key() + "/{" + sorted.join(",") + "}";
  }
  get rank() {
    return this.dt.rank;
  }
  get centralRank() {
    return this.marked.length;
  }
  get leviRank() {
    return this.unmarked.length;
  }
  key() {
    return this._key;
  }
  toString() {
    return `${this.dt} / P_{${this.marked.join(",")}}`;
  }
};
function _identifyCartanSubMatrix(Csub, n) {
  if (n === 0) return { type: null, perm: [] };
  const adj = new Array(n).fill(null).map(() => []);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (Csub[i * n + j] !== 0) {
        adj[i].push(j);
        adj[j].push(i);
      }
    }
  }
  const visited = new Uint8Array(n);
  const components = [];
  for (let start = 0; start < n; start++) {
    if (visited[start]) continue;
    const comp = [];
    const stack = [start];
    visited[start] = 1;
    while (stack.length > 0) {
      const v = stack.pop();
      comp.push(v);
      for (const u of adj[v]) {
        if (!visited[u]) {
          visited[u] = 1;
          stack.push(u);
        }
      }
    }
    comp.sort((a, b) => a - b);
    components.push(comp);
  }
  const factors = [];
  const globalPerm = new Array(n);
  let offset = 0;
  for (const comp of components) {
    const m = comp.length;
    const Cc = new Array(m * m);
    for (let i = 0; i < m; i++)
      for (let j = 0; j < m; j++)
        Cc[i * m + j] = Csub[comp[i] * n + comp[j]];
    const { type: sdt, perm: sperm } = _identifySimpleCartan(Cc, m);
    factors.push(sdt);
    for (let i = 0; i < m; i++) {
      globalPerm[offset + i] = comp[sperm[i]];
    }
    offset += m;
  }
  if (factors.length === 1) {
    return { type: factors[0], perm: globalPerm };
  }
  return { type: new ProductDynkinType(factors), perm: globalPerm };
}
function _identifySimpleCartan(C, n) {
  const id = Array.from({ length: n }, (_, i) => i);
  if (n === 1) return { type: new SimpleDynkinType("A", 1), perm: id };
  const hasDouble = [];
  const hasTriple = [];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i === j) continue;
      if (C[i * n + j] === -2) hasDouble.push([i, j]);
      if (C[i * n + j] === -3) hasTriple.push([i, j]);
    }
  }
  if (n === 2 && hasTriple.length > 0) {
    const [i, j] = hasTriple[0];
    return { type: new SimpleDynkinType("G", 2), perm: [i, j] };
  }
  if (hasDouble.length > 0 && n === 4) {
    const [a, b] = hasDouble[0];
    const degs2 = new Array(n).fill(0);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        if (i !== j && C[i * n + j] !== 0) degs2[i]++;
    if (Math.max(...degs2) <= 2) {
      const chain = _linearChain(n, C);
      if (chain) {
        let isF4 = false;
        for (let k = 0; k < n - 1; k++) {
          if (C[chain[k + 1] * n + chain[k]] === -2 || C[chain[k] * n + chain[k + 1]] === -2) {
            if (C[chain[k + 1] * n + chain[k]] === -2 && k === 1) {
              isF4 = true;
              return { type: new SimpleDynkinType("F", 4), perm: chain };
            } else if (C[chain[k] * n + chain[k + 1]] === -2 && k === 1) {
              isF4 = true;
              const rev = [...chain].reverse();
              return { type: new SimpleDynkinType("F", 4), perm: rev };
            }
            if (!isF4) {
              const rev = [...chain].reverse();
              if (C[rev[2] * n + rev[1]] === -2) {
                return { type: new SimpleDynkinType("F", 4), perm: rev };
              }
            }
          }
        }
      }
    }
  }
  if (hasDouble.length > 0 && hasTriple.length === 0) {
    const chain = _linearChain(n, C);
    if (chain) {
      if (C[chain[n - 1] * n + chain[n - 2]] === -2) {
        return { type: new SimpleDynkinType("B", n), perm: chain };
      }
      if (C[chain[n - 2] * n + chain[n - 1]] === -2) {
        return { type: new SimpleDynkinType("C", n), perm: chain };
      }
      const rev = [...chain].reverse();
      if (C[rev[n - 1] * n + rev[n - 2]] === -2) {
        return { type: new SimpleDynkinType("B", n), perm: rev };
      }
      if (C[rev[n - 2] * n + rev[n - 1]] === -2) {
        return { type: new SimpleDynkinType("C", n), perm: rev };
      }
    }
  }
  const degs = new Array(n).fill(0);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++)
      if (i !== j && C[i * n + j] !== 0) degs[i]++;
  const maxDeg = Math.max(...degs);
  if (maxDeg <= 2) {
    return { type: new SimpleDynkinType("A", n), perm: _linearChain(n, C) };
  }
  const branchNode = degs.findIndex((d) => d === 3);
  if (branchNode < 0) throw new Error(`Cannot identify Cartan matrix of rank ${n}`);
  const nbrs = [];
  for (let j = 0; j < n; j++) {
    if (j !== branchNode && C[branchNode * n + j] !== 0) nbrs.push(j);
  }
  const arms = nbrs.map((start) => {
    const arm = [start];
    let prev = branchNode;
    let cur = start;
    while (true) {
      let next = -1;
      for (let j = 0; j < n; j++) {
        if (j !== prev && j !== cur && C[cur * n + j] !== 0) {
          next = j;
          break;
        }
      }
      if (next < 0) break;
      arm.push(next);
      prev = cur;
      cur = next;
    }
    return arm;
  });
  arms.sort((a, b) => a.length - b.length);
  if (n >= 4 && arms.length === 3 && arms[0].length === 1 && arms[1].length === 1) {
    const longArm = arms[2];
    const shortArms = [arms[0], arms[1]];
    const perm2 = [];
    for (let i = longArm.length - 1; i >= 0; i--) perm2.push(longArm[i]);
    perm2.push(branchNode);
    perm2.push(shortArms[0][0]);
    perm2.push(shortArms[1][0]);
    return { type: new SimpleDynkinType("D", n), perm: perm2 };
  }
  const shortArmIdx = arms.findIndex((a) => a.length === 1);
  if (shortArmIdx < 0) throw new Error("Cannot identify E-type");
  const otherArms = arms.filter((_, i) => i !== shortArmIdx);
  otherArms.sort((a, b) => a.length - b.length);
  const shorterArm = otherArms[0];
  const longerArm = otherArms[1];
  const branchArm = arms[shortArmIdx];
  const perm = [];
  for (let i = shorterArm.length - 1; i >= 0; i--) perm.push(shorterArm[i]);
  perm.splice(1, 0, branchArm[0]);
  perm.push(branchNode);
  for (const nd of longerArm) perm.push(nd);
  const permE = new Array(n);
  permE[3] = branchNode;
  permE[1] = branchArm[0];
  permE[2] = shorterArm[0];
  if (shorterArm.length >= 2) permE[0] = shorterArm[1];
  else permE[0] = shorterArm[0];
  for (let i = 0; i < longerArm.length; i++) {
    permE[4 + i] = longerArm[i];
  }
  return { type: new SimpleDynkinType("E", n), perm: permE };
}
function _linearChain(n, C) {
  const degs = new Array(n).fill(0);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++)
      if (i !== j && C[i * n + j] !== 0) degs[i]++;
  let start = degs.findIndex((d) => d === 1);
  if (start < 0) start = 0;
  const chain = [start];
  const visited = new Uint8Array(n);
  visited[start] = 1;
  let cur = start;
  for (let step = 1; step < n; step++) {
    for (let j = 0; j < n; j++) {
      if (!visited[j] && C[cur * n + j] !== 0) {
        chain.push(j);
        visited[j] = 1;
        cur = j;
        break;
      }
    }
  }
  return chain;
}
function leviType(mdt) {
  const k = mdt.key();
  if (_leviTypeCache.has(k)) return _leviTypeCache.get(k);
  const unmarked = mdt.unmarked;
  const n = unmarked.length;
  if (n === 0) {
    _leviTypeCache.set(k, null);
    return null;
  }
  const R = mdt.rank;
  const C = cartanMatrix(mdt.dt);
  const Csub = new Array(n * n);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++)
      Csub[i * n + j] = C[(unmarked[i] - 1) * R + (unmarked[j] - 1)];
  const { type: lt } = _identifyCartanSubMatrix(Csub, n);
  _leviTypeCache.set(k, lt);
  return lt;
}
function leviPermutation(mdt) {
  const k = mdt.key();
  if (_leviPermCache.has(k)) return _leviPermCache.get(k);
  const unmarked = mdt.unmarked;
  const n = unmarked.length;
  if (n === 0) {
    _leviPermCache.set(k, []);
    return [];
  }
  const R = mdt.rank;
  const C = cartanMatrix(mdt.dt);
  const Csub = new Array(n * n);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++)
      Csub[i * n + j] = C[(unmarked[i] - 1) * R + (unmarked[j] - 1)];
  const { perm } = _identifyCartanSubMatrix(Csub, n);
  _leviPermCache.set(k, perm);
  return perm;
}
function decompositionMatrix(mdt) {
  const k = mdt.key();
  if (_decompMatrixCache.has(k)) return _decompMatrixCache.get(k);
  const R = mdt.rank;
  const Cinv = cartanMatrixInverse(mdt.dt);
  const M = new Array(R * R);
  for (let i = 0; i < R * R; i++) M[i] = rat(0);
  for (const u of mdt.unmarked) {
    const i = u - 1;
    M[i * R + i] = rat(1);
  }
  for (const m of mdt.marked) {
    const j = m - 1;
    for (let col = 0; col < R; col++) {
      M[j * R + col] = Cinv[j * R + col];
    }
  }
  _decompMatrixCache.set(k, M);
  return M;
}
function decompositionMatrixInv(mdt) {
  const k = mdt.key();
  if (_decompInvCache.has(k)) return _decompInvCache.get(k);
  const R = mdt.rank;
  const M = decompositionMatrix(mdt);
  const aug = new Array(R);
  for (let i = 0; i < R; i++) {
    aug[i] = new Array(2 * R);
    for (let j = 0; j < R; j++) aug[i][j] = M[i * R + j];
    for (let j = 0; j < R; j++) aug[i][R + j] = i === j ? rat(1) : rat(0);
  }
  for (let col = 0; col < R; col++) {
    let pivot = -1;
    for (let row = col; row < R; row++) {
      if (!aug[row][col].isZero()) {
        pivot = row;
        break;
      }
    }
    if (pivot < 0) throw new Error("Decomposition matrix is singular");
    if (pivot !== col) [aug[col], aug[pivot]] = [aug[pivot], aug[col]];
    const scale = aug[col][col];
    for (let j = 0; j < 2 * R; j++) aug[col][j] = aug[col][j].div(scale);
    for (let row = 0; row < R; row++) {
      if (row === col) continue;
      const f = aug[row][col];
      if (f.isZero()) continue;
      for (let j = 0; j < 2 * R; j++) aug[row][j] = aug[row][j].sub(f.mul(aug[col][j]));
    }
  }
  const inv = new Array(R * R);
  for (let i = 0; i < R; i++)
    for (let j = 0; j < R; j++)
      inv[i * R + j] = aug[i][R + j];
  _decompInvCache.set(k, inv);
  return inv;
}
function centralScalingFactor(mdt) {
  const k = mdt.key();
  if (_centralScaleCache.has(k)) return _centralScaleCache.get(k);
  const R = mdt.rank;
  const Cinv = cartanMatrixInverse(mdt.dt);
  let sf = 1;
  for (const m of mdt.marked) {
    const j = m - 1;
    for (let col = 0; col < R; col++) {
      sf = lcm(sf, Cinv[j * R + col].q);
    }
  }
  _centralScaleCache.set(k, sf);
  return sf;
}
function anticanonicalDegrees(mdt) {
  const k = mdt.key();
  if (_anticanCache.has(k)) return _anticanCache.get(k);
  const R = mdt.rank;
  const marked = mdt.marked;
  const unmarked = mdt.unmarked;
  const K = marked.length;
  const n = unmarked.length;
  const C = cartanMatrix(mdt.dt);
  if (n === 0) {
    const degs2 = new Array(K).fill(2);
    _anticanCache.set(k, degs2);
    return degs2;
  }
  const CL = new Array(n * n);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++)
      CL[i * n + j] = rat(C[(unmarked[i] - 1) * R + (unmarked[j] - 1)]);
  const ones = new Array(n).fill(rat(1));
  const aug = new Array(n);
  for (let i = 0; i < n; i++) {
    aug[i] = new Array(n + 1);
    for (let j = 0; j < n; j++) aug[i][j] = CL[i * n + j];
    aug[i][n] = ones[i];
  }
  for (let col = 0; col < n; col++) {
    let pivot = -1;
    for (let row = col; row < n; row++) {
      if (!aug[row][col].isZero()) {
        pivot = row;
        break;
      }
    }
    if (pivot !== col) [aug[col], aug[pivot]] = [aug[pivot], aug[col]];
    const s = aug[col][col];
    for (let j = 0; j <= n; j++) aug[col][j] = aug[col][j].div(s);
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const f = aug[row][col];
      if (f.isZero()) continue;
      for (let j = 0; j <= n; j++) aug[row][j] = aug[row][j].sub(f.mul(aug[col][j]));
    }
  }
  const xL = aug.map((r) => r[n]);
  const degs = new Array(K);
  for (let idx = 0; idx < K; idx++) {
    const i = marked[idx] - 1;
    let val = rat(2);
    for (let q = 0; q < n; q++) {
      val = val.sub(rat(2 * C[i * R + (unmarked[q] - 1)]).mul(xL[q]));
    }
    degs[idx] = val.toNumber();
  }
  _anticanCache.set(k, degs);
  return degs;
}
function positiveNonparabolicRoots(mdt) {
  const rs = rootSystem(mdt.dt);
  const markedSet = new Set(mdt.marked.map((m) => m - 1));
  return rs.positiveRoots.filter((alpha) => {
    for (const m of markedSet) if (alpha[m] !== 0) return true;
    return false;
  });
}
function positiveParabolicRoots(mdt) {
  const rs = rootSystem(mdt.dt);
  const markedSet = new Set(mdt.marked.map((m) => m - 1));
  return rs.positiveRoots.filter((alpha) => {
    for (const m of markedSet) if (alpha[m] !== 0) return false;
    return true;
  });
}
function tangentWeights(mdt) {
  const k = mdt.key();
  if (_tangentWtCache.has(k)) return _tangentWtCache.get(k);
  const R = mdt.rank;
  const C = cartanMatrix(mdt.dt);
  const rs = rootSystem(mdt.dt);
  const nonpar = positiveNonparabolicRoots(mdt);
  const simpleParRoots = mdt.unmarked.map((u) => {
    const v = new Int32Array(R);
    v[u - 1] = 1;
    return v;
  });
  const nonparSet = new Set(nonpar.map((alpha) => weightKey(alpha)));
  const maxRoots = nonpar.filter((alpha) => {
    for (const beta of simpleParRoots) {
      const sum = new Int32Array(R);
      for (let i = 0; i < R; i++) sum[i] = alpha[i] + beta[i];
      if (nonparSet.has(weightKey(sum))) return false;
    }
    return true;
  });
  const result = maxRoots.map((alpha) => rootToWeight(alpha, C, R));
  _tangentWtCache.set(k, result);
  return result;
}
function isGeneralizedGrassmannian(mdt) {
  return mdt.marked.length === 1;
}
function isFullFlag(mdt) {
  return mdt.marked.length === mdt.rank;
}
function isCominuscule(mdt) {
  const dt = mdt.dt;
  if (dt.factors) return false;
  if (mdt.marked.length !== 1) return false;
  const m = mdt.marked[0];
  const R = dt.rank;
  switch (dt.series) {
    case "A":
      return true;
    case "B":
      return m === 1;
    case "C":
      return m === R;
    case "D":
      return m === 1 || m === R - 1 || m === R;
    case "E":
      return R === 6 && (m === 1 || m === 6) || R === 7 && m === 7;
    default:
      return false;
  }
}
function isMinuscule(mdt) {
  const dt = mdt.dt;
  if (dt.factors) return false;
  if (mdt.marked.length !== 1) return false;
  const m = mdt.marked[0];
  const R = dt.rank;
  switch (dt.series) {
    case "A":
      return true;
    case "B":
      return m === R;
    case "C":
      return m === 1;
    case "D":
      return m === 1 || m === R - 1 || m === R;
    case "E":
      return R === 6 && (m === 1 || m === 6) || R === 7 && m === 7;
    default:
      return false;
  }
}
function isAdjoint(mdt) {
  const dt = mdt.dt;
  if (dt.factors) return false;
  const marked = mdt.marked;
  const R = dt.rank;
  if (dt.series === "A") return marked.length === 2 && marked[0] === 1 && marked[1] === R;
  if (marked.length !== 1) return false;
  const m = marked[0];
  switch (dt.series) {
    case "B":
      return m === 2;
    case "C":
      return m === 1;
    case "D":
      return R >= 4 && m === 2;
    case "E":
      return R === 6 && m === 2 || R === 7 && m === 1 || R === 8 && m === 8;
    case "F":
      return m === 1;
    case "G":
      return m === 2;
    default:
      return false;
  }
}
function isCoadjoint(mdt) {
  const dt = mdt.dt;
  if (dt.factors) return false;
  const marked = mdt.marked;
  const R = dt.rank;
  if (dt.series === "A") return marked.length === 2 && marked[0] === 1 && marked[1] === R;
  if (marked.length !== 1) return false;
  const m = marked[0];
  switch (dt.series) {
    case "B":
      return m === 1;
    case "C":
      return m === 2;
    case "D":
      return R >= 4 && m === 2;
    case "E":
      return R === 6 && m === 2 || R === 7 && m === 1 || R === 8 && m === 8;
    case "F":
      return m === 4;
    case "G":
      return m === 1;
    default:
      return false;
  }
}

// ../HomogeneousTools/Lie.js/src/PartialFlagVariety.js
var _bettiCache = /* @__PURE__ */ new Map();
var PartialFlagVariety = class {
  /**
   * @param {MarkedDynkinType} mdt
   * @param {string} [name]  Optional display name.
   */
  constructor(mdt, name) {
    if (typeof mdt === "string") {
      throw new Error("Use partialFlagVariety(dt, marked) factory");
    }
    this.mdt = mdt;
    this.name = name || mdt.toString();
  }
  get rank() {
    return this.mdt.rank;
  }
  get markedNodes() {
    return this.mdt.marked;
  }
  get unmarkedNodes() {
    return this.mdt.unmarked;
  }
  get centralRank() {
    return this.mdt.centralRank;
  }
  get picardRank() {
    return this.mdt.centralRank;
  }
  key() {
    return this.mdt.key();
  }
  toString() {
    return this.name;
  }
};
function partialFlagVariety(dt, marked) {
  if (typeof dt === "string") dt = parseDynkinType(dt);
  const mdt = new MarkedDynkinType(dt, marked);
  return new PartialFlagVariety(mdt);
}
function fullFlagVariety(dt) {
  if (typeof dt === "string") dt = parseDynkinType(dt);
  const marked = Array.from({ length: dt.rank }, (_, i) => i + 1);
  return partialFlagVariety(dt, marked);
}
function dimension(X) {
  return positiveNonparabolicRoots(X.mdt).length;
}
function eulerCharacteristic(X) {
  const wG = weylOrder(X.mdt.dt);
  const lt = leviType(X.mdt);
  const wL = lt ? weylOrder(lt) : 1n;
  return wG / wL;
}
function fanoIndex(X) {
  const degs = anticanonicalDegrees(X.mdt);
  if (degs.length === 1) return degs[0];
  let g = degs[0];
  for (let i = 1; i < degs.length; i++) {
    g = _gcdNum(g, degs[i]);
  }
  return g;
}
function _gcdNum(a, b) {
  a = Math.abs(a);
  b = Math.abs(b);
  while (b) {
    [a, b] = [b, a % b];
  }
  return a;
}
function bettiNumbers(X) {
  const k = X.key();
  if (_bettiCache.has(k)) return _bettiCache.get(k);
  const dt = X.mdt.dt;
  const lt = leviType(X.mdt);
  const degsG = _allFundamentalInvariantDegrees(dt);
  const degsL = lt ? _allFundamentalInvariantDegrees(lt) : [];
  let poly = [1];
  for (const d of degsG) {
    const factor = new Array(d).fill(1);
    poly = _polyMul(poly, factor);
  }
  for (const d of degsL) {
    poly = _polyDiv(poly, new Array(d).fill(1));
  }
  _bettiCache.set(k, poly);
  return poly;
}
function hilbertSeries(X, maxDegree = 20) {
  const mdt = X.mdt;
  const marked = mdt.marked;
  const R = mdt.rank;
  const dt = mdt.dt;
  if (marked.length !== 1) {
    throw new Error("hilbertSeries requires a generalized Grassmannian (single marked node)");
  }
  const result = [];
  const m = marked[0] - 1;
  for (let t = 0; t <= maxDegree; t++) {
    const w = new Int32Array(R);
    w[m] = t;
    const bwb = borelWeylBott(dt, w);
    if (bwb === null) {
      result.push(0n);
    } else {
      const sign = bwb.degree % 2 === 0 ? 1n : -1n;
      result.push(sign * weylDimension(dt, bwb.weight));
    }
  }
  return result;
}
function _allFundamentalInvariantDegrees(dt) {
  if (dt.factors) {
    const all = [];
    for (const f of dt.factors) all.push(...degreesFundamentalInvariants(f));
    return all;
  }
  return [...degreesFundamentalInvariants(dt)];
}
function _polyMul(a, b) {
  const lenA = a.length, lenB = b.length;
  const result = new Array(lenA + lenB - 1).fill(0);
  for (let i = 0; i < lenA; i++)
    for (let j = 0; j < lenB; j++)
      result[i + j] += a[i] * b[j];
  return result;
}
function _polyDiv(num, den) {
  const result = new Array(num.length - den.length + 1);
  const work = [...num];
  for (let i = result.length - 1; i >= 0; i--) {
    const c = work[i + den.length - 1] / den[den.length - 1];
    result[i] = c;
    for (let j = 0; j < den.length; j++) {
      work[i + j] -= c * den[j];
    }
  }
  return result;
}

// ../HomogeneousTools/Lie.js/src/IrrepLevi.js
var IrrepLevi = class {
  /**
   * @param {MarkedDynkinType} mdt   Marked Dynkin type.
   * @param {Int32Array} lambda       P-dominant weight (ambient ω-basis, length = R).
   * @param {Int32Array} central      Scaled central character (length = K = # marked).
   * @param {Int32Array} semisimple   Highest weight of semisimple part (Levi ω-basis, length = levi rank).
   */
  constructor(mdt, lambda, central, semisimple) {
    this.mdt = mdt;
    this.lambda = lambda;
    this.central = central;
    this.semisimple = semisimple;
  }
  key() {
    return this.mdt.key() + ":" + weightKey(this.lambda);
  }
};
function irrepLevi(mdt, lambda) {
  const R = mdt.rank;
  const marked = mdt.marked;
  const unmarked = mdt.unmarked;
  const K = marked.length;
  const lt = leviType(mdt);
  const lr = lt ? lt.rank : 0;
  const sf = centralScalingFactor(mdt);
  const Cinv = cartanMatrixInverse(mdt.dt);
  const central = new Int32Array(K);
  for (let idx = 0; idx < K; idx++) {
    const j = marked[idx] - 1;
    let val = 0;
    for (let col = 0; col < R; col++) {
      val += Math.round(Cinv[j * R + col].p * sf / Cinv[j * R + col].q) * lambda[col];
    }
    central[idx] = val;
  }
  const semisimple = new Int32Array(lr);
  if (lr > 0) {
    const perm = leviPermutation(mdt);
    for (let i = 0; i < lr; i++) {
      semisimple[i] = lambda[unmarked[perm[i]] - 1];
    }
  }
  return new IrrepLevi(mdt, Int32Array.from(lambda), central, semisimple);
}
function irrepLeviFromParts(mdt, central, semisimple) {
  const R = mdt.rank;
  const marked = mdt.marked;
  const unmarked = mdt.unmarked;
  const K = marked.length;
  const lr = unmarked.length;
  const sf = centralScalingFactor(mdt);
  const Minv = decompositionMatrixInv(mdt);
  const coords = new Array(R).fill(rat(0));
  for (let idx = 0; idx < K; idx++) {
    coords[marked[idx] - 1] = rat(central[idx], sf);
  }
  if (lr > 0) {
    const perm = leviPermutation(mdt);
    const invPerm = new Array(lr);
    for (let i = 0; i < lr; i++) invPerm[perm[i]] = i;
    for (let i = 0; i < lr; i++) {
      coords[unmarked[i] - 1] = rat(semisimple[invPerm[i]]);
    }
  }
  const lambda = new Int32Array(R);
  for (let i = 0; i < R; i++) {
    let val = rat(0);
    for (let j = 0; j < R; j++) {
      val = val.add(Minv[i * R + j].mul(coords[j]));
    }
    lambda[i] = Math.round(val.toNumber());
  }
  return new IrrepLevi(mdt, lambda, Int32Array.from(central), Int32Array.from(semisimple));
}
function centralPart(rep) {
  const sf = centralScalingFactor(rep.mdt);
  return Array.from(rep.central, (c) => rat(c, sf));
}
function fiberDimension(rep) {
  const lt = leviType(rep.mdt);
  if (!lt) return 1n;
  const ss = rep.semisimple;
  let allZero = true;
  for (let i = 0; i < ss.length; i++) if (ss[i] !== 0) {
    allZero = false;
    break;
  }
  if (allZero) return 1n;
  for (let i = 0; i < ss.length; i++) if (ss[i] < 0) return 0n;
  return weylDimension(lt, ss);
}
function dualIrrepLevi(rep) {
  const mdt = rep.mdt;
  const lt = leviType(mdt);
  const newCentral = Int32Array.from(rep.central, (c) => -c);
  let newSS;
  if (!lt || rep.semisimple.every((c) => c === 0)) {
    newSS = Int32Array.from(rep.semisimple);
  } else {
    newSS = dualWeight(lt, rep.semisimple);
  }
  return irrepLeviFromParts(mdt, newCentral, newSS);
}
function tensorProductIrrepLevi(a, b) {
  const mdt = a.mdt;
  const lt = leviType(mdt);
  const newCentral = new Int32Array(a.central.length);
  for (let i = 0; i < newCentral.length; i++) newCentral[i] = a.central[i] + b.central[i];
  if (!lt || a.semisimple.every((c) => c === 0) && b.semisimple.every((c) => c === 0)) {
    return [irrepLeviFromParts(mdt, newCentral, a.semisimple)];
  }
  if (a.semisimple.every((c) => c === 0)) {
    return [irrepLeviFromParts(mdt, newCentral, b.semisimple)];
  }
  if (b.semisimple.every((c) => c === 0)) {
    return [irrepLeviFromParts(mdt, newCentral, a.semisimple)];
  }
  const chi = tensorProduct(lt, a.semisimple, b.semisimple);
  const result = [];
  for (const [wk, mult] of chi) {
    const hw = Int32Array.from(wk.split(",").map(Number));
    for (let i = 0; i < mult; i++) {
      result.push(irrepLeviFromParts(mdt, newCentral, hw));
    }
  }
  return result;
}
function exteriorPowerIrrepLevi(rep, k) {
  if (k < 0) return [];
  const mdt = rep.mdt;
  const lt = leviType(mdt);
  const K = rep.central.length;
  if (k === 0) {
    return [irrepLeviFromParts(mdt, new Int32Array(K), new Int32Array(lt ? lt.rank : 0))];
  }
  if (k === 1) return [rep];
  const newCentral = Int32Array.from(rep.central, (c) => c * k);
  if (!lt || rep.semisimple.every((c) => c === 0)) {
    return [];
  }
  const dimSS = Number(weylDimension(lt, rep.semisimple));
  if (k > dimSS) return [];
  const chi = exteriorPower(lt, rep.semisimple, k);
  const result = [];
  for (const [wk, mult] of chi) {
    const hw = Int32Array.from(wk.split(",").map(Number));
    for (let i = 0; i < mult; i++) {
      result.push(irrepLeviFromParts(mdt, newCentral, hw));
    }
  }
  return result;
}
function symmetricPowerIrrepLevi(rep, k) {
  if (k < 0) return [];
  const mdt = rep.mdt;
  const lt = leviType(mdt);
  const K = rep.central.length;
  if (k === 0) {
    return [irrepLeviFromParts(mdt, new Int32Array(K), new Int32Array(lt ? lt.rank : 0))];
  }
  if (k === 1) return [rep];
  const newCentral = Int32Array.from(rep.central, (c) => c * k);
  if (!lt || rep.semisimple.every((c) => c === 0)) {
    return [irrepLeviFromParts(mdt, newCentral, rep.semisimple)];
  }
  const chi = symmetricPower(lt, rep.semisimple, k);
  const result = [];
  for (const [wk, mult] of chi) {
    const hw = Int32Array.from(wk.split(",").map(Number));
    for (let i = 0; i < mult; i++) {
      result.push(irrepLeviFromParts(mdt, newCentral, hw));
    }
  }
  return result;
}

// ../HomogeneousTools/Lie.js/src/CompletelyReducibleBundle.js
var CompletelyReducibleBundle = class {
  /**
   * @param {PartialFlagVariety} variety
   * @param {IrrepLevi[]} components
   */
  constructor(variety, components) {
    this.variety = variety;
    this.components = components;
  }
};
function bundleFromWeight(X, lambda) {
  const rep = irrepLevi(X.mdt, lambda);
  return new CompletelyReducibleBundle(X, [rep]);
}
function structureSheaf(X) {
  const R = X.rank;
  return bundleFromWeight(X, new Int32Array(R));
}
function zeroBundle(X) {
  return new CompletelyReducibleBundle(X, []);
}
function lineBundle(X, degrees) {
  const marked = X.markedNodes;
  const R = X.rank;
  if (typeof degrees === "number") degrees = [degrees];
  const w = new Int32Array(R);
  for (let j = 0; j < marked.length; j++) {
    w[marked[j] - 1] = degrees[j];
  }
  return bundleFromWeight(X, w);
}
function tangentBundle(X) {
  const tw = tangentWeights(X.mdt);
  const reps = tw.map((w) => irrepLevi(X.mdt, w));
  return new CompletelyReducibleBundle(X, reps);
}
function cotangentBundle(X) {
  const tw = tangentWeights(X.mdt);
  const reps = tw.map((w) => dualIrrepLevi(irrepLevi(X.mdt, w)));
  return new CompletelyReducibleBundle(X, reps);
}
function canonicalBundle(X) {
  const degs = anticanonicalDegrees(X.mdt);
  return lineBundle(X, degs.map((d) => -d));
}
function anticanonicalBundle(X) {
  const degs = anticanonicalDegrees(X.mdt);
  return lineBundle(X, degs);
}
function rankBundle(E) {
  let r = 0n;
  for (const c of E.components) r += fiberDimension(c);
  return r;
}
function directSum(E, F) {
  return new CompletelyReducibleBundle(E.variety, [...E.components, ...F.components]);
}
function dualBundle(E) {
  return new CompletelyReducibleBundle(E.variety, E.components.map((c) => dualIrrepLevi(c)));
}
function tensorProductBundle(E, F) {
  const eCounts = _countComponents(E.components);
  const fCounts = _countComponents(F.components);
  const result = [];
  for (const [aKey, { rep: a, count: ma }] of eCounts) {
    for (const [bKey, { rep: b, count: mb }] of fCounts) {
      const tp = tensorProductIrrepLevi(a, b);
      const total = ma * mb;
      for (let i = 0; i < total; i++) result.push(...tp);
    }
  }
  return new CompletelyReducibleBundle(E.variety, result);
}
function exteriorPowerBundle(E, k) {
  if (k < 0) return zeroBundle(E.variety);
  if (k === 0) return structureSheaf(E.variety);
  if (k === 1) return E;
  const comps = _deduplicateComponents(E.components);
  const n = comps.length;
  const ranks = comps.map(({ rep }) => Number(fiberDimension(rep)));
  const result = [];
  for (const alpha of _multiexponents(n, k)) {
    if (alpha.some((ai, i) => ai > ranks[i])) continue;
    let factors = null;
    for (let i = 0; i < n; i++) {
      const { rep, count } = comps[i];
      let pieces;
      if (count === 1) {
        pieces = exteriorPowerIrrepLevi(rep, alpha[i]);
      } else {
        pieces = _exteriorPowerMultiple(rep, count, alpha[i]);
      }
      if (pieces.length === 0) {
        factors = null;
        break;
      }
      if (factors === null) {
        factors = pieces;
      } else {
        const newFactors = [];
        for (const a of factors) {
          for (const b of pieces) {
            newFactors.push(...tensorProductIrrepLevi(a, b));
          }
        }
        factors = newFactors;
      }
    }
    if (factors) result.push(...factors);
  }
  return new CompletelyReducibleBundle(E.variety, result);
}
function symmetricPowerBundle(E, k) {
  if (k < 0) return zeroBundle(E.variety);
  if (k === 0) return structureSheaf(E.variety);
  if (k === 1) return E;
  const comps = _deduplicateComponents(E.components);
  const n = comps.length;
  const result = [];
  for (const alpha of _multiexponents(n, k)) {
    let factors = null;
    for (let i = 0; i < n; i++) {
      const { rep, count } = comps[i];
      let pieces;
      if (count === 1) {
        pieces = symmetricPowerIrrepLevi(rep, alpha[i]);
      } else {
        pieces = _symmetricPowerMultiple(rep, count, alpha[i]);
      }
      if (pieces.length === 0 && alpha[i] > 0) {
        factors = null;
        break;
      }
      if (pieces.length === 0) continue;
      if (factors === null) {
        factors = pieces;
      } else {
        const newFactors = [];
        for (const a of factors) {
          for (const b of pieces) {
            newFactors.push(...tensorProductIrrepLevi(a, b));
          }
        }
        factors = newFactors;
      }
    }
    if (factors) result.push(...factors);
  }
  return new CompletelyReducibleBundle(E.variety, result);
}
function twist(E, degree, nodeIdx = 0) {
  const X = E.variety;
  const marked = X.markedNodes;
  const R = X.rank;
  const w = new Int32Array(R);
  w[marked[nodeIdx] - 1] = degree;
  const twistRep = irrepLevi(X.mdt, w);
  const result = [];
  for (const c of E.components) {
    result.push(...tensorProductIrrepLevi(c, twistRep));
  }
  return new CompletelyReducibleBundle(X, result);
}
function detBundle(E) {
  const X = E.variety;
  const R = X.rank;
  const markedSet = new Set(X.markedNodes.map((m) => m - 1));
  const w = new Int32Array(R);
  for (const c of E.components) {
    const d = Number(fiberDimension(c));
    for (let i = 0; i < R; i++) {
      if (markedSet.has(i)) {
        w[i] += d * c.lambda[i];
      }
    }
  }
  return bundleFromWeight(X, w);
}
function universalSubbundle(X) {
  if (!isGeneralizedGrassmannian(X.mdt))
    throw new Error("universalSubbundle requires a generalized Grassmannian");
  const R = X.rank;
  return bundleFromWeight(X, fundamentalWeight(R, 0));
}
function universalQuotientBundle(X) {
  if (!isGeneralizedGrassmannian(X.mdt))
    throw new Error("universalQuotientBundle requires a generalized Grassmannian");
  const R = X.rank;
  const series = X.mdt.dt.series;
  if (series === "A") {
    return bundleFromWeight(X, fundamentalWeight(R, R - 1));
  }
  return dualBundle(universalSubbundle(X));
}
function tautologicalBundles(X) {
  const R = X.rank;
  return X.markedNodes.map(
    (m) => bundleFromWeight(X, fundamentalWeight(R, m - 1))
  );
}
function quotientBundles(X) {
  return tautologicalBundles(X).map((U) => dualBundle(U));
}
function spinorBundle(X, half) {
  const mdt = X.mdt;
  const series = mdt.dt.series;
  const marked = mdt.marked;
  if (marked.length !== 1 || marked[0] !== 1 || series !== "B" && series !== "D")
    throw new Error("spinorBundle requires a quadric (B_m/P_1 or D_m/P_1)");
  const R = X.rank;
  if (series === "B") {
    return bundleFromWeight(X, fundamentalWeight(R, R - 1));
  }
  if (half === void 0) {
    const reps = [
      irrepLevi(X.mdt, fundamentalWeight(R, R - 2)),
      irrepLevi(X.mdt, fundamentalWeight(R, R - 1))
    ];
    return new CompletelyReducibleBundle(X, reps);
  }
  if (half === "plus") {
    return bundleFromWeight(X, fundamentalWeight(R, R - 2));
  }
  if (half === "minus") {
    return bundleFromWeight(X, fundamentalWeight(R, R - 1));
  }
  throw new Error("half must be 'plus' or 'minus'");
}
function scalarMultiplyBundle(n, E) {
  if (n < 0) throw new RangeError("Cannot multiply a bundle by a negative integer");
  if (n === 0) return zeroBundle(E.variety);
  const repeated = [];
  for (let i = 0; i < n; i++) repeated.push(...E.components);
  return new CompletelyReducibleBundle(E.variety, repeated);
}
function _countComponents(components) {
  const map = /* @__PURE__ */ new Map();
  for (const c of components) {
    const k = c.key();
    if (map.has(k)) {
      map.get(k).count++;
    } else {
      map.set(k, { rep: c, count: 1 });
    }
  }
  return map;
}
function _deduplicateComponents(components) {
  const map = _countComponents(components);
  return [...map.values()];
}
function _multiexponents(n, k) {
  if (n === 0) return k === 0 ? [[]] : [];
  if (n === 1) return [[k]];
  const result = [];
  for (let i = 0; i <= k; i++) {
    for (const rest of _multiexponents(n - 1, k - i)) {
      result.push([i, ...rest]);
    }
  }
  return result;
}
function _exteriorPowerMultiple(rep, count, k) {
  if (k === 0) {
    return [irrepLeviFromParts(rep.mdt, new Int32Array(rep.central.length), new Int32Array(rep.semisimple.length))];
  }
  if (count === 1) return exteriorPowerIrrepLevi(rep, k);
  const result = [];
  const dimV = Number(fiberDimension(rep));
  for (let i = 0; i <= Math.min(k, dimV); i++) {
    const j = k - i;
    const left = exteriorPowerIrrepLevi(rep, i);
    const right = _exteriorPowerMultiple(rep, count - 1, j);
    if (left.length === 0 || right.length === 0) continue;
    for (const a of left) {
      for (const b of right) {
        result.push(...tensorProductIrrepLevi(a, b));
      }
    }
  }
  return result;
}
function _symmetricPowerMultiple(rep, count, k) {
  if (k === 0) {
    return [irrepLeviFromParts(rep.mdt, new Int32Array(rep.central.length), new Int32Array(rep.semisimple.length))];
  }
  if (count === 1) return symmetricPowerIrrepLevi(rep, k);
  const result = [];
  for (let i = 0; i <= k; i++) {
    const j = k - i;
    const left = symmetricPowerIrrepLevi(rep, i);
    const right = _symmetricPowerMultiple(rep, count - 1, j);
    if (left.length === 0 || right.length === 0) continue;
    for (const a of left) {
      for (const b of right) {
        result.push(...tensorProductIrrepLevi(a, b));
      }
    }
  }
  return result;
}

// ../HomogeneousTools/Lie.js/src/Cohomology.js
function cohomology(E) {
  const d = dimension(E.variety);
  const dt = E.variety.mdt.dt;
  const entries = [];
  for (let i = 0; i <= d; i++) entries.push(/* @__PURE__ */ new Map());
  for (const comp of E.components) {
    const lambda = comp.lambda;
    const bwb = borelWeylBott(dt, lambda);
    if (bwb === null) continue;
    const { degree, weight } = bwb;
    if (degree >= 0 && degree <= d) {
      const k = weightKey(weight);
      const m = entries[degree];
      m.set(k, (m.get(k) || 0) + 1);
    }
  }
  return entries;
}
function cohomologyDimensions(E) {
  const chars = cohomology(E);
  const dt = E.variety.mdt.dt;
  return chars.map((chi) => {
    let total = 0n;
    for (const [k, mult] of chi) {
      const w = Int32Array.from(k.split(",").map(Number));
      total += BigInt(mult) * weylDimension(dt, w);
    }
    return total;
  });
}
function eulerCharacteristicBundle(E) {
  const dt = E.variety.mdt.dt;
  const counts = /* @__PURE__ */ new Map();
  for (const comp of E.components) {
    const k = weightKey(comp.lambda);
    counts.set(k, counts.get(k) || { lambda: comp.lambda, count: 0 });
    counts.get(k).count++;
  }
  let result = 0n;
  for (const [, { lambda, count }] of counts) {
    const bwb = borelWeylBott(dt, lambda);
    if (bwb === null) continue;
    const { degree, weight } = bwb;
    const sign = degree % 2 === 0 ? 1n : -1n;
    result += sign * BigInt(count) * weylDimension(dt, weight);
  }
  return result;
}

// ../HomogeneousTools/Lie.js/src/Constructions.js
function Gr(k, n) {
  if (k < 1 || k >= n) throw new RangeError(`Gr(k, n) requires 1 \u2264 k < n, got k=${k}, n=${n}`);
  return partialFlagVariety(new SimpleDynkinType("A", n - 1), [k]);
}
function projectiveSpace(n) {
  return partialFlagVariety(new SimpleDynkinType("A", n), [1]);
}
function flagVariety(n, marked) {
  return partialFlagVariety(new SimpleDynkinType("A", n - 1), marked);
}
function OGr(k, n) {
  if (n % 2 === 1) {
    const r = (n - 1) / 2;
    return partialFlagVariety(new SimpleDynkinType("B", r), [k]);
  } else {
    const r = n / 2;
    return partialFlagVariety(new SimpleDynkinType("D", r), [k]);
  }
}
function SGr(k, n) {
  if (n % 2 !== 0) throw new RangeError(`SGr(k, n) requires even n, got n=${n}`);
  return partialFlagVariety(new SimpleDynkinType("C", n / 2), [k]);
}
function LGr(n) {
  return SGr(n, 2 * n);
}
var IGr = OGr;
function quadric(n) {
  if (n === 1) {
    return partialFlagVariety(new SimpleDynkinType("A", 1), [1]);
  }
  if (n === 2) {
    return partialFlagVariety(
      new ProductDynkinType([new SimpleDynkinType("A", 1), new SimpleDynkinType("A", 1)]),
      [1, 2]
    );
  }
  if (n % 2 === 1) {
    return partialFlagVariety(new SimpleDynkinType("B", (n + 1) / 2), [1]);
  }
  if (n === 4) {
    return partialFlagVariety(new SimpleDynkinType("A", 3), [1]);
  }
  return partialFlagVariety(new SimpleDynkinType("D", n / 2 + 1), [1]);
}
function cayleyPlane() {
  return partialFlagVariety(new SimpleDynkinType("E", 6), [1]);
}
function freudenthalVariety() {
  return partialFlagVariety(new SimpleDynkinType("E", 7), [7]);
}
function adjointVariety(dt) {
  const { series, n } = dt;
  let node;
  switch (series) {
    case "A":
      return partialFlagVariety(dt, [1, n]);
    case "B":
      node = 2;
      break;
    case "C":
      node = 1;
      break;
    case "D":
      node = 2;
      break;
    case "E":
      if (n === 6) node = 2;
      else if (n === 7) node = 1;
      else node = 8;
      break;
    case "F":
      node = 1;
      break;
    case "G":
      node = 2;
      break;
    default:
      throw new Error(`Unknown series ${series}`);
  }
  return partialFlagVariety(dt, [node]);
}
function coadjointVariety(dt) {
  const { series, n } = dt;
  let node;
  switch (series) {
    case "A":
      return partialFlagVariety(dt, [1, n]);
    case "B":
      node = 1;
      break;
    case "C":
      node = 2;
      break;
    case "D":
      node = 2;
      break;
    case "E":
      if (n === 6) node = 2;
      else if (n === 7) node = 1;
      else node = 8;
      break;
    case "F":
      node = 4;
      break;
    case "G":
      node = 1;
      break;
    default:
      throw new Error(`Unknown series ${series}`);
  }
  return partialFlagVariety(dt, [node]);
}
export {
  CompletelyReducibleBundle,
  Gr,
  IGr,
  IrrepLevi,
  LGr,
  MarkedDynkinType,
  OGr,
  ONE,
  PartialFlagVariety,
  ProductDynkinType,
  Rational,
  SGr,
  SimpleDynkinType,
  ZERO,
  adamsOperator,
  adjointVariety,
  anticanonicalBundle,
  anticanonicalDegrees,
  applyW0,
  bettiNumbers,
  bigGcd,
  borelWeylBott,
  bundleFromWeight,
  canonicalBundle,
  cartanBilinearForm,
  cartanDeterminant,
  cartanMatrix,
  cartanMatrixInverse,
  cartanSymmetrizer,
  cayleyPlane,
  centralPart,
  centralScalingFactor,
  charAdd,
  charDimension,
  charScale,
  charSub,
  coadjointVariety,
  cohomology,
  cohomologyDimensions,
  conjugateDominantWeight,
  conjugateDominantWeightWithLength,
  cotangentBundle,
  coxeterCoefficients,
  coxeterNumber,
  decompositionMatrix,
  decompositionMatrixInv,
  degreesFundamentalInvariants,
  detBundle,
  dimension,
  directSum,
  dominantCharacter,
  dominantWeights,
  dualBundle,
  dualChar,
  dualCoxeterCoefficients,
  dualCoxeterNumber,
  dualIrrepLevi,
  dualWeight,
  eulerCharacteristic,
  eulerCharacteristicBundle,
  exteriorPower,
  exteriorPowerBundle,
  exteriorPowerIrrepLevi,
  fanoIndex,
  fiberDimension,
  flagVariety,
  freudenthalVariety,
  fullFlagVariety,
  fundamentalWeight,
  gcd,
  highestRoot,
  highestShortRoot,
  highestWeight,
  hilbertSeries,
  irrepChar,
  irrepLevi,
  irrepLeviFromParts,
  isAdjoint,
  isCoadjoint,
  isCominuscule,
  isDominant,
  isEffective,
  isFullFlag,
  isGeneralizedGrassmannian,
  isIrreducible,
  isMinuscule,
  isSingular,
  lcm,
  leviPermutation,
  leviType,
  lineBundle,
  longestElement,
  lrTensorProduct,
  omegaBilinearFormScaled,
  parseDynkinType,
  partialFlagVariety,
  plethysm,
  positiveNonparabolicRoots,
  positiveParabolicRoots,
  projectiveSpace,
  quadric,
  quotientBundles,
  rankBundle,
  rat,
  reflect,
  reflectInPlace,
  rootHeight,
  rootSystem,
  rootToWeight,
  rootWeightPairing,
  scalarMultiplyBundle,
  spinorBundle,
  structureSheaf,
  symmetricPower,
  symmetricPowerBundle,
  symmetricPowerIrrepLevi,
  tangentBundle,
  tangentWeights,
  tautologicalBundles,
  tensorProduct,
  tensorProductBundle,
  tensorProductChars,
  tensorProductIrrepLevi,
  twist,
  universalQuotientBundle,
  universalSubbundle,
  weightAdd,
  weightKey,
  weightNeg,
  weightScale,
  weightSub,
  weylDimension,
  weylOrbit,
  weylOrder,
  weylVector,
  weylloop,
  zeroBundle
};
