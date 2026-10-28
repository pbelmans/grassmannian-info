/* Classical root-system pages outside the prebuilt catalogue. */
(function () {
  "use strict";

  function model(T, n) {
    const q = BigInt(n);
    const roots = T === "A" ? q * (q + 1n) : T === "D" ? 2n * q * (q - 1n) : 2n * q * q;
    let order = 1n;
    for (let i = 2n; i <= (T === "A" ? q + 1n : q); i++) order *= i;
    if (T !== "A") order *= 2n ** (T === "D" ? q - 1n : q);
    const exponents = Array.from({length: T === "D" ? n - 1 : n}, (_, i) => T === "A" ? i + 1 : 2 * i + 1);
    if (T === "D") exponents.push(n - 1);
    exponents.sort((a, b) => a - b);
    const cartan = Array.from({length: n}, (_, i) => Array.from({length: n}, (_, j) => i === j ? 2 : 0));
    for (let i = 0; i < (T === "D" ? n - 2 : n - 1); i++) cartan[i][i + 1] = cartan[i + 1][i] = -1;
    if (T === "B") cartan[n - 2][n - 1] = -2;
    if (T === "C") cartan[n - 1][n - 2] = -2;
    if (T === "D") cartan[n - 3][n - 1] = cartan[n - 1][n - 3] = -1;
    return {
      dimension: String(roots + q), roots: String(roots), order: String(order),
      coxeter: String(T === "A" ? q + 1n : T === "D" ? 2n * q - 2n : 2n * q),
      exponents, cartan,
    };
  }

  function diagram(T, n, siteRoot) {
    const line = points => `<polyline points="${points}" />`;
    const node = (i, x, y) => `<a href="${siteRoot}${T}${n}/${i}"><circle cx="${x}" cy="${y}" r="5" data-label="${T}${n}-${i}" /></a>`;
    let drawing = "";
    if (T === "D") {
      drawing += line(`0,10 ${(n - 3) * 20},10`);
      drawing += line(`${(n - 3) * 20},10 ${(n - 2) * 20},0`);
      drawing += line(`${(n - 3) * 20},10 ${(n - 2) * 20},20`);
      for (let i = 1; i <= n - 2; i++) drawing += node(i, (i - 1) * 20, 10);
      drawing += node(n - 1, (n - 2) * 20, 0) + node(n, (n - 2) * 20, 20);
    } else {
      drawing += line(`0,10 ${(n - (T === "A" ? 1 : 2)) * 20},10`);
      if (T !== "A") {
        const x = (n - 2) * 20;
        drawing += line(`${x},8 ${x + 20},8`) + line(`${x},12 ${x + 20},12`);
        drawing += T === "B" ? line(`${x + 7},17 ${x + 13},10 ${x + 7},3`) : line(`${x + 13},17 ${x + 7},10 ${x + 13},3`);
      }
      for (let i = 1; i <= n; i++) drawing += node(i, (i - 1) * 20, 10);
    }
    return `<svg class="dynkin" height="100" viewBox="-10 -10 ${20 * (T === "D" ? n - 1 : n)} 40"><g style="stroke-width:1;stroke:black;fill:white">${drawing}</g></svg>`;
  }

  function render(T, n, templates, siteRoot = "/") {
    const D = model(T, n);
    const values = {
      RANK: n, RANK_PLUS_1: n + 1, RANK_MINUS_1: n - 1,
      DIMENSION: D.dimension, ROOTS: D.roots, WEYL_ORDER: D.order,
      COXETER: D.coxeter, EXPONENTS: D.exponents.join(","),
      CARTAN: D.cartan.map(row => row.join(" & ")).join(" \\\\ "),
    };
    const fill = text => text.replace(/__([A-Z_0-9]+)__/g, (_, key) => String(values[key]));
    const card = (title, body) => `<div class="card mb-4"><div class="card-header">${title}</div><div class="card-body">${fill(body)}</div></div>`;
    return `<header class="text-center"><h2 class="mt-3">Dynkin type $\\mathrm{${T}}_{${n}}$</h2>
      <div class="row"><div class="col-2"><a class="h4" href="${siteRoot}${T}${n - 1}">&laquo;$\\mathrm{${T}}_{${n - 1}}$</a></div>
      <div class="col-8">${diagram(T, n, siteRoot)}</div><div class="col-2"><a class="h4" href="${siteRoot}${T}${n + 1}">$\\mathrm{${T}}_{${n + 1}}$&raquo;</a></div></div></header>
      <div class="row" id="type-information"><div class="col-12">
      ${card("Basic information", templates.basic)}${card("Description of the root system", templates.roots)}${card("Cartan matrix", templates.cartan)}</div></div>`;
  }

  // Export the arithmetic for comparison with Python's exact models.
  if (typeof module !== "undefined") module.exports = {model, render};
  if (typeof document === "undefined") return;
  const siteRoot = document.documentElement.dataset.siteRoot || "/";
  if (!location.pathname.startsWith(siteRoot)) return;
  const route = "/" + location.pathname.slice(siteRoot.length);
  const match = route.match(/^\/([ABCD])([1-9][0-9]*)\/?$/);
  if (!match) return;
  const T = match[1], n = Number(match[2]);
  if (!Number.isSafeInteger(n) || n < 21) return;
  const templates = JSON.parse(document.getElementById("type-templates").textContent);
  const target = document.getElementById("dynamic-type");
  target.innerHTML = render(T, n, templates[T], siteRoot);
  document.title = `Dynkin type ${T}${n}—grassmannian.info`;
  MathJax.Hub.Config({displayAlign: "left", "HTML-CSS": {linebreaks: {automatic: true}}});
  MathJax.Hub.Queue(["Typeset", MathJax.Hub, target]);
}());
