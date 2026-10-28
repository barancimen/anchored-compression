// Lab template of the quadtree candidate for context search (lab/qtsearch.py): 002 with the
// image transposed (pixel (x, y) of the scan is image pixel (y, x)), the 12 nearest stencil
// pixels available to contexts, and the contexts as macros CX0..CX4. Pixels are coded in three passes by class g = x % 2 + y % 2:
//   g = 0  anchors: the lattice x, y even, raster order over the lattice, predicted from
//          earlier anchors (stencil in lattice units);
//   g = 2  cell centres (x, y odd): their 4 diagonal neighbours are anchors;
//   g = 1  edge midpoints: their 4 axial neighbours are known.
// Each pass has its own stencil: the known pixels within radius 6, nearest first.
// One program: stdin of exactly 125000 bytes = image -> compress to stdout; otherwise the
// stdin is an archive -> decompress to stdout. Each 4-bit pixel is coded MSB-first as 4
// binary decisions. Each decision mixes (logistic mixer, one shared weight set):
//   5 adaptive counters in hashed contexts of the nearest stencil pixels v0..v5:
//   (v0,v1,v2,v3) (v0,v1) (v2,v3) (prediction, error bucket) (v0,v1,v4,v5), and
//   a discretized-Gaussian input centred on an online least-squares (RLS, forgetting .9997)
//   prediction from the stencil (as a residual from the mean of the 4 nearest; one
//   Gauss-Seidel sweep per pixel), with a variance learned per bucket of nearby errors.
// RLS weights and variances are per pass; counter contexts are hashed with the pass.
#include <stdio.h>
#include <math.h>
#define F(i, n) for (i = 0; i < n; i++)
#define Z 160  // max stencil length
#ifndef CX0  // defaults: the contexts found by lab/qtsearch.py (candidate 004)
#define CX0 r[9] << 8 | r[1] << 4 | r[0]
#define CX1 r[0] << 8 | r[1] << 4 | r[10]
#define CX2 r[2] << 4 | r[3]
#define CX3 (int)(m * 2 + .5) << 5 | q
#define CX4 r[7] << 12 | r[11] << 8 | r[4] << 4 | r[5]
#endif
#ifndef RB
#define RB 37  // stencil: offsets with dx^2 + dy^2 < RB inside the box |dx|, |dy| <= BOX
#endif
#ifndef BOX
#define BOX 6
#endif
#ifndef LR
#define LR .001
#endif
#ifndef LIM
#define LIM 80
#endif
#ifndef LAM
#define LAM .9997
#endif
#ifndef NU
#define NU 100
#endif
#ifndef KQ
#define KQ 3
#endif
#ifndef VR
#define VR 500
#endif
typedef double D;
unsigned char B[1 << 17], O[125000], N[1 << 22];
unsigned lo, hi = -1, X, M;
// n, sx, sy: stencil length and offsets per pass; r: nearest stencil pixels; h: contexts;
// H: counter slots (H[5] stays 0: a dummy slot so one loop updates everything)
int I[250000], n[3], sx[3][Z], sy[3][Z], r[12], h[5], H[6], d, u, i, j, k, x, y, c, e, v, q, g, p, L, bt, z;
// P: counter probabilities (N = update counts, rate 1/(n+1.5), n <= 80); E: |RLS error|;
// C: exponentially weighted RLS sums per pass, column L holds the right-hand side (o[L] =
// pixel and a[g][L] = -1 fold it into the Gauss-Seidel sum); V: 2 x error variance per bucket
D P[1 << 22], E[250000], w[6], s[6], o[Z + 1], a[3][Z + 1], C[3][Z][Z + 1], V[3][32], m, sg, t, f, b;

D st(D p) { return p = p * .9998 + .0001, log(p / (1 - p)); }
// P(pixel < (v << e) - 16) under N(m, sg^2 / 2); clipped tails folded into 0 and 15
D G(int v) { return v <<= e, v < 17 ? 0 : v > 31 ? 1 : erfc((m + 16.5 - v) / sg) / 2; }
int R(int v) { return v < 0 ? -v : v > 499 ? 998 - v : v; }  // reflect at the borders

int main() {
  d = fread(B, 1, 1 << 17, stdin) != 125000;
  F(i, 1 << 22) P[i] = .5, I[i % 250000] = 8, V[i % 3][i % 32] = 1;
  // stencils: g = 0 causal; g = 1 dx + dy odd, or causal; g = 2 dx, dy both odd, or both
  // even and causal. Offsets with dx^2 + dy^2 = e for e = 1, 2, ... (nearest first)
  F(g, 3) {
    for (e = 1; e < RB; e++)
      for (y = -BOX; y <= BOX; y++)
        for (x = -BOX; x <= BOX; x++)
          if (x * x + y * y == e && (c = y < 0 || !y && x < 0, g ? g & 1 ? (x + y) & 1 || c : x & y & 1 || !((x | y) & 1) && c : c))
            sx[g][n[g]] = x, sy[g][n[g]++] = y;
    a[g][n[g]] = -1;
  }
  F(i, 4) X = X << 8 | B[u++];
  F(p, 3) F(y, 500) F(x, 500) if ((g = x % 2 + y % 2) == p * 2 % 3) {
    k = y * 500 + x;
    L = n[g];
    F(i, L) o[i] = I[j = R(y + sy[g][i] * (2 - !!g)) * 500 + R(x + sx[g][i] * (2 - !!g))], i < 12 ? r[i] = j : 0;
    b = (o[0] + o[1] + o[2] + o[3]) / 4;
    m = b;
    F(i, L) m += a[g][i] * (o[i] -= b);
    q = fmin((E[r[0]] + E[r[1]] + E[r[2]] + E[r[3]]) * KQ, 31);
    sg = sqrt(V[g][q] * 2);
    F(i, 12) r[i] = I[r[i]];
    h[0] = CX0;
    h[1] = CX1;
    h[2] = CX2;
    h[3] = CX3;
    h[4] = CX4;
    z = x * 500 + y;
    v = B[z / 2] >> 4 - z % 2 * 4 & 15;
    // c = 2 x (bits so far, with a leading 1); the remaining values are (c << e) - 16 ...
    for (c = 1, e = 4; e--; c += bt) {
      c *= 2;
      t = G(c + 2), f = 0;
      s[5] = st((t - G(c + 1)) / (t - G(c) + 1e-30));
      F(i, 5) s[i] = st(P[H[i] = (h[i] << 5 | c | g << 27 | i << 29) * 2654435761u >> 10]);
      F(i, 6) f += w[i] * s[i];
      t = 1 / (1 + exp(-f));
      // binary arithmetic coder (carry-less, 32-bit range)
      M = lo + (hi - lo) * (t * .9999 + .00005);
      bt = d ? X <= M : v >> e & 1;
      bt ? hi = M : (lo = M + 1);
      for (; (lo ^ hi) < 1 << 24; lo <<= 8, hi = hi << 8 | 255) d ? X = X << 8 | B[u++] : putchar(hi >> 24);
      F(i, 6) w[i] += LR * (bt - t) * s[i], j = H[i], N[j] += N[j] < LIM, P[j] += (bt - P[j]) / (N[j] + .5);
    }
    I[k] = c -= 16;
    O[z / 2] |= c << 4 - z % 2 * 4;
    E[k] = fabs(c - m);
    V[g][q] += (E[k] * E[k] - V[g][q]) / VR;
    // RLS sums with forgetting .9997 (regulariser 100), then one Gauss-Seidel sweep
    o[L] = c - b;
    F(i, L) {
      m = 0;
      F(j, L + 1) C[g][i][j] = C[g][i][j] * LAM + o[i] * o[j], m -= j - i ? C[g][i][j] * a[g][j] : 0;
      a[g][i] = m / (C[g][i][i] + NU);
    }
  }
  d ? fwrite(O, 1, 125000, stdout) : putchar((lo >> 24) + 1);
}
