// 012: 005_quadtree_golf plus a non-local means input (readable source; p.c = this file minified
// by lab/minify2.py). Family B with the best model's extra input: compare with 007. As in 005: Lattice anchors + quadtree
// refinement; the image is coded transposed (scan pixel (x, y) is image pixel (y, x)).
// Pixels are coded in three passes by class g = x % 2 + y % 2:
//   g = 0  anchors: the lattice x, y even, raster order over the lattice, predicted from
//          earlier anchors (stencil offsets doubled: lattice units);
//   g = 2  cell centres (x, y odd): their 4 diagonal neighbours are anchors;
//   g = 1  edge midpoints: their 4 axial neighbours are known.
// Each pass has its own stencil: the known pixels within radius 6, nearest first.
// One program: stdin of exactly 125000 bytes = image -> compress to stdout; otherwise the
// stdin is an archive -> decompress to stdout. Each 4-bit pixel is coded MSB-first as 4
// binary decisions. Each decision mixes (logistic mixer, one shared weight set):
//   5 adaptive counters in hashed contexts of the nearest stencil pixels v0..v11:
//   (v9,v1,v0) (v0,v1,v10) (v2,v3) (prediction, error bucket) (v7,v11,v4,v5), and
//   a discretized-Gaussian input centred on an online least-squares (RLS, forgetting .9995)
//   prediction from the stencil (as a residual from the mean of the 4 nearest; one
//   Gauss-Seidel sweep per pixel), with a variance learned per (pass, error bucket), and
//   a second discretized-Gaussian input centred on a non-local means prediction: the pixels
//   of the same parity up to 8 pixels away in each direction that this pass has coded already
//   are weighted by exp(-d / 2), d = the squared difference between their 12 nearest stencil
//   pixels and ours (which this pass's stencil makes known on both sides); own error buckets.
#include <stdio.h>
#include <math.h>
#define F(i, n) for (i = 0; i < n; i++)
#define R(v) (499 - abs(499 - abs(v)))  // reflect at the borders
typedef double D;
unsigned char B[1 << 17], O[125000], N[1 << 22];
// l, U: coder low and high; X: decoder window
unsigned l, U = -1, X, M;
// n, J, T: stencil length and offsets per pass; r: nearest stencil pixels; h: contexts;
// H: counter slots (H[5] stays 0: a dummy slot so one loop updates everything); Y: bit
// Q: error bucket of the non-local means input; dx, dy: candidate offset
int I[250000], n[3], J[3][96], T[3][96], r[12], h[5], H[7], d, u, i, j, k, x, y, c, e, v, q, g, p, L, Y, z, Q, dx, dy;
// P: counter probabilities (N = update counts, rate 1/(n+1.5), n <= 80); E: |RLS error|;
// C: exponentially weighted RLS sums per pass, column L holds the right-hand side (o[L] =
// pixel and a[g][L] = -1 fold it into the Gauss-Seidel sum); V: 2 x error variance per
// bucket q (q includes the pass; see below)
// bucket q (q includes the pass; 96..191 for the non-local means input); E2: its |errors|;
// Mu, Sg: mean and sd of the two Gaussian inputs; Z, W: non-local means weight sums
D P[1 << 22], E[250000], E2[250000], w[7], s[7], o[97], a[3][97], C[3][96][97], V[192], m, K, t, f, b, Mu[2], Sg[2], Z, W;

D A(D p) { return p = p * .9998 + .0001, log(p / (1 - p)); }
// P(pixel < (v << e) - 16) under N(Mu[j], Sg[j]^2 / 2); clipped tails folded into 0 and 15
D G(int v, int j) { return v <<= e, v < 17 ? 0 : v > 31 ? 1 : erfc((Mu[j] + 16.5 - v) / Sg[j]) / 2; }

int main() {
  d = fread(B, 1, 1 << 17, stdin) != 125000;
  F(i, 1 << 22) P[i] = .5, I[i % 250000] = 8, V[i % 192] = 1;
  // stencils: g = 0 causal; g = 1 dx + dy odd, or causal; g = 2 dx, dy both odd, or both
  // even and causal. Offsets with dx^2 + dy^2 = e for e = 1, 2, ... (nearest first)
  F(g, 3) {
    for (e = 1; e < 37; e++)
      for (y = -6; y < 7; y++)
        for (x = -6; x < 7; x++)
          if (x * x + y * y == e && (c = y < 0 || !y && x < 0, g ? g & 1 ? x + y & 1 || c : x & y & 1 || ~(x | y) & 1 && c : c))
            J[g][n[g]] = x << !g, T[g][n[g]++] = y << !g;
    a[g][n[g]] = -1;
  }
  F(i, 4) X = X << 8 | B[u++];
  F(p, 3) F(y, 500) F(x, 500) if ((g = x % 2 + y % 2) == p * 2 % 3) {
    k = y * 500 + x;
    L = n[g], t = Z = 0;
    // stencil values; r = the 12 nearest; t, Z = their 4 nearest's errors (RLS, non-local means)
    F(i, L) o[i] = I[j = R(y + T[g][i]) * 500 + R(x + J[g][i])], i < 12 ? r[i] = o[i] : 0, t += i < 4 ? E[j] : 0, Z += i < 4 ? E2[j] : 0;
    Q = 96 + g * 32 + fmin(Z * 2, 31);
    Sg[1] = sqrt(V[Q] * 2);
    // non-local means: same-parity pixels (dx, dy even) up to 8 away, coded earlier in this pass
    Z = W = 0;
    for (dy = -8; dy < 1; dy += 2)
      for (dx = -8; dx < 9 && (dy || dx < 0); dx += 2)
        if (y + dy >= 0 && x + dx >= 0 && x + dx < 500) {
          t = 0;
          F(i, 12) f = r[i] - I[R(y + dy + T[g][i]) * 500 + R(x + dx + J[g][i])], t += f * f;
          Z += t = exp(-t / 2), W += t * I[(y + dy) * 500 + x + dx];
        }
    Mu[1] = W / (Z + 1e-99);
    t = 0;
    F(i, 4) t += E[R(y + T[g][i]) * 500 + R(x + J[g][i])];
    m = b = (o[0] + o[1] + o[2] + o[3]) / 4;
    F(i, L) m += a[g][i] * (o[i] -= b);
    q = g * 32 + fmin(t * 2, 31);
    Sg[0] = sqrt(V[q] * 2), Mu[0] = m;
    h[0] = r[9] << 8 | r[1] << 4 | r[0];
    h[1] = r[0] << 8 | r[1] << 4 | r[10];
    h[2] = r[2] << 4 | r[3];
    h[3] = (int)(m * 2 + .5) << 7 | q;
    h[4] = r[7] << 12 | r[11] << 8 | r[4] << 4 | r[5];
    z = x * 500 + y;  // image pixel (transposed)
    v = B[z / 2] >> 4 - z % 2 * 4 & 15;
    // c = 2 x (bits so far, with a leading 1); the remaining values are (c << e) - 16 ...
    for (c = 1, e = 4; e--; c += Y) {
      c *= 2;
      f = 0;
      F(j, 2) t = G(c + 2, j), s[5 + j] = A((t - G(c + 1, j)) / (t - G(c, j) + 1e-30));
      F(i, 5) s[i] = A(P[H[i] = (h[i] << 5 | c | g << 27 | i << 29) * 2654435761u >> 10]);
      F(i, 7) f += w[i] * s[i];
      t = 1 / (1 + exp(-f));
      // binary arithmetic coder (carry-less, 32-bit range)
      M = l + (U - l) * (t * .9999 + .00005);
      Y = d ? X <= M : v >> e & 1;
      Y ? U = M : (l = M + 1);
      for (; (l ^ U) < 1 << 24; l <<= 8, U = U << 8 | 255) d ? X = X << 8 | B[u++] : putchar(U >> 24);
      F(i, 7) w[i] += .001 * (Y - t) * s[i], j = H[i], N[j] += N[j] < 80, P[j] += (Y - P[j]) / (N[j] + .5);
    }
    I[k] = c -= 16;
    O[z / 2] |= c << 4 - z % 2 * 4;
    E[k] = fabs(c - m);
    V[q] += (E[k] * E[k] - V[q]) / 300;
    E2[k] = fabs(c - Mu[1]);
    V[Q] += (E2[k] * E2[k] - V[Q]) / 300;
    // RLS sums with forgetting .9995 (regulariser 100), then one Gauss-Seidel sweep
    o[L] = c - b;
    F(i, L) {
      m = 0;
      F(j, L + 1) C[g][i][j] = C[g][i][j] * .9995 + o[i] * o[j], m -= j - i ? C[g][i][j] * a[g][j] : 0;
      a[g][i] = m / (C[g][i][i] + 100);
    }
  }
  d ? fwrite(O, 1, 125000, stdout) : putchar((l >> 24) + 1);
}
