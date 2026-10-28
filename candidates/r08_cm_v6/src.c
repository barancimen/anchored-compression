// v6 (readable source; p.c = this file minified by lab/minify2.py).
// One program: stdin of exactly 125000 bytes = image -> compress to stdout; otherwise the
// stdin is an archive -> decompress to stdout. Each 4-bit pixel is coded MSB-first as 4
// binary decisions. Each decision mixes (logistic mixer, one shared weight set):
//   4 adaptive counters in hashed contexts of the causal neighbours, and
//   a discretized-Gaussian input centred on an online least-squares (RLS, forgetting .9997)
//   prediction from the 84 causal neighbours within radius 6 (one Gauss-Seidel sweep per
//   pixel), with a variance learned per bucket of nearby prediction errors.
#include <stdio.h>
#include <math.h>
#define F(i, n) for (i = 0; i < n; i++)
#define S 528  // row stride of the padded image I: 8-pixel left margin, 20-pixel right margin
#define L 84   // RLS inputs: causal neighbours within radius 6
typedef double D;
typedef unsigned U;
unsigned char B[1 << 17], O[125000], I[S * S], N[1 << 22];
float P[1 << 22];  // counter probabilities; N = update counts (rate 1/(n+1.5), n <= 60)
U h[4], H[4], x, z = -1, X, d, u;
D E[S * S], w[5], s[5], o[L], a[L], C[L][L], b[L], V[32], m, g;

D st(D p) { return p = p * .9998 + .0001, log(p / (1 - p)); }
// P(pixel < v) under N(m, g^2 / 2), with the clipped tails folded into 0 and 15
D G(int v) { return v < 1 ? 0 : v > 15 ? 1 : erfc((m + .5 - v) / g) / 2; }

int main() {
  int i, j, k, y, c, e, r, l, v, q, n = 0;
  d = fread(B, 1, 1 << 17, stdin) != 125000;
  F(i, 1 << 22) P[i] = .5;
  F(i, S * S) I[i] = 8;
  F(i, 32) V[i] = 1, w[i % 5] = .1;
  if (d) F(i, 4) X = X << 8 | B[u++];
  F(r, 500) F(l, 500) {
    k = (r + 8) * S + l + 8;
    if (!l) F(j, 8) I[k - j - 1] = I[k - S];  // replicate the left edge
    U W = I[k - 1], A = I[k - S], Q = I[k - S - 1], R = I[k - S + 1], T = I[k - 2], Y = I[k - 2 * S];
    e = m = 0;
    F(i, 7) for (j = -6; j < 7; j++) if (i || j < 0) m += a[e] * (o[e] = I[k - i * S + j]), e++;
    q = fmin((E[k - 1] + E[k - S] + (E[k - S - 1] + E[k - S + 1]) / 2) * 4, 31);
    g = sqrt(V[q] * 2);
    h[0] = W << 8 | A << 4 | R;
    h[1] = A << 8 | Y << 4 | R;
    h[2] = W << 8 | T << 4 | Q;
    h[3] = W << 20 | A << 16 | Q << 12 | R << 8 | T << 4 | Y;
    v = B[n / 2] >> (n & 1 ? 0 : 4) & 15;
    for (c = 1, e = 3; e >= 0; e--) {
      F(i, 4) s[i] = st(P[H[i] = (h[i] << 4 | c | i << 28) * 2654435761u >> 10]);
      j = (c << e + 1) - 16;  // remaining values j .. j+2^(e+1)-1, split at j+2^e
      s[4] = st((G(j + (2 << e)) - G(j + (1 << e))) / (G(j + (2 << e)) - G(j) + 1e-30));
      D f = 0, p;
      F(i, 5) f += w[i] * s[i];
      p = 1 / (1 + exp(-f));
      // binary arithmetic coder (carry-less, 32-bit range)
      U M = x + (z - x) * (p * .9999 + .00005);
      y = d ? X <= M : v >> e & 1;
      y ? (z = M) : (x = M + 1);
      while ((x ^ z) < 1 << 24) {
        if (d) X = X << 8 | B[u++]; else putchar(z >> 24);
        x <<= 8, z = z << 8 | 255;
      }
      F(i, 5) w[i] += .001 * (y - p) * s[i];
      F(i, 4) N[H[i]] += N[H[i]] < 60, P[H[i]] += (y - P[H[i]]) / (N[H[i]] + .5);
      c = c * 2 + y;
    }
    I[k] = c -= 16;
    O[n / 2] |= c << (n & 1 ? 0 : 4);
    if (l > 498) F(j, 20) I[k + j] = c;  // replicate the right edge
    E[k] = fabs(c - m);
    V[q] += (E[k] * E[k] - V[q]) / 500;
    // RLS with forgetting .9997, then one Gauss-Seidel sweep for the weights
    F(i, L) {
      F(j, L) C[i][j] = .9997 * C[i][j] + .0003 * o[i] * o[j];
      b[i] = .9997 * b[i] + .0003 * o[i] * c;
      D f = b[i];
      F(j, L) if (j != i) f -= C[i][j] * a[j];
      a[i] = f / (C[i][i] + .03);
    }
    n++;
  }
  d ? fwrite(O, 1, n / 2, stdout) : putchar((x >> 24) + 1);
}
