// v4: size-conscious rewrite of v3 (+ radius-4 OLS solved by one Gauss-Seidel sweep per
// pixel, no APM). One program: stdin of exactly 125000 bytes = image -> compress, else
// decompress. Pixels are coded MSB-first as 4 binary decisions; each decision mixes
// 4 context counters + a discretized-Gaussian input from the OLS prediction + bias.
#include <stdio.h>
#include <math.h>
#define F(i, n) for (i = 0; i < n; i++)
#define S 512
#define L 41
typedef double D;
unsigned char B[1 << 17], I[S * S];
unsigned T[1 << 22], h[4], x = 0, z = -1, X, d, u, t[1024];
D E[S * S], w[6] = {.1, .1, .1, .1, .1, .1}, s[6], o[L], a[L], C[L][L], b[L], V[32], K[32];

D st(D p) { return log(p / (1 - p)); }
// P(value < v) under the Gaussian N(mu, sg), with the clipped tails folded into 0 and 15
D G(int v, D mu, D sg) { return v < 1 ? 0 : v > 15 ? 1 : erfc((mu + .5 - v) / sg) / 2; }

int main() {
  int i, j, k, y, c0, e, r, c, v, q, n = 0, l;
  d = fread(B, 1, 1 << 17, stdin) != 125000;
  F(i, 1024) t[i] = 16384 / (i + i + 3);
  F(i, 1 << 22) T[i] = 1u << 31;
  F(i, 32) V[i] = K[i] = 1;
  if (d) F(i, 4) X = X << 8 | B[u++];
  F(r, 500) F(c, 500) {
    k = (r + 4) * S + c + 4;
    if (!c) F(j, 4) I[k - j - 1] = I[k - S];  // replicate the left edge
    unsigned W = I[k - 1], N = I[k - S], NW = I[k - S - 1], NE = I[k - S + 1], WW = I[k - 2],
             NN = I[k - 2 * S];
    D m = 0, g, p, f;
    // OLS inputs: bias + causal neighbours within radius 4
    l = 0;
    o[l++] = 1;
    F(i, 5) for (j = -4; j < 5; j++) if (i || j < 0) o[l++] = I[k - i * S + j];
    F(i, L) m += a[i] * o[i];
    m = m < -.5 ? -.5 : m > 15.5 ? 15.5 : m;
    q = (fabs(E[k - 1]) + fabs(E[k - S]) + (fabs(E[k - S - 1]) + fabs(E[k - S + 1])) / 2) * 4;
    q = q > 31 ? 31 : q;
    g = sqrt(V[q] / K[q] * 2) + .03;
    h[0] = W << 8 | N << 4 | NE;
    h[1] = N << 8 | NN << 4 | NE | 4096;
    h[2] = W << 8 | WW << 4 | NW | 8192;
    h[3] = W << 20 | N << 16 | NW << 12 | NE << 8 | WW << 4 | NN;
    v = B[n / 2] >> (n & 1 ? 0 : 4) & 15;
    for (c0 = 1, e = 3; e >= 0; e--) {
      F(i, 4) s[i] = st(((T[((h[i] << 4 | c0) * 0x9E3779B1u + i * 0x7F4A7C15u) >> 10] >> 16) + .5) / 65536);
      j = (c0 << (e + 1)) - 16;  // values j .. j + 2^(e+1) - 1 remain; bit splits at j + 2^e
      p = G(j + (2 << e), m, g) - G(j, m, g);
      p = (G(j + (2 << e), m, g) - G(j + (1 << e), m, g)) / (p + 1e-30);
      s[4] = st(p < 1e-4 ? 1e-4 : p > 1 - 1e-4 ? 1 - 1e-4 : p);
      s[5] = .3;
      f = 0;
      F(i, 6) f += w[i] * s[i];
      p = 1 / (1 + exp(-f));
      // binary arithmetic coder
      unsigned P = p * 65536, M;
      P = P < 32 ? 32 : P > 65504 ? 65504 : P;
      M = x + ((unsigned long long)(z - x) * P >> 16);
      y = d ? X <= M : v >> e & 1;
      y ? (z = M) : (x = M + 1);
      while (!((x ^ z) >> 24)) {
        if (d) X = X << 8 | B[u++]; else putchar(z >> 24);
        x <<= 8, z = z << 8 | 255;
      }
      F(i, 6) w[i] += .001 * (y - p) * s[i];
      F(i, 4) {
        unsigned *A = &T[((h[i] << 4 | c0) * 0x9E3779B1u + i * 0x7F4A7C15u) >> 10], M = *A & 1023;
        if (M < 255) ++*A;
        *A += ((int)((y << 22) - (*A >> 10)) >> 3) * t[M] & ~1023u;
      }
      c0 = c0 * 2 + y;
    }
    I[k] = c0 -= 16;
    if (!r) for (j = k; j < 5 * S; j++) F(i, 4) I[j - i * S - S] = c0;  // top edge
    if (c > 498) F(j, 8) I[k + j] = c0;                                // right edge
    E[k] = c0 - m;
    V[q] += E[k] * E[k], K[q]++;
    if (K[q] > 500) V[q] /= 2, K[q] /= 2;
    // recursive least squares (forgetting .999) + one Gauss-Seidel sweep for the weights
    F(i, L) {
      F(j, L) C[i][j] = .999 * C[i][j] + .001 * o[i] * o[j];
      b[i] = .999 * b[i] + .001 * o[i] * c0;
      f = b[i];
      F(j, L) if (j != i) f -= C[i][j] * a[j];
      a[i] = f / (C[i][i] + .01);
    }
    n++;
  }
  if (d) F(i, 125000) putchar(I[(i / 250 + 4) * S + i % 250 * 2 + 4] << 4 | I[(i / 250 + 4) * S + i % 250 * 2 + 5]);
  else putchar((x >> 24) + 1);
  return 0;
}
