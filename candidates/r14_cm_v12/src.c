// v12 (readable source; p.c = this file minified by lab/minify2.py). v11 golfed.
// One program: stdin of exactly 125000 bytes = image -> compress to stdout; otherwise the
// stdin is an archive -> decompress to stdout. Each 4-bit pixel is coded MSB-first as 4
// binary decisions. Each decision mixes (logistic mixer, one shared weight set):
//   5 adaptive counters in hashed contexts (W,WW,NEE) (N,NE,NN) (W,N,NW) (N,NE,NN,W,N,NW)
//   (WWW,NEEE,NE), and
//   a discretized-Gaussian input centred on an online least-squares (RLS, forgetting .9997)
//   prediction from the 84 causal neighbours within radius 6 (one Gauss-Seidel sweep per
//   pixel), with a variance learned per bucket of nearby prediction errors.
#include <stdio.h>
#include <math.h>
#define F(i, n) for (i = 0; i < n; i++)
#define S 528  // row stride of the padded image I: 8-pixel left margin, 20-pixel right margin
#define L 84   // RLS inputs: causal neighbours within radius 6
typedef double D;
unsigned char B[1 << 17], O[125000], I[S * S], N[1 << 22];
unsigned x, z = -1, X, M;
// h: contexts; H: counter slots (H[5] stays 0: a dummy slot so one loop updates everything)
int h[5], H[6], d, u, i, j, k, y, c, e, r, l, v, q, n;
// P: counter probabilities (N = update counts, rate 1/(n+1.5), n <= 60); E: |RLS error|;
// C: exponentially weighted RLS sums, column L holds the right-hand side (o[L] = pixel and
// a[L] = -1 fold it into the Gauss-Seidel sum); V: 2 x error variance per bucket
D P[1 << 22], E[S * S], w[6], s[6], o[L + 1], a[L + 1], C[L][L + 1], V[32], m, g, t, f, p;

D st(D p) { return p = p * .9998 + .0001, log(p / (1 - p)); }
// P(pixel < (v << e) - 16) under N(m, g^2 / 2); clipped tails folded into 0 and 15
D G(int v) { return v <<= e, v < 17 ? 0 : v > 31 ? 1 : erfc((m + 16.5 - v) / g) / 2; }

int main() {
  d = fread(B, 1, 1 << 17, stdin) != 125000;
  F(i, 1 << 22) P[i] = .5, I[i % (S * S)] = 8, V[i % 32] = 1, a[L] = -1;
  F(i, 4) X = X << 8 | B[u++];
  F(r, 500) F(l, 500) {
    k = (r + 8) * S + l + 8;
    if (!l) F(j, 8) I[k - j - 1] = I[k - S];  // replicate the left edge
    m = 0;
    // the 84 causal neighbours = first 84 cells (row-major) of the 7x13 block ending at k
    F(e, L) m += a[e] * (o[e] = I[k - 6 * S - 6 + e / 13 * S + e % 13]);
    q = fmin((E[k - 1] + E[k - S] + E[k - S - 1] + E[k - S + 1]) * 3, 31);
    g = sqrt(V[q] * 2);
    h[0] = I[k - 1] << 8 | I[k - 2] << 4 | I[k - S + 2];
    h[1] = I[k - S] << 8 | I[k - S + 1] << 4 | I[k - 2 * S];
    h[2] = I[k - 1] << 8 | I[k - S] << 4 | I[k - S - 1];
    h[3] = h[1] << 12 | h[2];
    h[4] = I[k - 3] << 8 | I[k - S + 3] << 4 | I[k - S + 1];
    v = B[n / 2] >> 4 - n % 2 * 4 & 15;
    // c = 2 x (bits so far, with a leading 1); the remaining values are (c << e) - 16 ...
    for (c = 1, e = 4; e--; c += y) {
      c *= 2;
      t = G(c + 2), f = 0;
      s[5] = st((t - G(c + 1)) / (t - G(c) + 1e-30));
      F(i, 5) s[i] = st(P[H[i] = (h[i] << 5 | c | i << 29) * 2654435761u >> 10]);
      F(i, 6) f += w[i] * s[i];
      p = 1 / (1 + exp(-f));
      // binary arithmetic coder (carry-less, 32-bit range)
      M = x + (z - x) * (p * .9999 + .00005);
      y = d ? X <= M : v >> e & 1;
      y ? z = M : (x = M + 1);
      for (; (x ^ z) < 1 << 24; x <<= 8, z = z << 8 | 255) d ? X = X << 8 | B[u++] : putchar(z >> 24);
      F(i, 6) w[i] += .001 * (y - p) * s[i], j = H[i], N[j] += N[j] < 60, P[j] += (y - P[j]) / (N[j] + .5);
    }
    I[k] = c -= 16;
    O[n / 2] |= c << 4 - n % 2 * 4;
    if (l > 498) F(j, 20) I[k + j] = c;  // replicate the right edge
    E[k] = fabs(c - m);
    V[q] += (E[k] * E[k] - V[q]) / 500;
    // RLS sums with forgetting .9997 (regulariser 100), then one Gauss-Seidel sweep
    o[L] = c;
    F(i, L) {
      m = 0;
      F(j, L + 1) C[i][j] = C[i][j] * .9997 + o[i] * o[j], m -= j - i ? C[i][j] * a[j] : 0;
      a[i] = m / (C[i][i] + 100);
    }
    n++;
  }
  d ? fwrite(O, 1, n / 2, stdout) : putchar((x >> 24) + 1);
}
