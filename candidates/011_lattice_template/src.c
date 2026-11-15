// 011: 003_lattice_lookahead_t128 plus the non-local means input of 007_template_ref (readable
// source; p.c = minified by lab/minify2.py). Family A with the best model: compare with 007.
// As in 003: Built from the raster reference (000_raster_ref): same scan order (scan
// position (r, l) codes image pixel (row 499 - l, column r)), same model and coder.
// Anchors: the scan positions with r % T == 0 and l % T == 0, stored verbatim (one byte each)
// at the start of the archive. All other cells of the image buffer start as the bilinear
// interpolation of the anchors, so the RLS stencil can read ahead of the scan: it covers the
// 112 causal neighbours within radius 7 plus the current cell, the 7 cells to its right and the
// 15 cells of the next scan row (fill estimates until they are coded). Anchors are not coded
// again (the model still learns from them).
// One program: stdin of exactly 125000 bytes = image -> compress to stdout; otherwise the
// stdin is an archive -> decompress to stdout. Each 4-bit pixel is coded MSB-first as 4
// binary decisions. Each decision mixes (logistic mixer, one shared weight set):
//   5 adaptive counters in hashed contexts (W,WW,NEE) (N,NE,NNE) (W,N,NW) (N,NE,NNE,W,N,NW)
//   (WWW,NEEE,NE), and
//   a discretized-Gaussian input centred on an online least-squares (RLS, forgetting .9997)
//   prediction from the stencil (one Gauss-Seidel sweep per pixel), with a variance learned
//   per bucket of nearby prediction errors, and
//   a second discretized-Gaussian input centred on a non-local means prediction: every
//   position up to 4 scan rows back and 4 columns either side (causal ones only) is weighted
//   by exp(-d / 2), d = the squared difference between its 12 nearest causal neighbours and
//   ours; with its own error buckets.
#include <stdio.h>
#include <math.h>
#define F(i, n) for (i = 0; i < n; i++)
#define S 528  // row stride of the padded image I: 8-pixel left margin, 20-pixel right margin
#define L 135  // RLS inputs: the first 135 cells (row-major) of the 9x15 block from (r-7, l-7)
#define T 128  // anchor lattice step (scan coordinates)
// offset of the i-th of the 12 nearest causal neighbours: 2 rows x 5 columns above, then W, WW
#define TT(i) (i < 10 ? (i / 5 - 2) * S + i % 5 - 2 : 9 - i)
typedef double D;
unsigned char B[1 << 17], O[125000], I[S * S], N[1 << 22];
unsigned x, z = -1, X, M;
// h: contexts; H: counter slots (H[5], H[6] stay 0: dummy slots so one loop updates all)
// A: current pixel is an anchor; n: candidate offset; q: error bucket of each Gaussian input
int h[5], H[7], d, u, i, j, k, y, c, e, r, l, v, q[2], b, A, n;
// P: counter probabilities (N = update counts, rate 1/(n+1.5), n <= 80); Y, R: mean and sd of
// the Gaussian inputs (0: RLS, 1: non-local means); E: their |errors|; C: exponentially
// weighted RLS sums, column L holds the right-hand side (o[L] = pixel and a[L] = -1 fold it
// into the Gauss-Seidel sum); V: 2 x error variance per bucket (0..31 RLS, 32..63 non-local
// means); Z, W: weight sum and weighted sum of the non-local means
D P[1 << 22], E[2][S * S], w[7], s[7], o[L + 1], a[L + 1], C[L][L + 1], V[64], Y[2], R[2], m, t, f, p, Z, W;

D st(D p) { return p = p * .9998 + .0001, log(p / (1 - p)); }
// P(pixel < (v << e) - 16) under N(Y[j], R[j]^2 / 2); clipped tails folded into 0 and 15
D G(int v, int j) { return v <<= e, v < 17 ? 0 : v > 31 ? 1 : erfc((Y[j] + 16.5 - v) / R[j]) / 2; }
// value of the anchor at lattice point (i, j), or of the last one in range
D Q(int i, int j) { return I[((i < 499 / T * T ? i : 499 / T * T) + 8) * S + (j < 499 / T * T ? j : 499 / T * T) + 8]; }

int main() {
  d = fread(B, 1, 1 << 17, stdin) != 125000;
  F(i, 1 << 22) P[i] = .5, I[i % (S * S)] = 8, V[i % 64] = 1, a[L] = -1;
  // anchors: read from the image and written out (compress), or read from the archive
  F(r, 500) F(l, 500) if (!(r % T | l % T)) b = (499 - l) * 500 + r,
      I[(r + 8) * S + l + 8] = v = d ? B[u++] : B[b / 2] >> 4 - b % 2 * 4 & 15, d || putchar(v);
  // bilinear fill of every other cell from the 4 surrounding anchors
  F(r, 500) F(l, 500) i = r / T * T, j = l / T * T, f = r % T / (D)T, p = l % T / (D)T,
      I[(r + 8) * S + l + 8] = (1 - f) * ((1 - p) * Q(i, j) + p * Q(i, j + T)) +
                               f * ((1 - p) * Q(i + T, j) + p * Q(i + T, j + T)) + .5;
  F(i, 4) X = X << 8 | B[u++];
  F(r, 500) F(l, 500) {
    k = (r + 8) * S + l + 8;
    A = !(r % T | l % T);
    if (!l) F(j, 8) I[k - j - 1] = I[k - S];  // replicate the left edge
    F(e, L) Y[0] += a[e] * (o[e] = I[k - 7 * S - 7 + e / 15 * S + e % 15]);
    // non-local means over the positions n = -4 rows .. this row, -4 .. +4 columns, before k
    Z = W = 0;
    for (n = -4 * S - 4; n < 0; n++)
      if ((n + 2116) % S < 9) {  // 2116 = 4 * S + 4
        t = 0;
        F(i, 12) f = I[k + TT(i)] - I[k + n + TT(i)], t += f * f;
        Z += t = exp(-t / 2), W += t * I[k + n];
      }
    Y[1] = W / (Z + 1e-99);
    F(j, 2) q[j] = j * 32 + fmin((E[j][k - 1] + E[j][k - S] + E[j][k - S - 1] + E[j][k - S + 1]) * 3, 31), R[j] = sqrt(V[q[j]] * 2);
    h[0] = I[k - 1] << 8 | I[k - 2] << 4 | I[k - S + 2];
    h[1] = I[k - S] << 8 | I[k - S + 1] << 4 | I[k - 2 * S + 1];
    h[2] = I[k - 1] << 8 | I[k - S] << 4 | I[k - S - 1];
    h[3] = h[1] << 12 | h[2];
    h[4] = I[k - 3] << 8 | I[k - S + 3] << 4 | I[k - S + 1];
    b = (499 - l) * 500 + r;
    v = A ? I[k] : B[b / 2] >> 4 - b % 2 * 4 & 15;
    // c = 2 x (bits so far, with a leading 1); the remaining values are (c << e) - 16 ...
    for (c = 1, e = 4; e--; c += y) {
      c *= 2;
      f = 0;
      F(j, 2) t = G(c + 2, j), s[5 + j] = st((t - G(c + 1, j)) / (t - G(c, j) + 1e-30));
      F(i, 5) s[i] = st(P[H[i] = (h[i] << 5 | c | i << 29) * 2654435761u >> 10]);
      F(i, 7) f += w[i] * s[i];
      p = 1 / (1 + exp(-f));
      y = v >> e & 1;
      // binary arithmetic coder (carry-less, 32-bit range); anchors are not coded
      if (!A) {
        M = x + (z - x) * (p * .9999 + .00005);
        d ? y = X <= M : 0;
        y ? z = M : (x = M + 1);
        for (; (x ^ z) < 1 << 24; x <<= 8, z = z << 8 | 255) d ? X = X << 8 | B[u++] : putchar(z >> 24);
      }
      F(i, 7) w[i] += .001 * (y - p) * s[i], j = H[i], N[j] += N[j] < 80, P[j] += (y - P[j]) / (N[j] + .5);
    }
    I[k] = c -= 16;
    O[b / 2] |= c << 4 - b % 2 * 4;
    if (l > 498) F(j, 20) I[k + j] = c;  // replicate the right edge
    F(j, 2) t = E[j][k] = fabs(c - Y[j]), V[q[j]] += (t * t - V[q[j]]) / 500;
    Y[0] = 0;
    // RLS sums with forgetting .9997 (regulariser 100), then one Gauss-Seidel sweep
    o[L] = c;
    F(i, L) {
      m = 0;
      F(j, L + 1) C[i][j] = C[i][j] * .9997 + o[i] * o[j], m -= j - i ? C[i][j] * a[j] : 0;
      a[i] = m / (C[i][i] + 100);
    }
  }
  d ? fwrite(O, 1, 125000, stdout) : putchar((x >> 24) + 1);
}
