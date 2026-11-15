// Cost map and probes on the raster reference (lab tool). The model and coder are
// candidates/000_raster_ref/src.c unchanged; without -D options the archive is identical.
// Compressing also writes, per image pixel (float32, row-major, for lab/where_bits.py):
//   cost.f32   bits the pixel cost
//   gcost.f32  bits under the Gaussian input alone (no counters, no mixer)
//   pred.f32   the RLS prediction, var.f32 the Gaussian's variance
// Build: cc -O2 [-D...] lab/costmap.c -lm; run: ./a.out<image>archive. Probes:
//   -DVQ=.167  Gaussian variance V - 1/12 (V is measured on rounded pixels and the Gaussian is
//              rounded again, so V alone counts the rounding twice)
//   -DGSN=4    4 Gauss-Seidel sweeps of the RLS system per pixel instead of 1
//   -DROB=4    robust RLS: each update weighted by 1 / (1 + e^2 / (ROB V))
//   -DTMPL     oracle: reads tmpl.f32 (lab/gen_labeled.py, the true additive objects in grey
//              levels); the RLS predicts the image minus that layer and adds it back
//   -DCLS      oracle: reads cls.u8 (a region class 0..2 per pixel) and keeps one set of RLS
//              weights per class
//
// The reference's own description follows.
// v15 (readable source; p.c = this file minified by lab/minify2.py). v12 with the scan
// order chosen per image: scan position (r, l) codes image pixel (row 499 - l, column r),
// i.e. the image is coded transposed and mirrored (best of the 8 orientations).
// One program: stdin of exactly 125000 bytes = image -> compress to stdout; otherwise the
// stdin is an archive -> decompress to stdout. Each 4-bit pixel is coded MSB-first as 4
// binary decisions. Each decision mixes (logistic mixer, one shared weight set):
//   5 adaptive counters in hashed contexts (W,WW,NEE) (N,NE,NNE) (W,N,NW) (N,NE,NNE,W,N,NW)
//   (WWW,NEEE,NE), and
//   a discretized-Gaussian input centred on an online least-squares (RLS, forgetting .9997)
//   prediction from the 112 causal neighbours within radius 7 (one Gauss-Seidel sweep per
//   pixel), with a variance learned per bucket of nearby prediction errors.
#include <stdio.h>
#include <math.h>
#define F(i, n) for (i = 0; i < n; i++)
#define S 528  // row stride of the padded image I: 8-pixel left margin, 20-pixel right margin
#define L 112  // RLS inputs: causal neighbours within radius 7
typedef double D;
unsigned char B[1 << 17], O[125000], I[S * S], N[1 << 22];
unsigned x, z = -1, X, M;
// h: contexts; H: counter slots (H[5] stays 0: a dummy slot so one loop updates everything)
int h[5], H[6], d, u, i, j, k, y, c, e, r, l, v, q, b;
// P: counter probabilities (N = update counts, rate 1/(n+1.5), n <= 80); E: |RLS error|;
// C: exponentially weighted RLS sums, column L holds the right-hand side (o[L] = pixel and
// a[L] = -1 fold it into the Gauss-Seidel sum); V: 2 x error variance per bucket
D P[1 << 22], E[S * S], w[6], s[6], o[L + 1], a0[3][L + 1], C0[3][L][L + 1], V[32], m, g, t, f, p;
D *a, (*C)[L + 1];  // the RLS weights and sums in use (one set unless -DCLS)
#ifdef CLS
unsigned char K[250000];
#endif

float Q[250000], Q2[250000], PM[250000], PV[250000];  // also the prediction and variance used
D pg;
#ifdef TMPL
float T[250000];
D Z[S * S];  // image minus the additive layer, padded like I
#define ZI Z
#define TB T[b]
#else
#define ZI I
#define TB 0
#endif

D st(D p) { return p = p * .9998 + .0001, log(p / (1 - p)); }
// P(pixel < (v << e) - 16) under N(m, g^2 / 2); clipped tails folded into 0 and 15
D G(int v) { return v <<= e, v < 17 ? 0 : v > 31 ? 1 : erfc((m + 16.5 - v) / g) / 2; }

int main() {
  d = fread(B, 1, 1 << 17, stdin) != 125000;
  F(i, 1 << 22) P[i] = .5, I[i % (S * S)] = 8, V[i % 32] = 1;
  F(i, 3) a0[i][L] = -1;
#ifdef CLS
  fread(K, 1, 250000, fopen("cls.u8", "rb"));
#endif
#ifdef TMPL
  fread(T, 4, 250000, fopen("tmpl.f32", "rb"));
  F(i, S * S) Z[i] = 8;
#endif
  F(i, 4) X = X << 8 | B[u++];
  F(r, 500) F(l, 500) {
    k = (r + 8) * S + l + 8;
    if (!l) F(j, 8) I[k - j - 1] = I[k - S], ZI[k - j - 1] = ZI[k - S];  // replicate the left edge
    b = (499 - l) * 500 + r;
#ifdef CLS
    a = a0[K[b]], C = C0[K[b]];
#else
    a = a0[0], C = C0[0];
#endif
    m = TB;
    // the 112 causal neighbours = first 112 cells (row-major) of the 8x15 block ending at k
    F(e, L) m += a[e] * (o[e] = ZI[k - 7 * S - 7 + e / 15 * S + e % 15]);
    q = fmin((E[k - 1] + E[k - S] + E[k - S - 1] + E[k - S + 1]) * 3, 31);
#ifdef VQ  // probe: the Gaussian's variance minus the rounding variance (VQ = 2/12), floored at the noise
    g = sqrt(fmax(V[q] * 2 - VQ, .157));
#else
    g = sqrt(V[q] * 2);
#endif
    h[0] = I[k - 1] << 8 | I[k - 2] << 4 | I[k - S + 2];
    h[1] = I[k - S] << 8 | I[k - S + 1] << 4 | I[k - 2 * S + 1];
    h[2] = I[k - 1] << 8 | I[k - S] << 4 | I[k - S - 1];
    h[3] = h[1] << 12 | h[2];
    h[4] = I[k - 3] << 8 | I[k - S + 3] << 4 | I[k - S + 1];
    b = (499 - l) * 500 + r;
    v = B[b / 2] >> 4 - b % 2 * 4 & 15;
    PM[b] = m, PV[b] = g * g / 2;
    // c = 2 x (bits so far, with a leading 1); the remaining values are (c << e) - 16 ...
    for (c = 1, e = 4; e--; c += y) {
      c *= 2;
      t = G(c + 2), f = 0;
      s[5] = st(pg = (t - G(c + 1)) / (t - G(c) + 1e-30));
      F(i, 5) s[i] = st(P[H[i] = (h[i] << 5 | c | i << 29) * 2654435761u >> 10]);
      F(i, 6) f += w[i] * s[i];
      p = 1 / (1 + exp(-f));
      // binary arithmetic coder (carry-less, 32-bit range)
      M = x + (z - x) * (p * .9999 + .00005);
      y = d ? X <= M : v >> e & 1;
      Q[b] -= log2(y ? p * .9999 + .00005 : 1 - (p * .9999 + .00005));
      Q2[b] -= log2(y ? pg * .9999 + .00005 : 1 - (pg * .9999 + .00005));
      y ? z = M : (x = M + 1);
      for (; (x ^ z) < 1 << 24; x <<= 8, z = z << 8 | 255) d ? X = X << 8 | B[u++] : putchar(z >> 24);
      F(i, 6) w[i] += .001 * (y - p) * s[i], j = H[i], N[j] += N[j] < 80, P[j] += (y - P[j]) / (N[j] + .5);
    }
    I[k] = c -= 16;
    O[b / 2] |= c << 4 - b % 2 * 4;
#ifdef TMPL
    Z[k] = c - TB;
#endif
    if (l > 498) F(j, 20) I[k + j] = c, ZI[k + j] = ZI[k];  // replicate the right edge
    E[k] = fabs(c - m);
    V[q] += (E[k] * E[k] - V[q]) / 500;
    // RLS sums with forgetting .9997 (regulariser 100), then one Gauss-Seidel sweep
    o[L] = c - TB;
#ifdef ROB  // probe: robust RLS, each update weighted by 1 / (1 + e^2 / (ROB V)) (Cauchy weight)
    pg = sqrt(1 / (1 + E[k] * E[k] / (ROB * V[q])));
    F(i, L + 1) o[i] *= pg;
#endif
    F(i, L) {
      m = 0;
      F(j, L + 1) C[i][j] = C[i][j] * .9997 + o[i] * o[j], m -= j - i ? C[i][j] * a[j] : 0;
      a[i] = m / (C[i][i] + 100);
    }
#ifdef GSN  // probe: GSN Gauss-Seidel sweeps per pixel instead of one
    for (int n = 1; n < GSN; n++) F(i, L) {
      m = 0;
      F(j, L + 1) m -= j - i ? C[i][j] * a[j] : 0;
      a[i] = m / (C[i][i] + 100);
    }
#endif
  }
  d ? fwrite(O, 1, 125000, stdout) : putchar((x >> 24) + 1);
  d || fwrite(Q, 4, 250000, fopen("cost.f32", "wb"));
  d || fwrite(Q2, 4, 250000, fopen("gcost.f32", "wb"));
  d || fwrite(PM, 4, 250000, fopen("pred.f32", "wb"));
  d || fwrite(PV, 4, 250000, fopen("var.f32", "wb"));
}
