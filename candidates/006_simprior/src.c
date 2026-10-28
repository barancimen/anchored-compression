// 006: source-simulation prior on the raster reference (readable source; p.c = this file
// minified by lab/minify2.py). Not an anchor method: it tests "move the learning out of the
// archive" by simulating the source. Both directions first generate 3 images with a C port of
// gen_image.py (fixed seed, own PRNG; see lab/sim.c for the readable port) and run the model
// over them in training mode, then code the input with the trained state.
// The model, scan order and coder are 000_raster_ref's, unchanged: scan position (r, l) codes
// image pixel (row 499 - l, column r), i.e. the image is coded transposed and mirrored.
// One program: stdin of exactly 125000 bytes = image -> compress to stdout; otherwise the
// stdin is an archive -> decompress to stdout. Each 4-bit pixel is coded MSB-first as 4
// binary decisions. Each decision mixes (logistic mixer, one shared weight set):
//   5 adaptive counters in hashed contexts (W,WW,NEE) (N,NE,NNE) (W,N,NW) (N,NE,NNE,W,N,NW)
//   (WWW,NEEE,NE), and
//   a discretized-Gaussian input centred on an online least-squares (RLS, forgetting .9997)
//   prediction from the 112 causal neighbours within radius 7 (one Gauss-Seidel sweep per
//   pixel), with a variance learned per bucket of nearby prediction errors.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define F(i, n) for (i = 0; i < n; i++)
#define S 528  // row stride of the padded image I: 8-pixel left margin, 20-pixel right margin
#define L 112  // RLS inputs: causal neighbours within radius 7
typedef double D;
typedef double _Complex Q;
unsigned char B[1 << 17], O[125000], I[S * S], N[1 << 22], K[250000];  // K: simulated image
unsigned x, z = -1, X, M;
unsigned long long R = 1;  // simulator PRNG state (fixed seed)
// h: contexts; H: counter slots (H[5] stays 0: a dummy slot so one loop updates everything)
int h[5], H[6], d, u, i, j, k, y, c, e, r, l, v, q, b, n;
// P: counter probabilities (N = update counts, rate 1/(n+1.5), n <= 80); E: |RLS error|;
// C: exponentially weighted RLS sums, column L holds the right-hand side (o[L] = pixel and
// a[L] = -1 fold it into the Gauss-Seidel sum); V: 2 x error variance per bucket
D P[1 << 22], E[S * S], w[6], s[6], o[L + 1], a[L + 1], C[L][L + 1], V[32], m, g, t, f, p;
// simulator: W one random field, A the image being built, J a sorted copy; Z the FFT grid
D W[250000], A[250000], J[250000];
Q Z[1 << 18];

D st(D p) { return p = p * .9998 + .0001, log(p / (1 - p)); }
// P(pixel < (v << e) - 16) under N(m, g^2 / 2); clipped tails folded into 0 and 15
D G(int v) { return v <<= e, v < 17 ? 0 : v > 31 ? 1 : erfc((m + 16.5 - v) / g) / 2; }

// ---- source simulator (C port of gen_image.py) ----
D U() { return ((R = R * 6364136223846793005ull + 1) >> 11) * 0x1p-53; }  // uniform [0, 1)
D Y() { return sqrt(-2 * log(1 - U())) * cos(6.283185 * U()); }            // N(0, 1)
// in-place FFT of 512 points with stride s (bit reversal, then butterflies)
void T(Q *a, int s) {
  int i, j = 0, k, n;
  Q t;
  for (i = 1; i < 512; i++) {
    for (k = 256; j & k; k /= 2) j ^= k;
    if (i < (j ^= k)) t = a[i * s], a[i * s] = a[j * s], a[j * s] = t;
  }
  for (n = 1; n < 512; n *= 2)
    for (i = 0; i < 512; i += 2 * n)
      for (j = i; j < i + n; j++)
        t = a[(j + n) * s] * (cos(M_PI * (j - i) / n) + 1i * sin(M_PI * (j - i) / n)), a[(j + n) * s] = a[j * s] - t, a[j * s] += t;
}
// unit-variance Gaussian random field with power spectrum ~ 1/f^b on a 512x512 grid,
// cropped to 500x500, into W
void Gf(D b) {
  F(i, 1 << 18) {
    D fx = i % 512 < 256 ? i % 512 : i % 512 - 512, fy = i / 512 < 256 ? i / 512 : i / 512 - 512;
    Z[i] = fx || fy ? (Y() + 1i * Y()) * pow(hypot(fx, fy), -b / 2) : 0;
  }
  F(i, 512) T(Z + i * 512, 1);
  F(i, 512) T(Z + i, 512);
  t = f = 0;
  F(i, 250000) W[i] = __real__ Z[i / 500 * 512 + i % 500], t += W[i], f += W[i] * W[i];
  t /= 250000, f = sqrt(f / 250000 - t * t);
  F(i, 250000) W[i] = (W[i] - t) / f;
}
int cm(const void *a, const void *b) { return *(D *)a > *(D *)b ? 1 : *(D *)a < *(D *)b ? -1 : 0; }
// one simulated image into K: three random fields, a light gradient, 30-44 objects (disc,
// ellipse, rectangle, ring, line, striped disc; half of the non-striped ones opaque), then
// scaled by the 0.5 / 99.5 percentiles, sensor noise .28 and quantized to 0..15
void sim() {
  D pk, cx, cy, ra, co, si, bb, ww, lg, pe, am, lv, uu, vv, dd, mm;
  int kd, op;
  Gf(3.4);
  F(i, 250000) A[i] = W[i];
  Gf(2.2);
  F(i, 250000) A[i] += .3 * W[i];
  Gf(1);
  F(i, 250000) A[i] += .08 * W[i];
  t = Y() * .6, f = Y() * .6;
  F(i, 250000) A[i] += t * (i % 500 - 250) / 500 + f * (i / 500 - 250) / 500;
  for (j = 30 + U() * 15; j--;) {
    pk = U(), cx = U() * 500, cy = U() * 500, ra = 8 + U() * 52, dd = U() * M_PI, co = cos(dd), si = sin(dd);
    kd = (pk > .22) + (pk > .4) + (pk > .62) + (pk > .72) + (pk > .86);
    bb = ra * (kd - 1 ? .3 + U() * .7 : .3 + U() * .6);
    ww = kd - 3 ? .8 + U() * 2.2 : 2 + U() * 4;
    lg = 60 + U() * 240, pe = 5 + U() * 13, am = Y() * 1.1, op = kd < 5 && U() < .5, lv = Y() * 1.2;
    if (op) Gf(3);
    F(i, 250000) {
      uu = (i % 500 - cx) * co + (i / 500 - cy) * si, vv = (i / 500 - cy) * co - (i % 500 - cx) * si, dd = hypot(uu, vv);
      dd = kd == 1 ? (hypot(uu / ra, vv / bb) - 1) * bb : kd == 2 ? fmax(fabs(uu) - ra, fabs(vv) - bb) : kd == 3 ? fabs(dd - ra) - ww
         : kd == 4 ? fmax(fabs(vv) - ww, fabs(uu) - lg) : dd - ra;
      mm = fmin(fmax(.5 - dd, 0), 1) * (kd > 4 ? sin(2 * M_PI * uu / pe) : 1);
      A[i] = op ? A[i] * (1 - mm) + (lv + .25 * W[i]) * mm : A[i] + am * mm;
    }
  }
  F(i, 250000) J[i] = A[i];
  qsort(J, 250000, 8, cm);
  t = J[1250], f = J[248750];
  F(i, 250000) p = (A[i] - t) / (f - t) * 16 + Y() * .28, K[i] = p < 0 ? 0 : p > 15 ? 15 : p;
}

int main() {
  d = fread(B, 1, 1 << 17, stdin) != 125000;
  F(i, 1 << 22) P[i] = .5, V[i % 32] = 1, a[L] = -1;
  F(i, 4) X = X << 8 | B[u++];
  // n = 0, 1, 2: train on a simulated image; n = 3: code the input
  F(n, 4) {
  if (n < 3) sim();
  F(i, S * S) I[i] = 8;
  F(r, 500) F(l, 500) {
    k = (r + 8) * S + l + 8;
    if (!l) F(j, 8) I[k - j - 1] = I[k - S];  // replicate the left edge
    m = 0;
    // the 112 causal neighbours = first 112 cells (row-major) of the 8x15 block ending at k
    F(e, L) m += a[e] * (o[e] = I[k - 7 * S - 7 + e / 15 * S + e % 15]);
    q = fmin((E[k - 1] + E[k - S] + E[k - S - 1] + E[k - S + 1]) * 3, 31);
    g = sqrt(V[q] * 2);
    h[0] = I[k - 1] << 8 | I[k - 2] << 4 | I[k - S + 2];
    h[1] = I[k - S] << 8 | I[k - S + 1] << 4 | I[k - 2 * S + 1];
    h[2] = I[k - 1] << 8 | I[k - S] << 4 | I[k - S - 1];
    h[3] = h[1] << 12 | h[2];
    h[4] = I[k - 3] << 8 | I[k - S + 3] << 4 | I[k - S + 1];
    b = (499 - l) * 500 + r;
    v = n < 3 ? K[b] : B[b / 2] >> 4 - b % 2 * 4 & 15;
    // c = 2 x (bits so far, with a leading 1); the remaining values are (c << e) - 16 ...
    for (c = 1, e = 4; e--; c += y) {
      c *= 2;
      t = G(c + 2), f = 0;
      s[5] = st((t - G(c + 1)) / (t - G(c) + 1e-30));
      F(i, 5) s[i] = st(P[H[i] = (h[i] << 5 | c | i << 29) * 2654435761u >> 10]);
      F(i, 6) f += w[i] * s[i];
      p = 1 / (1 + exp(-f));
      y = v >> e & 1;
      // binary arithmetic coder (carry-less, 32-bit range); not used while training
      if (n > 2) {
        M = x + (z - x) * (p * .9999 + .00005);
        d ? y = X <= M : 0;
        y ? z = M : (x = M + 1);
        for (; (x ^ z) < 1 << 24; x <<= 8, z = z << 8 | 255) d ? X = X << 8 | B[u++] : putchar(z >> 24);
      }
      F(i, 6) w[i] += .001 * (y - p) * s[i], j = H[i], N[j] += N[j] < 80, P[j] += (y - P[j]) / (N[j] + .5);
    }
    I[k] = c -= 16;
    O[b / 2] |= n / 3 * c << 4 - b % 2 * 4;  // decoded pixels only
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
  }
  }
  d ? fwrite(O, 1, 125000, stdout) : putchar((x >> 24) + 1);
}
