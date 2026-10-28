// Priming probe: how much of the archive is the model's learning cost?
// The model is the raster reference (candidates/000_raster_ref/src.c: same scan order, model,
// constants). It first runs over the priming images in training mode (no output, same
// updates as when coding), then codes the target with the trained state.
// Priming with real images from gen_image.py measures what a perfect source simulator inside
// the decoder could save at most.
//
// usage: prime TARGET [PRIME ...]   prints the ideal cost of every image in order
//        (bytes, under the state left by the images before it) and the target's archive size
// knobs: -DLIM (counter limit), -DLAM (RLS forgetting), -DLR (mixer rate)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef LIM
#define LIM 80
#endif
#ifndef LAM
#define LAM .9997
#endif
#ifndef LR
#define LR .001
#endif
#define F(i, n) for (i = 0; i < n; i++)
#define S 528
#define L 112
typedef double D;
unsigned char B[125000], I[S * S], N[1 << 22], out[1 << 17];
unsigned x, z, M;
int h[5], H[6], i, j, k, y, c, e, r, l, v, q, b, op;
D P[1 << 22], E[S * S], w[6], s[6], o[L + 1], a[L + 1], C[L][L + 1], V[32], m, g, t, f, p, bits;

D st(D p) { return p = p * .9998 + .0001, log(p / (1 - p)); }
D G(int v) { return v <<= e, v < 17 ? 0 : v > 31 ? 1 : erfc((m + 16.5 - v) / g) / 2; }

// one pass over the image in B; code = 1: arithmetic-code it into out (op = bytes)
void pass(int code) {
  F(i, S * S) I[i] = 8;
  x = 0, z = -1, op = 0, bits = 0;
  F(r, 500) F(l, 500) {
    k = (r + 8) * S + l + 8;
    if (!l) F(j, 8) I[k - j - 1] = I[k - S];
    m = 0;
    F(e, L) m += a[e] * (o[e] = I[k - 7 * S - 7 + e / 15 * S + e % 15]);
    q = fmin((E[k - 1] + E[k - S] + E[k - S - 1] + E[k - S + 1]) * 3, 31);
    g = sqrt(V[q] * 2);
    h[0] = I[k - 1] << 8 | I[k - 2] << 4 | I[k - S + 2];
    h[1] = I[k - S] << 8 | I[k - S + 1] << 4 | I[k - 2 * S + 1];
    h[2] = I[k - 1] << 8 | I[k - S] << 4 | I[k - S - 1];
    h[3] = h[1] << 12 | h[2];
    h[4] = I[k - 3] << 8 | I[k - S + 3] << 4 | I[k - S + 1];
    b = (499 - l) * 500 + r;
    v = B[b / 2] >> 4 - b % 2 * 4 & 15;
    for (c = 1, e = 4; e--; c += y) {
      c *= 2;
      t = G(c + 2), f = 0;
      s[5] = st((t - G(c + 1)) / (t - G(c) + 1e-30));
      F(i, 5) s[i] = st(P[H[i] = (h[i] << 5 | c | i << 29) * 2654435761u >> 10]);
      F(i, 6) f += w[i] * s[i];
      p = 1 / (1 + exp(-f));
      y = v >> e & 1;
      bits -= log2(y ? p * .9999 + .00005 : 1 - (p * .9999 + .00005));
      if (code) {
        M = x + (z - x) * (p * .9999 + .00005);
        y ? z = M : (x = M + 1);
        for (; (x ^ z) < 1 << 24; x <<= 8, z = z << 8 | 255) out[op++] = z >> 24;
      }
      F(i, 6) w[i] += LR * (y - p) * s[i], j = H[i], N[j] += N[j] < LIM, P[j] += (y - P[j]) / (N[j] + .5);
    }
    I[k] = c -= 16;
    if (l > 498) F(j, 20) I[k + j] = c;
    E[k] = fabs(c - m);
    V[q] += (E[k] * E[k] - V[q]) / 500;
    o[L] = c;
    F(i, L) {
      m = 0;
      F(j, L + 1) C[i][j] = C[i][j] * LAM + o[i] * o[j], m -= j - i ? C[i][j] * a[j] : 0;
      a[i] = m / (C[i][i] + 100);
    }
  }
  if (code) out[op++] = (x >> 24) + 1;
}

void load(const char *fn) {
  FILE *fp = fopen(fn, "rb");
  if (!fp || fread(B, 1, 125000, fp) != 125000) fprintf(stderr, "bad image %s\n", fn), exit(1);
  fclose(fp);
}

int main(int argc, char **argv) {
  F(i, 1 << 22) P[i] = .5;
  F(i, 32) V[i] = 1;
  a[L] = -1;
  for (int n = 2; n < argc; n++) load(argv[n]), pass(0), printf("prime %-28s %8.0f\n", argv[n], bits / 8);
  load(argv[1]), pass(1);
  printf("target %-27s %8.0f  archive %d\n", argv[1], bits / 8, op);
  return 0;
}
