// Oracle probe on the 007 model (lab tool): candidates/007_template_ref/src.c with a known
// additive layer. Without -D options the archive is identical to 007's.
//   -DTMPL  read the layer from the file named by $LAYER (default tmpl.f32; float32 grey levels
//           per image pixel, row-major, as written by lab/gen_labeled.py); the RLS predicts the
//           image minus the layer and adds the layer back
//   -DNLMT  non-local means matches and averages the image minus the layer, plus the layer
//   -DCTXT  counter contexts read the image minus the rounded layer
// Build: cc -O2 [-D...] lab/objects/oracle.c -lm; run: ./a.out<image>archive
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define F(i, n) for (i = 0; i < n; i++)
#define S 528
#define L 112
#define T(i) (i < 10 ? (i / 5 - 2) * S + i % 5 - 2 : 9 - i)
typedef double D;
unsigned char B[1 << 17], O[125000], I[S * S], N[1 << 22];
unsigned x, z = -1, X, M;
int h[5], H[7], d, u, i, j, k, y, c, e, r, l, v, q[2], b, n;
D P[1 << 22], E[2][S * S], w[7], s[7], o[L + 1], a[L + 1], C[L][L + 1], V[64], Y[2], R[2], m, t, f, p, Zs, W;
float TL[250000];
D Z[S * S];  // image minus layer, padded like I
#if defined(TMPL) || defined(NLMT) || defined(CTXT)
#define TB TL[b]
#else
#define TB 0
#endif
#ifdef TMPL
#define ZR Z
#define TR TB
#else
#define ZR I
#define TR 0
#endif
#ifdef NLMT
#define ZN Z
#define TN TB
#else
#define ZN I
#define TN 0
#endif
#ifdef CTXT
unsigned char J[S * S];
#define IC J
#else
#define IC I
#endif

D st(D p) { return p = p * .9998 + .0001, log(p / (1 - p)); }
D G(int v, int j) { return v <<= e, v < 17 ? 0 : v > 31 ? 1 : erfc((Y[j] + 16.5 - v) / R[j]) / 2; }

int main() {
  d = fread(B, 1, 1 << 17, stdin) != 125000;
  F(i, 1 << 22) P[i] = .5, I[i % (S * S)] = 8, V[i % 64] = 1, a[L] = -1;
  F(i, S * S) Z[i] = 8;
#ifdef CTXT
  F(i, S * S) J[i] = 8;
#endif
#if defined(TMPL) || defined(NLMT) || defined(CTXT)
  {
    char *fn = getenv("LAYER");
    FILE *fp = fopen(fn ? fn : "tmpl.f32", "rb");
    if (!fp || fread(TL, 4, 250000, fp) != 250000) return fprintf(stderr, "no layer\n"), 1;
  }
#endif
  F(i, 4) X = X << 8 | B[u++];
  F(r, 500) F(l, 500) {
    k = (r + 8) * S + l + 8;
    b = (499 - l) * 500 + r;
    if (!l) F(j, 8) I[k - j - 1] = I[k - S], Z[k - j - 1] = Z[k - S];
#ifdef CTXT
    if (!l) F(j, 8) J[k - j - 1] = J[k - S];
#endif
    Y[0] = TR;
    F(e, L) Y[0] += a[e] * (o[e] = ZR[k - 7 * S - 7 + e / 15 * S + e % 15]);
    Zs = W = 0;
    for (n = -4 * S - 4; n < 0; n++)
      if ((n + 2116) % S < 9) {
        t = 0;
        F(i, 12) f = ZN[k + T(i)] - ZN[k + n + T(i)], t += f * f;
        Zs += t = exp(-t / 2), W += t * ZN[k + n];
      }
    Y[1] = W / (Zs + 1e-99) + TN;
    F(j, 2) q[j] = j * 32 + fmin((E[j][k - 1] + E[j][k - S] + E[j][k - S - 1] + E[j][k - S + 1]) * 3, 31), R[j] = sqrt(V[q[j]] * 2);
    h[0] = IC[k - 1] << 8 | IC[k - 2] << 4 | IC[k - S + 2];
    h[1] = IC[k - S] << 8 | IC[k - S + 1] << 4 | IC[k - 2 * S + 1];
    h[2] = IC[k - 1] << 8 | IC[k - S] << 4 | IC[k - S - 1];
    h[3] = h[1] << 12 | h[2];
    h[4] = IC[k - 3] << 8 | IC[k - S + 3] << 4 | IC[k - S + 1];
    v = B[b / 2] >> 4 - b % 2 * 4 & 15;
    for (c = 1, e = 4; e--; c += y) {
      c *= 2;
      f = 0;
      F(j, 2) t = G(c + 2, j), s[5 + j] = st((t - G(c + 1, j)) / (t - G(c, j) + 1e-30));
      F(i, 5) s[i] = st(P[H[i] = (h[i] << 5 | c | i << 29) * 2654435761u >> 10]);
      F(i, 7) f += w[i] * s[i];
      p = 1 / (1 + exp(-f));
      M = x + (z - x) * (p * .9999 + .00005);
      y = d ? X <= M : v >> e & 1;
      y ? z = M : (x = M + 1);
      for (; (x ^ z) < 1 << 24; x <<= 8, z = z << 8 | 255) d ? X = X << 8 | B[u++] : putchar(z >> 24);
      F(i, 7) w[i] += .001 * (y - p) * s[i], j = H[i], N[j] += N[j] < 80, P[j] += (y - P[j]) / (N[j] + .5);
    }
    I[k] = c -= 16;
    Z[k] = c - TB;
#ifdef CTXT
    J[k] = (int)fmin(fmax(floor(c - TB + .5), 0), 15);
#endif
    O[b / 2] |= c << 4 - b % 2 * 4;
    if (l > 498) F(j, 20) I[k + j] = c, Z[k + j] = Z[k];
#ifdef CTXT
    if (l > 498) F(j, 20) J[k + j] = J[k];
#endif
    F(j, 2) t = E[j][k] = fabs(c - Y[j]), V[q[j]] += (t * t - V[q[j]]) / 500;
    o[L] = c - TR;
    F(i, L) {
      m = 0;
      F(j, L + 1) C[i][j] = C[i][j] * .9997 + o[i] * o[j], m -= j - i ? C[i][j] * a[j] : 0;
      a[i] = m / (C[i][i] + 100);
    }
  }
  d ? fwrite(O, 1, 125000, stdout) : putchar((x >> 24) + 1);
}
