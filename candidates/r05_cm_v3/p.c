// v3: pruned CM. Inputs: 4 neighbour-context counters, a discretized-Gaussian input from
// an online OLS predictor, and a bias; one shared mixer weight set; APM on (node, mu).
// Same program compresses (stdin = 125000-byte image) and decompresses (stdin = archive).
#include <stdio.h>
#include <math.h>
#define NI 14
typedef unsigned U;
unsigned char B[1 << 17], I[504][504];
U T[1 << 22], h[4], x1, x2 = -1, X, d, ip, dt[1024];
float E[504][504], w[6], s[6], A[16 * 64 * 33];
double ox[NI], ow[NI], oc[NI][NI], ob[NI], oL[NI][NI], sv[32], sn[32];

// binary arithmetic coder: codes bit y (or decodes it when d) with P(1) = p / 65536
int code(int y, U p) {
  U m = x1 + ((unsigned long long)(x2 - x1) * p >> 16);
  if (d) y = X <= m;
  y ? (x2 = m) : (x1 = m + 1);
  while (!((x1 ^ x2) >> 24)) {
    if (d) X = X << 8 | B[ip++]; else putchar(x2 >> 24);
    x1 <<= 8, x2 = x2 << 8 | 255;
  }
  return y;
}
double st(double p) { return log(p / (1 - p)); }
double sq(double x) { return 1 / (1 + exp(-x)); }

// recursive least squares with forgetting: refit weights ow to predict v from ox
void ols(double v) {
  int i, j, k;
  double t, y[NI];
  for (i = 0; i < NI; i++) {
    for (j = 0; j < NI; j++) oc[i][j] = .998 * oc[i][j] + .002 * ox[i] * ox[j];
    ob[i] = .998 * ob[i] + .002 * ox[i] * v;
  }
  for (i = 0; i < NI; i++)  // Cholesky factorisation of oc + .001 I
    for (j = 0; j <= i; j++) {
      t = oc[i][j] + (i == j) * .001;
      for (k = 0; k < j; k++) t -= oL[i][k] * oL[j][k];
      if (i == j) { if (t <= 0) return; oL[i][i] = sqrt(t); }
      else oL[i][j] = t / oL[j][j];
    }
  for (i = 0; i < NI; i++) { t = ob[i]; for (k = 0; k < i; k++) t -= oL[i][k] * y[k]; y[i] = t / oL[i][i]; }
  for (i = NI; i--;) { t = y[i]; for (k = i + 1; k < NI; k++) t -= oL[k][i] * ow[k]; ow[i] = t / oL[i][i]; }
}

int main() {
  int i, j, y, c0, b, r, c, k, v, a, q, mq;
  d = fread(B, 1, 1 << 17, stdin) != 125000;
  for (i = 0; i < 1024; i++) dt[i] = 16384 / (i + i + 3);
  for (i = 0; i < 1 << 22; i++) T[i] = 1u << 31;
  for (i = 0; i < 16 * 64 * 33; i++) A[i] = sq((i % 33 - 16) / 2.);
  for (i = 0; i < 6; i++) w[i] = .1;
  for (i = 0; i < 32; i++) sv[i] = sn[i] = 1;
  if (d) for (i = 0; i < 4; i++) X = X << 8 | B[ip++];
  for (r = 2; r < 502; r++)
    for (c = 2; c < 502; c++) {
      if (c == 2) I[r][1] = I[r][0] = I[r - 1][2];  // replicate edges
      U W = I[r][c - 1], N = I[r - 1][c], NW = I[r - 1][c - 1], NE = I[r - 1][c + 1],
        WW = I[r][c - 2], NN = I[r - 2][c];
      k = (r - 2) * 500 + c - 2;
      v = B[k / 2] >> (k & 1 ? 0 : 4) & 15;
      double xs[NI] = {W, N, NW, NE, WW, NN, I[r - 2][c + 1], I[r - 2][c - 1], I[r - 1][c - 2],
                       I[r - 1][c + 2], I[r - 2][c + 2], I[r - 2][c - 2], I[r][c - 3], 1},
             mu = 0, sc;
      for (i = 0; i < NI; i++) ox[i] = xs[i], mu += ow[i] * xs[i];
      mu = mu < -.5 ? -.5 : mu > 15.5 ? 15.5 : mu;
      q = (E[r][c - 1] + E[r - 1][c] + (E[r - 1][c - 1] + E[r - 1][c + 1]) / 2) * 4;
      q = q > 31 ? 31 : q;
      sc = sqrt(sv[q] / sn[q]) + .02;
      mq = mu * 4 + 2;
      h[0] = W << 8 | N << 4 | NE;
      h[1] = N << 8 | NN << 4 | NE | 1 << 12;
      h[2] = W << 8 | WW << 4 | NW | 2 << 12;
      h[3] = (W << 20 | N << 16 | NW << 12 | NE << 8 | WW << 4 | NN) * 2654435761u >> 8;
      for (c0 = 1, b = 3; b >= 0; b--) {
        for (i = 0; i < 4; i++) s[i] = st(((T[(h[i] * 16 + c0) & 0x3fffff] >> 16) + .5) / 65536);
        // discretized Gaussian: P(bit = 1 | bits so far) from the CDF at bin edges
        int lo = (c0 << (b + 1)) - 16, md = lo + (1 << b), hi = md + (1 << b);
#define F(v) ((v) <= 0 ? 0 : (v) >= 16 ? 1 : .5 * erfc((mu + .5 - (v)) / sc / 1.4142))
        double pl = F(lo), pm = F(md), ph = F(hi), pp = (ph - pm) / (ph - pl + 1e-30);
        s[4] = st(pp < 1e-4 ? 1e-4 : pp > 1 - 1e-4 ? 1 - 1e-4 : pp);
        s[5] = .3;
        float t = 0, p, z, f, e;
        for (i = 0; i < 6; i++) t += w[i] * s[i];
        p = sq(t);
        z = (t < -7.99 ? -7.99 : t > 7.99 ? 7.99 : t) * 2 + 16;
        a = (c0 * 64 + (mq > 63 ? 63 : mq)) * 33 + (int)z;
        f = z - (int)z;
        U pi = (p + 3 * (A[a] * (1 - f) + A[a + 1] * f)) * 16384;
        y = code(v >> b & 1, pi < 32 ? 32 : pi > 65504 ? 65504 : pi);
        e = y - p;
        for (i = 0; i < 6; i++) w[i] += .001 * e * s[i];
        for (i = 0; i < 4; i++) {
          U *e2 = &T[(h[i] * 16 + c0) & 0x3fffff], m = *e2 & 1023;
          if (m < 255) ++*e2;
          *e2 += ((int)((y << 22) - (*e2 >> 10)) >> 3) * dt[m] & ~1023u;
        }
        A[a] += (y - A[a]) * (1 - f) / 32;
        A[a + 1] += (y - A[a + 1]) * f / 32;
        c0 = c0 * 2 + y;
      }
      I[r][c] = c0 -= 16;
      if (r == 2) for (j = c; j < 504; j++) I[1][j] = I[0][j] = c0;
      if (c == 501) I[r][502] = I[r][503] = c0;
      double er = c0 - mu;
      E[r][c] = fabs(er);
      sv[q] += er * er, sn[q]++;
      if (sn[q] > 500) sv[q] *= .5, sn[q] *= .5;
      ols(c0);
    }
  if (d) for (i = 0; i < 125000; i++) putchar(I[2 + i / 250][2 + i % 250 * 2] << 4 | I[2 + i / 250][3 + i % 250 * 2]);
  else putchar((x1 >> 24) + 1);
  return 0;
}
