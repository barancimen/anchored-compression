// v2: v1 + online OLS predictor (latent estimate mu) -> discretized-logistic input,
// plus mu-keyed contexts. Same program compresses/decompresses (stdin size decides).
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned U;
typedef unsigned long long Q;
#ifndef NC
#define NC 10
#endif
#ifndef LAM
#define LAM .998
#endif
#ifndef NU
#define NU .001
#endif
#ifndef MLR
#define MLR .02
#endif
#define NI 14
#ifndef WI
#define WI .2
#endif
#ifndef MASK
#define MASK 0xfff
#endif
#ifndef SCF
#define SCF 1
#endif
unsigned char B[1 << 17], I[504][504];
float E[504][504];  // |error| of OLS prediction per pixel
U *T[NC], M[NC], h[NC], x1, x2 = -1, X, d, ip, dt[1024];
float w[256][NC + 2], s[NC + 2];
float A[16 * 64 * 33];
double cost[NC+4]; float PC[250000];
double ox[NI], ow[NI], oc[NI][NI], ob[NI], oL[NI][NI], sv[32], sn[32];

int code(int y, U p) {
  U m = x1 + ((Q)(x2 - x1) * p >> 16);
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

void ols(double v) {
  int i, j, k;
  for (i = 0; i < NI; i++) {
    for (j = 0; j < NI; j++) oc[i][j] = LAM * oc[i][j] + (1 - LAM) * ox[i] * ox[j];
    ob[i] = LAM * ob[i] + (1 - LAM) * ox[i] * v;
  }
  for (i = 0; i < NI; i++)
    for (j = 0; j <= i; j++) {
      double t = oc[i][j] + (i == j) * NU;
      for (k = 0; k < j; k++) t -= oL[i][k] * oL[j][k];
      if (i == j) { if (t <= 0) return; oL[i][i] = sqrt(t); }
      else oL[i][j] = t / oL[j][j];
    }
  double y[NI];
  for (i = 0; i < NI; i++) { double t = ob[i]; for (k = 0; k < i; k++) t -= oL[i][k] * y[k]; y[i] = t / oL[i][i]; }
  for (i = NI - 1; i >= 0; i--) { double t = y[i]; for (k = i + 1; k < NI; k++) t -= oL[k][i] * ow[k]; ow[i] = t / oL[i][i]; }
}

int main() {
  int n = fread(B, 1, 1 << 17, stdin), i, j, y, c0, b;
  d = n != 125000;
  for (i = 0; i < 1024; i++) dt[i] = 16384 / (i + i + 3);
  for (i = 0; i < NC; i++) {
    int z = i == 7 ? 22 : 20;
    M[i] = (1 << z) - 1;
    T[i] = malloc(4 << z);
    for (j = 0; j < 1 << z; j++) T[i][j] = 1u << 31;
  }
  for (i = 0; i < 16 * 64 * 33; i++) A[i] = sq((i % 33 - 16) / 2.0);
  for (i = 0; i < 256; i++) for (j = 0; j < NC + 2; j++) w[i][j] = WI;
  for (i = 0; i < 32; i++) sv[i] = 1, sn[i] = 1;
  if (d) for (i = 0; i < 4; i++) X = X << 8 | B[ip++];
  for (int r = 2; r < 502; r++)
    for (int c = 2; c < 502; c++) {
      U W = I[r][c - 1], N = I[r - 1][c], NW = I[r - 1][c - 1], NE = I[r - 1][c + 1],
        WW = I[r][c - 2], NN = I[r - 2][c], NNE = I[r - 2][c + 1];
      int k = (r - 2) * 500 + c - 2, v = d ? 0 : (B[k / 2] >> (k & 1 ? 0 : 4)) & 15;
      // OLS inputs: causal neighbourhood + bias
      double xs[NI] = {W, N, NW, NE, WW, NN, NNE, I[r - 2][c - 1], I[r - 1][c - 2], I[r - 1][c + 2],
                       I[r - 2][c + 2], I[r - 2][c - 2], I[r][c - 3] , 1};
      double mu = 0;
      for (i = 0; i < NI; i++) ox[i] = xs[i], mu += ow[i] * xs[i];
      mu = mu < -.5 ? -.5 : mu > 15.5 ? 15.5 : mu;
      float ec = E[r][c - 1] + E[r - 1][c] + (E[r - 1][c - 1] + E[r - 1][c + 1]) / 2;
      int q = ec * 4; q = q > 31 ? 31 : q;
      #ifdef GAUSS
      double sc = sqrt(sv[q] / sn[q]) * SCF + .02;
#else
      double sc = sqrt(sv[q] / sn[q] * 3) / 3.1416 * SCF + .02;
#endif
      int mq = mu * 4 + 2;
      h[0] = 0;
      h[1] = W << 4 | N;
      h[2] = W << 8 | N << 4 | NW;
      h[3] = W << 8 | N << 4 | NE;
      h[4] = W << 12 | N << 8 | NW << 4 | NE;
      h[5] = N << 8 | NN << 4 | NE;
      h[6] = W << 8 | WW << 4 | NW;
      h[7] = (W << 20 | N << 16 | NW << 12 | NE << 8 | WW << 4 | NN) * 2654435761u >> 14;
      h[8] = mq << 5 | q;
      h[9] = mq << 8 | W << 4 | N;
      for (c0 = 1, b = 3; b >= 0; b--) {
        for (i = 0; i < NC; i++) s[i] = (MASK >> i & 1) ? st(((T[i][(h[i] * 16 + c0) & M[i]] >> 16) + .5) / 65536) : 0;
        // discretized logistic: P(bit=1 | values so far) from CDF at bin edges
        int lo = (c0 << (b + 1)) - 16, md = lo + (1 << b), hi = md + (1 << b);
#ifdef GAUSS
#define F(v) ((v) <= 0 ? 0 : (v) >= 16 ? 1 : .5 * erfc(-((v) - .5 - mu) / (sc * 1.41421)))
#else
#define F(v) ((v) <= 0 ? 0 : (v) >= 16 ? 1 : sq(((v) - .5 - mu) / sc))
#endif
        double pl = F(lo), pm = F(md), ph = F(hi), pp = (ph - pm) / (ph - pl + 1e-30);
        pp = pp < 1e-4 ? 1e-4 : pp > 1 - 1e-4 ? 1 - 1e-4 : pp;
        s[NC] = (MASK >> NC & 1) ? st(pp) : 0;
        s[NC + 1] = .3;
        #ifndef WSEL
#define WSEL c0 * 16 + (q > 15 ? 15 : q)
#endif
        float t = 0, *wt = w[WSEL];
        for (i = 0; i < NC + 2; i++) t += wt[i] * s[i];
        float p = sq(t);
        float z = (t < -7.99 ? -7.99 : t > 7.99 ? 7.99 : t) * 2 + 16;
        int a = (c0 * 64 + (mq > 63 ? 63 : mq)) * 33 + (int)z;
        float f = z - (int)z;
        float pa = A[a] * (1 - f) + A[a + 1] * f;
        float pf = (p + 3 * pa) / 4;
        U pi = pf * 65536;
        pi = pi < 32 ? 32 : pi > 65504 ? 65504 : pi;
        y = code(v >> b & 1, pi);
#ifdef DBG
        for (i = 0; i <= NC; i++) { double q1 = sq(s[i]); cost[i] -= log2(y ? q1 : 1 - q1); }
        cost[NC+1] -= log2(y ? p : 1 - p); cost[NC+2] -= log2(y ? pf : 1 - pf); PC[k] -= log2(y ? pf : 1 - pf);
#endif
        float e = y - p;
        for (i = 0; i < NC + 2; i++) wt[i] += MLR * e * s[i];
        for (i = 0; i < NC; i++) {
          U *e2 = &T[i][(h[i] * 16 + c0) & M[i]], m = *e2 & 1023;
          if (m < 255) ++*e2;
          *e2 += ((int)((y << 22) - (*e2 >> 10)) >> 3) * dt[m] & ~1023u;
        }
        A[a] += (y - A[a]) * (1 - f) / 32;
        A[a + 1] += (y - A[a + 1]) * f / 32;
        c0 = c0 * 2 + y;
      }
      I[r][c] = c0 -= 16;
      double er = c0 - mu;
      E[r][c] = fabs(er);
      sv[q] += er * er, sn[q]++;
      if (sn[q] > 500) sv[q] *= .5, sn[q] *= .5;
      ols(c0);
    }
#ifdef DBG
  for (i = 0; i < NC + 3; i++) fprintf(stderr, "in%-2d %8.0f bytes\n", i, cost[i] / 8);
#endif
  if (d) {
    for (i = 0; i < 125000; i++) putchar(I[2 + i / 250][2 + i % 250 * 2] << 4 | I[2 + i / 250][3 + i % 250 * 2]);
  } else putchar((x1 >> 24) + 1);
  return 0;
}
