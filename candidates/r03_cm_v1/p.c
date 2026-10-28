// v1: binary-decomposition context mixing. Same program compresses (stdin = image.bin,
// 125000 bytes) and decompresses (stdin = archive). Pixels coded MSB-first as 4 binary
// decisions; each decision mixes neighbour-context counters with a logistic mixer + APM.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned U;
typedef unsigned long long Q;
#define NC 8
unsigned char B[1 << 17], I[504][504];
U *T[NC], M[NC], h[NC], x1, x2 = -1, X, d, ip, dt[1024];
float w[16][NC + 1], s[NC + 1];
float A[16 * 16 * 33];  // APM: [node][N][bucket] -> P(1)

// binary arithmetic coder, encodes or decodes bit y with P(1) = p/65536
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
float st(float p) { return log(p / (1 - p)); }
float sq(float x) { return 1 / (1 + exp(-x)); }

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
  for (i = 0; i < 16 * 16 * 33; i++) A[i] = sq((i % 33 - 16) / 2.0);
  for (i = 0; i < 16; i++) for (j = 0; j <= NC; j++) w[i][j] = .3;
  if (d) for (i = 0; i < 4; i++) X = X << 8 | B[ip++];
  for (int r = 2; r < 502; r++)
    for (int c = 2; c < 502; c++) {
      U W = I[r][c - 1], N = I[r - 1][c], NW = I[r - 1][c - 1], NE = I[r - 1][c + 1],
        WW = I[r][c - 2], NN = I[r - 2][c], NNE = I[r - 2][c + 1];
      int k = (r - 2) * 500 + c - 2, v = d ? 0 : (B[k / 2] >> (k & 1 ? 0 : 4)) & 15;
      h[0] = 0;
      h[1] = W << 4 | N;
      h[2] = W << 8 | N << 4 | NW;
      h[3] = W << 8 | N << 4 | NE;
      h[4] = W << 12 | N << 8 | NW << 4 | NE;
      h[5] = N << 8 | NN << 4 | NE;
      h[6] = W << 8 | WW << 4 | NW;
      h[7] = (W << 20 | N << 16 | NW << 12 | NE << 8 | WW << 4 | NN) * 2654435761u >> 14;
      for (c0 = 1, b = 3; b >= 0; b--) {
        for (i = 0; i < NC; i++) s[i] = st(((T[i][(h[i] * 16 + c0) & M[i]] >> 16) + .5) / 65536);
        s[NC] = .3;
        float t = 0;
        for (i = 0; i <= NC; i++) t += w[c0][i] * s[i];
        float p = sq(t);
        // APM on (node, N) interpolating between 33 buckets of stretch(p)
        float z = (t < -7.99 ? -7.99 : t > 7.99 ? 7.99 : t) * 2 + 16;
        int a = (c0 * 16 + N) * 33 + (int)z;
        float f = z - (int)z;
        float pa = A[a] * (1 - f) + A[a + 1] * f;
        float pf = (p + 3 * pa) / 4;
        U pi = pf * 65536;
        pi = pi < 32 ? 32 : pi > 65504 ? 65504 : pi;
        y = code(v >> b & 1, pi);
        float e = y - p;
        for (i = 0; i <= NC; i++) w[c0][i] += .02 * e * s[i];
        for (i = 0; i < NC; i++) {
          U *e2 = &T[i][(h[i] * 16 + c0) & M[i]], m = *e2 & 1023;
          if (m < 255) ++*e2;
          *e2 += ((int)((y << 22) - (*e2 >> 10)) >> 3) * dt[m] & ~1023u;
        }
        A[a] += (y - A[a]) * (1 - f) / 32;
        A[a + 1] += (y - A[a + 1]) * f / 32;
        c0 = c0 * 2 + y;
      }
      I[r][c] = c0 - 16;
    }
  if (d) {
    for (i = 0; i < 125000; i++) putchar(I[2 + i / 250][2 + i % 250 * 2] << 4 | I[2 + i / 250][3 + i % 250 * 2]);
  } else putchar((x1 >> 24) + 1);
  return 0;
}
