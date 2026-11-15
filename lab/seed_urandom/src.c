// Variant of candidates/016_seed_replay whose 32-byte archive is what gen_image.py read from
// os.urandom, not the generator state derived from it (lab probe: it measures what numpy's
// seed hashing costs in program bytes). The hash is numpy's SeedSequence for a 256-bit
// integer; a seed whose top 32 bits are zero would hash differently and is not handled.
// Everything else is 016.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// Z: the seed, then PCG64 state and increment; M: its multiplier
__uint128_t Z[2], M = (__uint128_t)2549297995355413924 << 64 | 4865540595714422341;
// E: Z as 32-bit words; p: SeedSequence pool; h, m: its hash constant and multiplier
unsigned *E = (void *)Z, p[4], h = 0x43b0d7e5, m = 0x931e8875, v;
long K[257];  // ziggurat: accept at once if the 52-bit draw is below K
int i, j, k, x, y, o, l, O;
// w, f: ziggurat strip widths / 2^52 and heights; I: the image before quantizing; G: one
// random field, later the sorted image; R, V: ziggurat tail start and strip area
double w[256], f[256], I[250000], G[250000], R = 3.654152885361009, V = .004928673233974652, T = 0x1p52, d, t, a, b, r, e, X, Y, C, S;
// A: half spectrum (500 rows x 251 columns); B: after the pass over rows; W: 500th roots of 1
_Complex double A[125500], B[125500], W[500];

unsigned H(unsigned v) { return v ^= h, h *= m, v *= h, v ^ v >> 16; }
// PCG64 (XSL-RR 128/64)
size_t Q() {
  *Z = *Z * M + Z[1];
  size_t a = *Z >> 64 ^ *Z;
  int r = *Z >> 122;
  return a >> r | a << (-r & 63);
}
double u() { return (Q() >> 11) * 0x1p-53; }
// standard normal, numpy's random_standard_normal
double n() {
  size_t r = Q(), a = r << 3 >> 12;
  int i = r & 255;
  double x = a * w[i] * (r & 256 ? -1 : 1);
  if (a < K[i]) return x;
  if (i) return (f[i - 1] - f[i]) * u() + f[i] < exp(-x * x / 2) ? x : n();
  do t = -log1p(-u()) / R, e = -log1p(-u());  // the tail beyond R; t, e are free here
  while (e + e <= t * t);
  return a & 256 ? -R - t : R + t;
}
// G = e * (random field with power spectrum f^(4b), unit variance): numpy's grf(-4b). While
// O is 0 (the three terrain fields) G is also added to the image.
// Real parts of the half spectrum are drawn first, then imaginary parts. The inverse transform
// is the direct sum: B[y][k] = sum_j A[j][k] W^(jy), then G[y][x] = Re sum_k B[y][k] W^(kx),
// with columns 1..249 counted twice for their mirror images (folded into A).
void g(double b, double e) {
  double q = 0;
  for (i = 0; i < 125500; i++) A[i] = n();
  for (i = 0; i < 125500; i++) x = i % 251, y = (i / 251 + 250) % 500 - 250, A[i] = (A[i] + n() * 1i) * pow(x * x + y * y, b) * (x % 250 ? 2 : 1);
  *A = 0;
  for (i = 0; i < 125500; i++) {
    B[i] = 0;
    for (j = 0; j < 500; j++) B[i] += A[j * 251 + i % 251] * W[j * (i / 251) % 500];
  }
  for (i = 0; i < 250000; i++) {
    d = 0;
    for (k = 0; k < 251; k++) d += B[i / 500 * 251 + k] * W[k * i % 500];
    q += d * d;
    G[i] = d;
  }
  for (i = 0; i < 250000; i++) I[i] += !O * (G[i] *= e / sqrt(q / 250000));
}
int P(double *a, double *b) { return *a > *b ? 1 : -1; }
int main() {
  // numpy SeedSequence: 8 seed words -> pool of 4 -> 8 state words; then PCG64's own seeding
  fread(Z, 16, 2, stdin);
  for (i = 0; i < 4; i++) p[i] = H(E[i]);
  for (i = 0; i < 8; i++)
    for (j = 0; j < 4; j++)
      if (i - j) v = 0xca01f9dd * p[j] - 0x4973f715 * H(i < 4 ? p[i] : E[i]), p[j] = v ^ v >> 16;
  h = 0x8b51f9dd;
  m = 0x58f38ded;
  for (i = 0; i < 8; i++) E[i ^ 2] = H(p[i % 4]);
  Z[1] = Z[1] * 2 + 1;
  *Z = (Z[1] + *Z) * M + Z[1];
  // ziggurat tables, from the widest strip (x = R) up
  t = d = R;
  for (i = 255; i; i--) w[i] = d / T, f[i] = exp(-d * d / 2), K[i + 1] = d / t * T, t = d, d = sqrt(-2 * log(V / d + f[i]));
  *K = R * f[255] / V * T;
  *w = V / f[255] / T;
  *f = 1;
  for (i = 0; i < 500; i++) t = i * M_PI / 250, W[i] = cos(t) + 1i * sin(t);
  // terrain: grf(3.4) + .3 grf(2.2) + .08 grf(1), plus a random linear gradient
  g(-.85, 1);
  g(-.55, .3);
  g(-.25, .08);
  a = .6 * n();
  b = .6 * n();
  for (i = 0; i < 250000; i++) I[i] += (a * (i % 500 - 250) + b * (i / 500 - 250)) / 500;
  // 30..44 objects; l: 0 disc, 1 ellipse, 2 rectangle, 3 ring, 4 line, 5 striped disc
  for (o = 30 + ((unsigned)Q() * 15L >> 32); o--;) {
    t = u();
    l = (t >= .22) + (t >= .4) + (t >= .62) + (t >= .72) + (t >= .86);
    X = 500 * u();
    Y = 500 * u();
    r = 8 + 52 * u();
    t = M_PI * u();
    C = cos(t);
    S = sin(t);
    // b: second half-axis, ring or line half-width, stripe period; a line is a rectangle
    // whose long half-axis replaces r
    l > 4 ? b = 5 + 13 * u() : l > 3 ? b = .8 + 2.2 * u(), r = 60 + 240 * u() : l > 2 ? b = 2 + 4 * u() : l ? b = r * (.3 + (.5 + l * .1) * u()) : 0;
    a = 1.1 * n();
    // O: the object is opaque and replaces the image with its own level a + G
    if (O = l < 5 && u() < .5) a = 1.2 * n(), g(-.75, .25);
    for (i = 0; i < 250000; i++) {
      e = i % 500 - X;
      d = i / 500 - Y;
      t = e * C + d * S;
      d = d * C - e * S;
      e = hypot(t, d);
      // signed distance to the object's edge, then the anti-aliased mask
      d = l == 1 ? (hypot(t / r, d / b) - 1) * b : l == 3 ? fabs(e - r) - b : l % 2 | !l ? e - r : fmax(fabs(t) - r, fabs(d) - b);
      d = fmin(fmax(.5 - d, 0), 1);
      if (l > 4) d *= sin(2 * M_PI * t / b);
      I[i] += d * (a + O * (G[i] - I[i]));
    }
  }
  // percentiles .5 and 99.5 map to levels 0 and 16; sensor noise; two pixels per byte
  for (i = 0; i < 250000; i++) G[i] = I[i];
  qsort(G, 250000, 8, (void *)P);
  a = G[1249] * .005 + G[1250] * .995;
  b = G[248749] * .995 + G[248750] * .005 - a;
  for (i = 0; i < 250000; i++) j = fmin(fmax((I[i] - a) / b * 16 + .28 * n(), 0), 15), i % 2 ? putchar(k * 16 + j) : (k = j);
}
