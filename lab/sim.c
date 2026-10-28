// C port of gen_image.py (the source simulator, readable lab version).
// Same recipe and parameter distributions; its own PRNG (splitmix64 + Box-Muller), and the
// Gaussian random fields are made by an inverse FFT on a 512x512 grid cropped to 500x500
// (numpy uses a 500x500 grid; same spectrum shape).
//
// usage: sim SEED N OUTPREFIX   writes N images OUTPREFIX000.bin ... in the image.bin format
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define W 500
#define G 512
typedef double complex Q;
static unsigned long long rs;
static double img[W * W], fld[W * W];
static Q buf[G * G];

static double unif(void) {  // splitmix64 -> [0, 1)
  unsigned long long z = rs += 0x9E3779B97F4A7C15ull;
  z = (z ^ z >> 30) * 0xBF58476D1CE4E5B9ull, z = (z ^ z >> 27) * 0x94D049BB133111EBull;
  return ((z ^ z >> 31) >> 11) * 0x1p-53;
}
static double gauss(void) { return sqrt(-2 * log(1 - unif())) * cos(2 * M_PI * unif()); }

static void fft(Q *a, int s) {  // in-place inverse-direction FFT of G points with stride s
  for (int i = 1, j = 0; i < G; i++) {
    int k = G >> 1;
    for (; j & k; k >>= 1) j ^= k;
    j ^= k;
    if (i < j) { Q t = a[i * s]; a[i * s] = a[j * s], a[j * s] = t; }
  }
  for (int len = 2; len <= G; len <<= 1) {
    Q wl = cexp(2 * M_PI * I / len);
    for (int i = 0; i < G; i += len) {
      Q w = 1;
      for (int j = 0; j < len / 2; j++, w *= wl) {
        Q u = a[(i + j) * s], v = a[(i + j + len / 2) * s] * w;
        a[(i + j) * s] = u + v, a[(i + j + len / 2) * s] = u - v;
      }
    }
  }
}

// unit-variance Gaussian random field with power spectrum ~ 1/f^beta, into fld
static void grf(double beta) {
  for (int y = 0; y < G; y++)
    for (int x = 0; x < G; x++) {
      double fy = (y < G / 2 ? y : y - G) / (double)G, fx = (x < G / 2 ? x : x - G) / (double)G;
      double f = hypot(fx, fy);
      buf[y * G + x] = f > 0 ? (gauss() + I * gauss()) * pow(f, -beta / 2) : 0;
    }
  for (int y = 0; y < G; y++) fft(buf + y * G, 1);
  for (int x = 0; x < G; x++) fft(buf + x, G);
  double s = 0, s2 = 0;
  for (int y = 0; y < W; y++)
    for (int x = 0; x < W; x++) {
      double v = creal(buf[y * G + x]);
      fld[y * W + x] = v, s += v, s2 += v * v;
    }
  double mu = s / (W * W), sd = sqrt(s2 / (W * W) - mu * mu);
  for (int i = 0; i < W * W; i++) fld[i] = (fld[i] - mu) / sd;
}

static double soft(double d) { return d < -.5 ? 1 : d > .5 ? 0 : .5 - d; }
static int cmpd(const void *a, const void *b) { double d = *(double *)a - *(double *)b; return (d > 0) - (d < 0); }

// one image -> pix (values 0..15, row-major)
void simulate(unsigned char *pix) {
  static double tmp[W * W];
  grf(3.4); for (int i = 0; i < W * W; i++) img[i] = fld[i];
  grf(2.2); for (int i = 0; i < W * W; i++) img[i] += .3 * fld[i];
  grf(1.0); for (int i = 0; i < W * W; i++) img[i] += .08 * fld[i];
  double gx = gauss() * .6, gy = gauss() * .6;
  for (int y = 0; y < W; y++)
    for (int x = 0; x < W; x++) img[y * W + x] += gx * (x - W / 2) / W + gy * (y - W / 2) / W;
  int nobj = 30 + (int)(unif() * 15);
  for (int o = 0; o < nobj; o++) {
    double pk = unif(), cx = unif() * W, cy = unif() * W, r = 8 + unif() * 52, th = unif() * M_PI;
    // kinds: disc .22, ellipse .18, rect .22, ring .10, line .14, texture .14
    int kind = pk < .22 ? 0 : pk < .40 ? 1 : pk < .62 ? 2 : pk < .72 ? 3 : pk < .86 ? 4 : 5;
    double c = cos(th), s = sin(th), a = r, b = r, w = 0, len = 0, period = 0;
    if (kind == 1) b = r * (.3 + unif() * .6);
    if (kind == 2) b = r * (.3 + unif() * .7);
    if (kind == 3) w = 2 + unif() * 4;
    if (kind == 4) w = .8 + unif() * 2.2, len = 60 + unif() * 240;
    if (kind == 5) period = 5 + unif() * 13;
    double amp = gauss() * 1.1;
    int opaque = kind != 5 && unif() < .5;
    double level = 0;
    if (opaque) level = gauss() * 1.2, grf(3.0);
    for (int y = 0; y < W; y++)
      for (int x = 0; x < W; x++) {
        double u = (x - cx) * c + (y - cy) * s, v = -(x - cx) * s + (y - cy) * c, m;
        if (kind == 0) m = soft(hypot(u, v) - r);
        else if (kind == 1) m = soft((hypot(u / a, v / b) - 1) * fmin(a, b));
        else if (kind == 2) m = soft(fmax(fabs(u) - a, fabs(v) - b));
        else if (kind == 3) m = soft(fabs(hypot(u, v) - r) - w);
        else if (kind == 4) m = soft(fmax(fabs(v) - w, fabs(u) - len));
        else m = soft(hypot(u, v) - r) * sin(2 * M_PI * u / period);
        int i = y * W + x;
        img[i] = opaque ? img[i] * (1 - m) + (level + .25 * fld[i]) * m : img[i] + amp * m;
      }
  }
  memcpy(tmp, img, sizeof img);
  qsort(tmp, W * W, sizeof(double), cmpd);
  // numpy percentile (linear interpolation) at 0.5 and 99.5
  double plo = .005 * (W * W - 1), phi = .995 * (W * W - 1);
  int il = (int)plo, ih = (int)phi;
  double lo = tmp[il] + (plo - il) * (tmp[il + 1] - tmp[il]), hi = tmp[ih] + (phi - ih) * (tmp[ih + 1] - tmp[ih]);
  for (int i = 0; i < W * W; i++) {
    double q = (img[i] - lo) / (hi - lo) * 16 + gauss() * .28;
    pix[i] = q < 0 ? 0 : q >= 15 ? 15 : (int)q;
  }
}

#ifndef NOMAIN
int main(int argc, char **argv) {
  if (argc < 4) return fprintf(stderr, "usage: sim SEED N OUTPREFIX\n"), 1;
  rs = strtoull(argv[1], 0, 10);
  int n = atoi(argv[2]);
  static unsigned char pix[W * W], packed[W * W / 2];
  for (int k = 0; k < n; k++) {
    simulate(pix);
    for (int i = 0; i < W * W / 2; i++) packed[i] = pix[2 * i] << 4 | pix[2 * i + 1];
    char fn[512];
    snprintf(fn, sizeof fn, "%s%03d.bin", argv[3], k);
    FILE *f = fopen(fn, "wb");
    fwrite(packed, 1, sizeof packed, f), fclose(f);
  }
  return 0;
}
#endif
