// Research version of the anchor coder (readable; knobs are -D defines). Not a submission.
//
// Order of coding:
//   1. anchors: the lattice of step T (x % T == 0 and y % T == 0), raster order over the
//      lattice, either verbatim (VERB=1, 4 bits each) or modelled from earlier anchors;
//   2. for t = T, T/2, ..., 2 (h = t/2), knowing the lattice of step t:
//      pass D codes the cell centres (x % t == h, y % t == h), whose 4 diagonal neighbours
//      at distance h are known; pass A then codes the edge midpoints (x % t == h, y % t == 0)
//      and (x % t == 0, y % t == h), whose 4 axial neighbours at distance h are known.
//      After pass A the lattice of step h is known.
//   T = 1 is plain raster coding (every pixel is an "anchor" coded from causal neighbours).
// Every pixel is coded MSB-first as 4 binary decisions, mixing (logistic mixer) hashed
// adaptive counters with a discretized Gaussian centred on an online least-squares (RLS)
// prediction from a stencil of known pixels. Stencils, RLS weights, variance buckets and
// (optionally) mixer weights are per group: anchors, then (level, pass type).
//
// usage: anchor IMAGE [ARCHIVE]   prints bits per group; with ARCHIVE writes the archive.
// -DVERIFY: decode the archive again and compare with the image.
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef T
#define T 16  // anchor lattice step (power of 2)
#endif
#ifndef VERB
#define VERB 0  // 1: anchors verbatim
#endif
#ifndef RA
#define RA 2.9  // anchor stencil radius, in lattice units (causal offsets)
#endif
#ifndef RD
#define RD 3.2  // centre-pass stencil radius, in units of h
#endif
#ifndef RX
#define RX 3.2  // axis-pass stencil radius, in units of h
#endif
#ifndef LAM
#define LAM .9997  // RLS forgetting
#endif
#ifndef NU
#define NU 100  // RLS regulariser
#endif
#ifndef RES
#define RES 0  // 1: RLS on residuals from the mean of the 4 nearest stencil pixels
#endif
#ifndef LR
#define LR .001  // mixer learning rate
#endif
#ifndef LIM
#define LIM 80  // counter count limit
#endif
#ifndef MIXG
#define MIXG 0  // 1: one mixer weight set per group
#endif
#ifndef KQ
#define KQ 3.0  // variance bucket: KQ x (sum of |errors| at the 4 nearest)
#endif
#ifndef KS
#define KS 0.0  // ... + KS x (max - min of the 4 nearest)
#endif
#ifndef VR
#define VR 500  // variance learning rate denominator
#endif
#ifndef NC
#define NC 5  // counter contexts used (of those defined in ctx())
#endif

#ifndef SHR
#define SHR 0  // 1: RLS weights shared across levels (one set per pass type)
#endif
#ifndef SHV
#define SHV 0  // 1: variance buckets shared across levels
#endif
#ifndef SHC
#define SHC 0  // 1: counter contexts shared across levels (level not hashed)
#endif

#ifndef MIXS
#define MIXS 0  // mixer weight set: 0 one, 1 per pass type, 2 per activity bucket q/4, 3 type x (q>7)
#endif
#ifndef DIR
#define DIR 0  // 1: two extra Gaussian inputs centred on the means of opposite pairs of the 4 nearest
#endif
#ifndef DSC
#define DSC 1.0  // their sigma = DSC x the main sigma
#endif
#ifndef RMODE
#define RMODE 0  // 1: after the anchors, code the rest in raster order (stencil reaches ahead)
#endif
#ifndef RL
#define RL 2.3  // RMODE: radius of the look-ahead part of the stencil (reads the fill estimate)
#endif
#ifndef FILL
#define FILL 0  // fill estimate of unknown pixels: 0 lattice point below-left, 1 bilinear
#endif
#ifndef ROWS
#define ROWS 0  // 1: anchors are whole rows (y % T == 0); refinement halves the row spacing
#endif
#define NI (NC + 1 + 2 * DIR)  // mixer inputs

#define W 500
#define NPIX (W * W)
#define NG 40
#define MAXL 140
#define HB 22

typedef struct { int n, dx[MAXL], dy[MAXL]; } Sten;
static Sten SA, SD, SX, SR, SW0, SW1;

static unsigned char X[NPIX];  // the image
static double F[NPIX];  // known value, or fill estimate of an unknown pixel
static int K[NPIX];     // known flag
static double E[NPIX];         // |error| of each coded pixel
static float COST[NPIX];       // bits spent on each pixel (for analysis)
static double P[1 << HB];
static unsigned char N[1 << HB];
static double C[NG][MAXL][MAXL], Bv[NG][MAXL], A[NG][MAXL], V[NG][32];
static double Wm[NG][NI];
static double gbits[NG];
static long gcount[NG];

// arithmetic coder (carry-less, 32-bit range), as in the reference compressor
static unsigned char *buf;
static long bp, blen;
static unsigned lo, hi, xw;
static int dec;

static int bitcode(int y, double p) {  // p = P(y = 1)
  p = p * .9999 + .00005;
  unsigned M = lo + (unsigned)((hi - lo) * p);
  if (dec) y = xw <= M;
  y ? (hi = M) : (lo = M + 1);
  while ((lo ^ hi) < 1u << 24) {
    if (dec) xw = xw << 8 | (bp < blen ? buf[bp] : 0), bp++;
    else buf[bp++] = hi >> 24;
    lo <<= 8, hi = hi << 8 | 255;
  }
  return y;
}

static int cmpd(const void *a, const void *b) {
  const int *p = a, *q = b;
  int d = p[0] - q[0];
  return d ? d : p[1] != q[1] ? p[1] - q[1] : p[2] - q[2];
}

// kind 0: anchor raster (causal offsets); 1: centres (known: odd/odd, or even/even causal);
// 3: raster after anchors (causal within r, plus look-ahead within RL);
// 4: anchor rows (dx in pixels, dy in units of T; causal); 5: refined rows (dy odd = known
// rows above/below, dy even = rows of the same pass, known if causal);
// 2: axis pixels (known: dx+dy odd, or dx+dy even causal). Sorted by distance.
static void mksten(Sten *s, int kind, double r) {
  int t[4 * MAXL * 4][3], n = 0, R = (int)r;
  for (int dy = -R; dy <= R; dy++)
    for (int dx = -R; dx <= R; dx++) {
      int d2 = dx * dx + dy * dy, causal = dy < 0 || (dy == 0 && dx < 0), ok;
      if (!d2 || d2 > r * r) continue;
      if (kind == 0) ok = causal;
      else if (kind == 1) ok = (dx & 1 && dy & 1) || (!(dx & 1) && !(dy & 1) && causal);
      else if (kind == 2) ok = ((dx + dy) & 1) || causal;
      else if (kind == 3) ok = causal || d2 <= RL * RL;
      else if (kind == 4) ok = causal;
      else ok = (dy & 1) || causal;
      // kind 3: causal offsets first, so the 4 nearest used by contexts are real pixels
      if (ok) t[n][0] = d2 + (kind == 3 && !causal) * 10000, t[n][1] = dy, t[n][2] = dx, n++;
    }
  qsort(t, n, sizeof t[0], cmpd);
  if (n > MAXL) n = MAXL;
  s->n = n;
  for (int i = 0; i < n; i++) s->dy[i] = t[i][1], s->dx[i] = t[i][2];
}

static double st(double p) { p = p * .9998 + .0001; return log(p / (1 - p)); }
static double cdf(int v, double m, double sg) {  // P(value < v) under N(m, sg^2), tails folded
  return v <= 0 ? 0 : v >= 16 ? 1 : .5 * erfc((m + .5 - v) / (sg * M_SQRT2));
}
static int clampi(int v) { return v < 0 ? 0 : v > W - 1 ? W - 1 : v; }

static int shg(int g, int on) { return on && g ? 1 + (g - 1) % 3 : g; }

static void code_pixel(int x, int y, int g, Sten *s, int ux, int uy) {
  int L = s->n, k = y * W + x, i, j, gr = shg(g, SHR), gv = shg(g, SHV), gc = shg(g, SHC);
  double o[MAXL], m = 0, base = 0;
  int nb[MAXL];
  for (i = 0; i < L; i++) {
    int kk = clampi(y + s->dy[i] * uy) * W + clampi(x + s->dx[i] * ux);
    nb[i] = kk, o[i] = F[kk];
  }
  if (RES) base = (o[0] + o[1] + o[2] + o[3]) / 4;
  for (i = 0; i < L; i++) o[i] -= base, m += A[gr][i] * o[i];
  m += base;
  int mn = 99, mx = -1;
  double es = 0;
  for (i = 0; i < 4; i++) {
    int v = (int)(F[nb[i]] + .5);
    mn = v < mn ? v : mn, mx = v > mx ? v : mx, es += E[nb[i]];
  }
  int q = (int)fmin(KQ * es + KS * (mx - mn), 31);
  double sg = sqrt(V[gv][q]);
#define Q(i) (int)(F[nb[i]] + .5)
  int a0 = Q(0), a1 = Q(1), a2 = Q(2), a3 = Q(3), a4 = Q(4), a5 = Q(5);
  int rm = (int)fmin(fmax(m * 2 + .5, 0), 31);
  unsigned h[8];
  h[0] = a0 << 12 | a1 << 8 | a2 << 4 | a3;
  h[1] = a0 << 4 | a1;
  h[2] = a2 << 4 | a3;
  h[3] = rm << 5 | q;
  h[4] = a0 << 12 | a1 << 8 | a4 << 4 | a5;
  h[5] = a0 << 8 | a2 << 4 | a4;
  double *w = Wm[MIXG ? g : MIXS == 1 ? shg(g, 1) : MIXS == 2 ? q / 4 : MIXS == 3 ? shg(g, 1) * 2 + (q > 7) : 0], sx[NI];
  double md[2] = {(a0 + a3) / 2., (a1 + a2) / 2.};
  int v = X[k], c = 1, lov = 0;
  unsigned H[NC];
  for (int e = 3; e >= 0; e--) {  // remaining values: lov .. lov + 2^(e+1) - 1
    int mid = lov + (1 << e), top = lov + (2 << e);
    double t = cdf(top, m, sg), pg = (t - cdf(mid, m, sg)) / (t - cdf(lov, m, sg) + 1e-30);
    sx[NC] = st(pg);
    double f = w[NC] * sx[NC];
    for (i = 0; i < 2 * DIR; i++) {
      t = cdf(top, md[i], sg * DSC);
      sx[NC + 1 + i] = st((t - cdf(mid, md[i], sg * DSC)) / (t - cdf(lov, md[i], sg * DSC) + 1e-30));
      f += w[NC + 1 + i] * sx[NC + 1 + i];
    }
    for (i = 0; i < NC; i++) {
      H[i] = ((h[i] << 4 | c) * 2654435761u ^ (unsigned)(gc * 0x9E3779B1u) ^ (unsigned)i * 0x85EBCA77u) * 2246822519u >> (32 - HB);
      sx[i] = st(P[H[i]]), f += w[i] * sx[i];
    }
    double p = 1 / (1 + exp(-f));
    int y1 = bitcode(v >> e & 1, p);
    double pp = p * .9999 + .00005;
    gbits[g] -= log2(y1 ? pp : 1 - pp), COST[k] -= log2(y1 ? pp : 1 - pp);
    for (i = 0; i < NI; i++) w[i] += LR * (y1 - p) * sx[i];
    for (i = 0; i < NC; i++) j = H[i], N[j] += N[j] < LIM, P[j] += (y1 - P[j]) / (N[j] + .5);
    c = c * 2 + y1, lov += y1 << e;
  }
  v = lov;
  X[k] = v;  // in decode mode this writes the decoded pixel
  gcount[g]++;
  double err = v - m;
  E[k] = fabs(err);
  V[gv][q] += (err * err - V[gv][q]) / VR;
  double yv = v - base;
  for (i = 0; i < L; i++) {
    double r = Bv[gr][i] = Bv[gr][i] * LAM + o[i] * yv;
    for (j = 0; j < L; j++) C[gr][i][j] = C[gr][i][j] * LAM + o[i] * o[j];
    for (j = 0; j < L; j++) if (j != i) r -= C[gr][i][j] * A[gr][j];
    A[gr][i] = r / (C[gr][i][i] + NU);
  }
  F[k] = v, K[k] = 1;
}

static void code_verbatim(int x, int y) {
  int k = y * W + x, v = 0;
  for (int e = 3; e >= 0; e--) v = v * 2 + bitcode(X[k] >> e & 1, .5);
  X[k] = v, F[k] = v, K[k] = 1, E[k] = 0;
  gbits[0] += 4, gcount[0]++;
}

static void fill(int u) {  // fill estimate of unknown pixels from the known lattice of step u
  int last = (W - 1) / u * u;
  for (int y = 0; y < W; y++)
    for (int x = 0; x < W; x++) {
      if (K[y * W + x]) continue;
      int x0 = x - x % u, y0 = y - y % u, x1 = x0 + u > last ? x0 : x0 + u, y1 = y0 + u > last ? y0 : y0 + u;
      double fx = x1 > x0 ? (double)(x - x0) / u : 0, fy = y1 > y0 ? (double)(y - y0) / u : 0;
      F[y * W + x] = FILL ? (1 - fy) * ((1 - fx) * F[y0 * W + x0] + fx * F[y0 * W + x1]) +
                                fy * ((1 - fx) * F[y1 * W + x0] + fx * F[y1 * W + x1])
                          : F[y0 * W + x0];
    }
}

static void fill_rows(int u) {  // unknown rows: linear interpolation between known rows
  int last = (W - 1) / u * u;
  for (int y = 0; y < W; y++)
    for (int x = 0; x < W; x++) {
      if (K[y * W + x]) continue;
      int y0 = y - y % u, y1 = y0 + u > last ? y0 : y0 + u;
      double fy = y1 > y0 ? (double)(y - y0) / u : 0;
      F[y * W + x] = FILL ? (1 - fy) * F[y0 * W + x] + fy * F[y1 * W + x] : F[y0 * W + x];
    }
}

static void reset(void) {
  memset(F, 0, sizeof F), memset(K, 0, sizeof K), memset(E, 0, sizeof E);
  memset(N, 0, sizeof N), memset(C, 0, sizeof C), memset(Bv, 0, sizeof Bv), memset(A, 0, sizeof A);
  memset(Wm, 0, sizeof Wm), memset(gbits, 0, sizeof gbits), memset(gcount, 0, sizeof gcount);
  for (int i = 0; i < 1 << HB; i++) P[i] = .5;
  for (int g = 0; g < NG; g++) for (int q = 0; q < 32; q++) V[g][q] = 1;
  for (int i = 0; i < NPIX; i++) F[i] = 8;
  lo = 0, hi = -1, bp = 0;
}

static void run(void) {
  if (ROWS) {
    for (int y = 0; y < W; y += T)
      for (int x = 0; x < W; x++) code_pixel(x, y, 0, &SW0, 1, T);
    fill_rows(T);
    for (int t = T, lv = 0; t >= 2; t /= 2, lv++) {
      int h = t / 2;
      for (int y = h; y < W; y += t)
        for (int x = 0; x < W; x++) code_pixel(x, y, 1 + 3 * lv, &SW1, 1, h);
      fill_rows(h);
    }
    if (!dec) buf[bp++] = (lo >> 24) + 1;
    return;
  }
  for (int y = 0; y < W; y += T)
    for (int x = 0; x < W; x += T)
      VERB && T > 1 ? code_verbatim(x, y) : code_pixel(x, y, 0, &SA, T, T);
  fill(T);
  if (RMODE) {
    // probe only (not transmitted): extra anchors from a mask file, given to the model for free
    static unsigned char mask[NPIX];
    static int nmask = -1;
    if (nmask < 0) {
      FILE *mf = getenv("ANCHMASK") ? fopen(getenv("ANCHMASK"), "rb") : 0;
      nmask = 0;
      if (mf && fread(mask, 1, NPIX, mf) != NPIX) nmask = 0;
      if (mf) fclose(mf);
      for (int i = 0; i < NPIX; i++) nmask += mask[i] && !K[i];
      if (nmask) fprintf(stderr, "free mask anchors: %d\n", nmask);
    }
    for (int i = 0; i < NPIX; i++)
      if (mask[i] && !K[i]) F[i] = X[i], K[i] = 1;
    for (int y = 0; y < W; y++)
      for (int x = 0; x < W; x++)
        if (!K[y * W + x]) code_pixel(x, y, 1, &SR, 1, 1);
    if (!dec) buf[bp++] = (lo >> 24) + 1;
    return;
  }
  for (int t = T, lv = 0; t >= 2; t /= 2, lv++) {
    int h = t / 2, g = 1 + 3 * lv;
    for (int y = h; y < W; y += t)
      for (int x = h; x < W; x += t) code_pixel(x, y, g, &SD, h, h);
    for (int y = 0; y < W; y += h)
      for (int x = y % t ? 0 : h; x < W; x += t) code_pixel(x, y, g + 1 + (y % t != 0), &SX, h, h);
    fill(h);
  }
  if (!dec) buf[bp++] = (lo >> 24) + 1;
}

int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "rb");
  unsigned char raw[NPIX / 2];
  if (!f || fread(raw, 1, NPIX / 2, f) != NPIX / 2) return fprintf(stderr, "bad image\n"), 1;
  fclose(f);
  for (int i = 0; i < NPIX; i++) X[i] = raw[i / 2] >> (i % 2 ? 0 : 4) & 15;
  mksten(&SA, 0, RA), mksten(&SD, 1, RD), mksten(&SX, 2, RX), mksten(&SR, 3, RA > RL ? RA : RL), mksten(&SW0, 4, RA), mksten(&SW1, 5, RD);
  buf = malloc(1 << 20);
  reset(), dec = 0, run();
  long alen = bp;
  double tot = 0;
  printf("stencils: anchor %d, centre %d, axis %d, raster-after-anchors %d\n", SA.n, SD.n, SX.n, SR.n);
  for (int g = 0; g < NG; g++)
    if (gcount[g]) {
      tot += gbits[g];
      static const char *kind[3] = {"centre", "axis-h", "axis-v"};
      char lab[32];
      if (g && RMODE) snprintf(lab, sizeof lab, "rest (raster)");
      else if (g && ROWS) snprintf(lab, sizeof lab, "rows t=%d", T >> (g - 1) / 3);
      else if (g) snprintf(lab, sizeof lab, "t=%d %s", T >> (g - 1) / 3, kind[(g - 1) % 3]);
      else snprintf(lab, sizeof lab, VERB && T > 1 ? "anchors verbatim" : "anchors");
      printf("  %-16s %7ld px %9.0f bits %6.3f bpp\n", lab, gcount[g], gbits[g], gbits[g] / gcount[g]);
    }
  printf("archive %ld bytes (ideal %.0f), %.4f bpp\n", alen, tot / 8, alen * 8.0 / NPIX);
  if (argc > 2) { FILE *o = fopen(argv[2], "wb"); fwrite(buf, 1, alen, o); fclose(o); }
  if (getenv("COSTMAP")) { FILE *o = fopen(getenv("COSTMAP"), "wb"); fwrite(COST, 4, NPIX, o); fclose(o); }
#ifdef VERIFY
  unsigned char orig[NPIX];
  memcpy(orig, X, NPIX);
  memset(X, 0, NPIX);
  reset(), dec = 1, blen = alen;
  for (int i = 0; i < 4; i++) xw = xw << 8 | buf[bp++];
  run();
  printf("verify: %s\n", memcmp(orig, X, NPIX) ? "MISMATCH" : "ok");
#endif
  return 0;
}
