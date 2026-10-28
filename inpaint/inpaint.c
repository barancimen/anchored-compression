// Inpainting from anchors with the quadtree anchor model.
//
// The model is the quadtree path of lab/anchor.c (the research anchor coder), copied with its
// model unchanged; the other coding orders of the lab coder are left out, and the lines that
// differ are marked INPAINT. The arithmetic coder is removed:
//   - Anchors are the lattice of step T (x % T == 0 and y % T == 0). The rest is filled coarse
//     to fine: for t = T, T/2, ..., 2 the cell centres, then the edge midpoints.
//   - Training: the model first codes each TRAIN image completely (true pixels, learning on),
//     then the anchors of the target (learning on). That is all it learns.
//   - Sampling: after the anchors the process forks once per sample (one at a time). Each
//     child draws every remaining bit from the model's probability with its own seed and
//     learns nothing: counters, mixer weights, RLS weights and error variances keep their
//     values. Only the error buckets follow the drawn pixels.
//   - Reference: one more child runs the same frozen model over the true pixels, and finally
//     the parent codes them with learning on, exactly as the compressor would. With no TRAIN
//     images its total cost equals lab/anchor.c's for the same knobs (port check).
//
// usage: inpaint TARGET NSAMPLES SEED OUTDIR [TRAIN ...]
// Writes OUTDIR/sample_NNN.bin (the image.bin format) and one tab-separated line per sample:
// index, seed, bits per drawn pixel (the model's own cost of what it drew).
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

// knobs: the best lab quadtree settings (lab/notes.md), anchor lattice step 8
#ifndef T
#define T 8  // anchor lattice step (power of 2)
#endif
#ifndef RA
#define RA 5.5  // anchor stencil radius, in lattice units (causal offsets)
#endif
#ifndef RD
#define RD 6.5  // centre-pass stencil radius, in units of h
#endif
#ifndef RX
#define RX 5.9  // axis-pass stencil radius, in units of h
#endif
#ifndef LAM
#define LAM .9997  // RLS forgetting
#endif
#ifndef NU
#define NU 100  // RLS regulariser
#endif
#ifndef RES
#define RES 1  // 1: RLS on residuals from the mean of the 4 nearest stencil pixels
#endif
#ifndef LR
#define LR .001  // mixer learning rate
#endif
#ifndef LIM
#define LIM 80  // counter count limit
#endif
#ifndef KQ
#define KQ 1.5  // variance bucket: KQ x (sum of |errors| at the 4 nearest)
#endif
#ifndef KS
#define KS 2.0  // ... + KS x (max - min of the 4 nearest)
#endif
#ifndef VR
#define VR 500  // variance learning rate denominator
#endif
#ifndef NC
#define NC 6  // counter contexts used
#endif
#ifndef SHR
#define SHR 0  // 1: RLS weights shared across levels (one set per pass type)
#endif
#ifndef SHV
#define SHV 1  // 1: variance buckets shared across levels
#endif
#ifndef SHC
#define SHC 1  // 1: counter contexts shared across levels (level not hashed)
#endif
#ifndef DIR
#define DIR 1  // 1: two extra Gaussian inputs centred on the means of opposite pairs of the 4 nearest
#endif
#ifndef DSC
#define DSC 0.7  // their sigma = DSC x the main sigma
#endif
#ifndef FILL
#define FILL 0  // fill estimate of unknown pixels: 0 lattice point below-left, 1 bilinear
#endif
#define NI (NC + 1 + 2 * DIR)  // mixer inputs

#define W 500
#define NPIX (W * W)
#define NG 40
#define MAXL 140
#define HB 22

typedef struct { int n, dx[MAXL], dy[MAXL]; } Sten;
static Sten SA, SD, SX;

static unsigned char X[NPIX];  // the image (a sample overwrites its own copy)
static double F[NPIX];         // known value, or fill estimate of an unknown pixel
static int K[NPIX];            // known flag
static double E[NPIX];         // |error| of each coded pixel
static double P[1 << HB];
static unsigned char N[1 << HB];
static double C[NG][MAXL][MAXL], Bv[NG][MAXL], A[NG][MAXL], V[NG][32];
static double Wm[NI];
static double gbits[NG];
static long gcount[NG];

// INPAINT: learn = update the model; smp = sample index in a child, -1 otherwise;
// rng = splitmix64 state
static int learn = 1, smp = -1;
static unsigned long long rng;

// INPAINT: uniform double in [0, 1) from splitmix64
static double uni(void) {
  unsigned long long z = rng += 0x9E3779B97F4A7C15ull;
  z = (z ^ z >> 30) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ z >> 27) * 0x94D049BB133111EBull;
  return ((z ^ z >> 31) >> 11) * 0x1p-53;
}

static int cmpd(const void *a, const void *b) {
  const int *p = a, *q = b;
  int d = p[0] - q[0];
  return d ? d : p[1] != q[1] ? p[1] - q[1] : p[2] - q[2];
}

// kind 0: anchor raster (causal offsets); 1: centres (known: odd/odd, or even/even causal);
// 2: axis pixels (known: dx+dy odd, or dx+dy even causal). Sorted by distance.
static void mksten(Sten *s, int kind, double r) {
  int t[4 * MAXL * 4][3], n = 0, R = (int)r;
  for (int dy = -R; dy <= R; dy++)
    for (int dx = -R; dx <= R; dx++) {
      int d2 = dx * dx + dy * dy, causal = dy < 0 || (dy == 0 && dx < 0), ok;
      if (!d2 || d2 > r * r) continue;
      if (kind == 0) ok = causal;
      else if (kind == 1) ok = (dx & 1 && dy & 1) || (!(dx & 1) && !(dy & 1) && causal);
      else ok = ((dx + dy) & 1) || causal;
      if (ok) t[n][0] = d2, t[n][1] = dy, t[n][2] = dx, n++;
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
  double *w = Wm, sx[NI];
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
    double pp = p * .9999 + .00005;  // the probability the coder would use
    // INPAINT: no coder. A sampling child draws the bit, otherwise it is the true bit.
    int y1 = smp >= 0 ? uni() < pp : v >> e & 1;
    gbits[g] -= log2(y1 ? pp : 1 - pp);
    if (learn) {  // INPAINT
      for (i = 0; i < NI; i++) w[i] += LR * (y1 - p) * sx[i];
      for (i = 0; i < NC; i++) j = H[i], N[j] += N[j] < LIM, P[j] += (y1 - P[j]) / (N[j] + .5);
    }
    c = c * 2 + y1, lov += y1 << e;
  }
  v = lov;
  X[k] = v;  // in a sampling child this writes the drawn pixel
  gcount[g]++;
  double err = v - m;
  E[k] = fabs(err);
  F[k] = v, K[k] = 1;
  if (!learn) return;  // INPAINT
  V[gv][q] += (err * err - V[gv][q]) / VR;
  double yv = v - base;
  for (i = 0; i < L; i++) {
    double r = Bv[gr][i] = Bv[gr][i] * LAM + o[i] * yv;
    for (j = 0; j < L; j++) C[gr][i][j] = C[gr][i][j] * LAM + o[i] * o[j];
    for (j = 0; j < L; j++) if (j != i) r -= C[gr][i][j] * A[gr][j];
    A[gr][i] = r / (C[gr][i][i] + NU);
  }
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

static void reset(void) {  // the whole model, once
  memset(N, 0, sizeof N), memset(C, 0, sizeof C), memset(Bv, 0, sizeof Bv), memset(A, 0, sizeof A);
  memset(Wm, 0, sizeof Wm);
  for (int i = 0; i < 1 << HB; i++) P[i] = .5;
  for (int g = 0; g < NG; g++) for (int q = 0; q < 32; q++) V[g][q] = 1;
}

// INPAINT: per image: known values, known flags, errors and the cost tally
static void reset_image(void) {
  memset(K, 0, sizeof K), memset(E, 0, sizeof E), memset(gbits, 0, sizeof gbits), memset(gcount, 0, sizeof gcount);
  for (int i = 0; i < NPIX; i++) F[i] = 8;
}

// INPAINT: lab/anchor.c's run(), split into the anchor stage and the refinement stage
static void anchors(void) {
  for (int y = 0; y < W; y += T)
    for (int x = 0; x < W; x += T) code_pixel(x, y, 0, &SA, T, T);
  fill(T);
}

static void refine(void) {
  for (int t = T, lv = 0; t >= 2; t /= 2, lv++) {
    int h = t / 2, g = 1 + 3 * lv;
    for (int y = h; y < W; y += t)
      for (int x = h; x < W; x += t) code_pixel(x, y, g, &SD, h, h);
    for (int y = 0; y < W; y += h)
      for (int x = y % t ? 0 : h; x < W; x += t) code_pixel(x, y, g + 1 + (y % t != 0), &SX, h, h);
    fill(h);
  }
}

static double total(int from) {  // bits of groups from..NG-1
  double s = 0;
  for (int g = from; g < NG; g++) s += gbits[g];
  return s;
}

static void load(const char *fn) {
  static unsigned char raw[NPIX / 2];
  FILE *f = fopen(fn, "rb");
  if (!f || fread(raw, 1, NPIX / 2, f) != NPIX / 2) fprintf(stderr, "cannot read %s\n", fn), exit(1);
  fclose(f);
  for (int i = 0; i < NPIX; i++) X[i] = raw[i / 2] >> (i % 2 ? 0 : 4) & 15;
}

int main(int argc, char **argv) {
  if (argc < 5) return fprintf(stderr, "usage: %s TARGET NSAMPLES SEED OUTDIR [TRAIN ...]\n", argv[0]), 2;
  int ns = atoi(argv[2]);
  unsigned long long seed = strtoull(argv[3], 0, 10);
  long drawn = NPIX - (long)((W - 1) / T + 1) * ((W - 1) / T + 1);  // pixels that are not anchors
  mksten(&SA, 0, RA), mksten(&SD, 1, RD), mksten(&SX, 2, RX);
  reset();
  for (int n = 5; n < argc; n++) {  // training images: coded completely, learning on
    load(argv[n]), reset_image(), anchors(), refine();
    fprintf(stderr, "train %s: %.0f bytes\n", argv[n], total(0) / 8);
  }
  load(argv[1]), reset_image(), anchors();  // the target's anchors, learning on
  double abits = total(0);
  learn = 0;
  fflush(stdout);
  for (int n = 0; n <= ns; n++) {
    pid_t pid = fork();
    if (pid < 0) perror("fork"), exit(1);
    if (pid) { waitpid(pid, 0, 0); continue; }
    // child n < ns: one sample; child ns: the frozen model over the true pixels
    if (n < ns) smp = n, rng = seed << 20 | n;
    refine();
    if (n == ns) {
      printf("true-frozen\t-\t%.4f\n", total(1) / drawn);
      return 0;
    }
    static unsigned char packed[NPIX / 2];
    for (int i = 0; i < NPIX / 2; i++) packed[i] = X[2 * i] << 4 | X[2 * i + 1];
    char fn[4096];
    snprintf(fn, sizeof fn, "%s/sample_%03d.bin", argv[4], n);
    FILE *f = fopen(fn, "wb");
    if (!f || fwrite(packed, 1, sizeof packed, f) != sizeof packed || fclose(f)) return perror(fn), 1;
    printf("%d\t%llu\t%.4f\n", n, seed << 20 | n, total(1) / drawn);
    return 0;
  }
  learn = 1, refine();  // the true pixels with learning on, as the compressor codes them
  printf("true-learning\t-\t%.4f\n", total(1) / drawn);
  fprintf(stderr, "target: anchors %.0f bytes, total %.0f bytes (ideal, as the compressor)\n", abits / 8, total(0) / 8);
  return 0;
}
