// Object fitter (lab tool). Reads a 125000-byte image (image.bin format) from stdin, finds
// objects of the generator's kinds and prints them in render.c's format, one per line:
//   K cx cy R th E A Q   (K 0 ellipse/disc, 1 rect/line, 2 ring, 3 striped disc; A grey levels;
//                         Q 1 = opaque, covers the layer below)
// Method: greedy. Residual r = image - layer (layer = the objects found so far, composed in
// order). h_s(f) = f minus its mean over the (2s+1)^2 window (clipped to the image), a high-pass.
// An object's gain at scale s = the drop of sum h_s(r)^2 when it is added on top with the best
// amplitude (closed form), as additive or as opaque, whichever drops more. Seeds = the pixel
// with most h_1 energy (5x5) not tried yet; proposals of each kind are placed through the seed
// using the local orientation (structure tensor) and scored at s = 4; the best NREF of each kind
// are refined by coordinate descent at s = 4, 2, 1 (step halving), and the best refined object
// is kept if its gain at s = 1 exceeds GMIN.
// Knobs: -DGMIN=.. -DNMAX=.. -DFAILS=.. (consecutive failed seeds before stopping) -DNREF=..
// Diagnostic: ./fit TRUE.txt <image scores the true objects in order and refines each from
// the truth (or from perturbed parameters with -DPERT=n).
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef GMIN
#define GMIN 60
#endif
#ifndef NMAX
#define NMAX 80
#endif
#ifndef FAILS
#define FAILS 30
#endif
#ifndef SMAX
#define SMAX 8
#endif
#ifndef NREF
#define NREF 2
#endif
#define W 500
typedef double D;
unsigned char B[125000];
D P[W * W], Ly[W * W], Rr[W * W], Hs[9][W * W], Mb[W * W], Ml[W * W], Hm[W * W], Hl[W * W], Sat[(W + 1) * (W + 1)];
char tried[W * W];
D ob[200][8];
int nob;
long nev;

// mask of object o at pixel (x, y), exactly as gen_image.py (c, s = cos, sin of o[4])
D msk(D *o, D c, D s, int x, int y) {
  D u = (x - o[1]) * c + (y - o[2]) * s, v = (y - o[2]) * c - (x - o[1]) * s, g = hypot(u, v), d, m;
  int K = o[0];
  d = K == 0 ? (hypot(u / o[3], v / o[5]) - 1) * o[5] : K == 1 ? fmax(fabs(u) - o[3], fabs(v) - o[5]) : K == 2 ? fabs(g - o[3]) - o[5] : g - o[3];
  m = fmin(fmax(.5 - d, 0), 1);
  return K == 3 ? m * sin(2 * M_PI * u / o[5]) : m;
}

// dst = src - mean of src over the (2s+1)^2 window clipped to the w x h buffer
void hp(D *src, D *dst, int w, int h, int s) {
  int x, y, w1 = w + 1;
  for (y = 0; y <= h; y++)
    for (x = 0; x <= w; x++) Sat[y * w1 + x] = x && y ? src[(y - 1) * w + x - 1] + Sat[(y - 1) * w1 + x] + Sat[y * w1 + x - 1] - Sat[(y - 1) * w1 + x - 1] : 0;
  for (y = 0; y < h; y++)
    for (x = 0; x < w; x++) {
      int a = x > s ? x - s : 0, b = x + s + 1 < w ? x + s + 1 : w, c = y > s ? y - s : 0, d = y + s + 1 < h ? y + s + 1 : h;
      dst[y * w + x] = src[y * w + x] - (Sat[d * w1 + b] - Sat[c * w1 + b] - Sat[d * w1 + a] + Sat[c * w1 + a]) / ((b - a) * (d - c));
    }
}

void update() {
  for (int i = 0; i < W * W; i++) Rr[i] = P[i] - Ly[i];
  for (int s = 1; s <= SMAX; s *= 2) hp(Rr, Hs[s], W, W, s);
}

int bx[4];
void box(D *o, D mg) {
  D rr = o[3] + o[5] + 1 + mg;
  bx[0] = fmax(floor(o[1] - rr), 0), bx[1] = fmin(ceil(o[1] + rr), W - 1);
  bx[2] = fmax(floor(o[2] - rr), 0), bx[3] = fmin(ceil(o[2] + rr), W - 1);
}

// gain at scale s of adding o on top; sets o[6] (amplitude) and o[7] (opaque)
D score(D *o, int s) {
  int K = o[0];
  if (o[3] < 1 || o[5] < .3 || o[3] > 400 || o[5] > 80 || (K == 0 && o[5] > o[3]) || (K == 3 && o[5] < 3)) return -1;
  nev++;
  box(o, 2 * s);
  int x0 = bx[0], y0 = bx[2], bw = bx[1] - x0 + 1, bh = bx[3] - y0 + 1, x, y;
  if (bw < 1 || bh < 1) return -1;
  box(o, s);
  D c = cos(o[4]), sn = sin(o[4]), a = 0, b = 0, rr = 0, ww = 0, wy = 0;
  for (y = 0; y < bh; y++)
    for (x = 0; x < bw; x++) {
      int i = y * bw + x;
      Mb[i] = msk(o, c, sn, x + x0, y + y0), Ml[i] = Mb[i] * Ly[(y + y0) * W + x + x0];
    }
  hp(Mb, Hm, bw, bh, s), hp(Ml, Hl, bw, bh, s);
  for (y = bx[2]; y <= bx[3]; y++)
    for (x = bx[0]; x <= bx[1]; x++) {
      int i = (y - y0) * bw + x - x0;
      D hm = Hm[i], hr = Hs[s][y * W + x], w = hr + Hl[i];
      a += hr * hm, b += hm * hm, rr += hr * hr, ww += w * w, wy += w * hm;
    }
  if (b < 1e-9) return -1;
  D ga = a * a / b, go = rr - ww + wy * wy / b;
  o[7] = K < 3 && go > ga;
  o[6] = o[7] ? wy / b : a / b;
  return o[7] ? go : ga;
}

void addobj(D *o) {
  int x, y;
  D c = cos(o[4]), s = sin(o[4]);
  box(o, 1);
  for (y = bx[2]; y <= bx[3]; y++)
    for (x = bx[0]; x <= bx[1]; x++) {
      D m = msk(o, c, s, x, y);
      Ly[y * W + x] = o[7] ? Ly[y * W + x] * (1 - m) + o[6] * m : Ly[y * W + x] + o[6] * m;
    }
  update();
  for (x = 0; x < 8; x++) ob[nob][x] = o[x];
  nob++;
}

// coordinate descent at scales 4, 2, 1: moves along the object's own axes u, v, then R, theta, E
D refine(D *o) {
  D best = 0, st[5], t[8];
  int j, k, s, h;
  for (s = SMAX; s; s /= 2) {
    best = score(o, s);
    st[0] = o[0] == 3 ? o[5] : s, st[1] = s, st[2] = s, st[3] = s / (o[3] + o[5] + 1), st[4] = o[0] == 3 ? s * o[5] / (o[3] + 1) / 2 : s;
    for (h = 0; h < 4;) {
      int imp = 0;
      for (j = 0; j < 5; j++)
        for (k = -1; k < 2; k += 2) {
          for (int q = 0; q < 8; q++) t[q] = o[q];
          D c = cos(o[4]), sn = sin(o[4]), dd = k * st[j];
          if (j == 0) t[1] += dd * c, t[2] += dd * sn;
          if (j == 1) t[1] -= dd * sn, t[2] += dd * c;
          if (j == 2) t[3] += dd, t[5] += o[0] == 0 && o[3] == o[5] ? dd : 0;
          if (j == 3) t[4] += dd;
          if (j == 4) t[5] += dd;
          D g = score(t, s);
          if (g > best) {
            best = g, imp = 1;
            for (int q = 0; q < 8; q++) o[q] = t[q];
          }
        }
      if (!imp) {
        for (j = 0; j < 5; j++) st[j] /= 2;
        h++;
      }
    }
  }
  return score(o, 1);
}

int main(int argc, char **argv) {
  int i, x, y, fails = 0;
  if (fread(B, 1, 125000, stdin) != 125000) return 1;
  for (i = 0; i < W * W; i++) P[i] = B[i / 2] >> (4 - i % 2 * 4) & 15;
  update();
  if (argc > 1) {  // diagnostic: score the true objects (render.c format) in order, refine each
    FILE *f = fopen(argv[1], "r");
    D o[8], r[8];
    int K, Q;
    while (fscanf(f, "%d %lf %lf %lf %lf %lf %lf %d", &K, o + 1, o + 2, o + 3, o + 4, o + 5, o + 6, &Q) == 8) {
      o[0] = K;
      D A = o[6];
      D g = score(o, 1);
      for (i = 0; i < 8; i++) r[i] = o[i];
#ifdef PERT  // start from perturbed parameters
      r[1] += PERT * 1.5, r[2] -= PERT, r[3] *= 1 + PERT * .05, r[4] += PERT * .03, r[5] *= K == 0 && r[3] == r[5] ? 1 + PERT * .05 : 1 - PERT * .05;
#endif
      D gr = refine(r);
      fprintf(stderr, "K%d X %.1f Y %.1f R %.1f E %.1f A %.2f Q %d | gain %.0f A %.2f Q %d | refined %.0f dX %.2f dY %.2f dR %.2f dth %.3f dE %.2f\n", K, o[1], o[2], o[3], o[5], A, Q, g, o[6], (int)o[7], gr, r[1] - o[1], r[2] - o[2], r[3] - o[3], r[4] - o[4], r[5] - o[5]);
      o[6] = A, o[7] = Q;
      addobj(o);
    }
    return 0;
  }
  while (nob < NMAX && fails < FAILS) {
    // seed: most high-pass energy in a 5x5 window, not tried yet
    D bs = -1;
    int sp = 0;
    for (y = 3; y < W - 3; y++)
      for (x = 3; x < W - 3; x++) {
        i = y * W + x;
        if (tried[i]) continue;
        D e = 0;
        for (int dy = -2; dy < 3; dy++)
          for (int dx = -2; dx < 3; dx++) e += Hs[1][i + dy * W + dx] * Hs[1][i + dy * W + dx];
        if (e > bs) bs = e, sp = i;
      }
    int px = sp % W, py = sp / W;
    for (y = -3; y < 4; y++)
      for (x = -3; x < 4; x++) tried[sp + y * W + x] = 1;
    // orientation of the residual's gradient around the seed (structure tensor, 7x7)
    D jxx = 0, jxy = 0, jyy = 0;
    for (y = py - 3; y <= py + 3; y++)
      for (x = px - 3; x <= px + 3; x++) {
        i = y * W + x;
        D gx = Rr[i + 1] - Rr[i - 1], gy = Rr[i + W] - Rr[i - W];
        jxx += gx * gx, jxy += gx * gy, jyy += gy * gy;
      }
    D al = atan2(2 * jxy, jxx - jyy) / 2, nc = cos(al), ns = sin(al);
    D bo[4][NREF][8], bg[4][NREF], t[8];
    for (int K = 0; K < 4; K++) {
      for (int q = 0; q < NREF; q++) bg[K][q] = -1;
      int npr = K == 0 ? 28 : K == 1 ? 28 : K == 2 ? 168 : 112;
      for (int j = 0; j < npr; j++) {
        int sg = j % 2 * 2 - 1;
        t[0] = K, t[4] = al;
        if (K == 0) t[3] = t[5] = 8 + j / 2 * 4, t[1] = px - sg * t[3] * nc, t[2] = py - sg * t[3] * ns;
        if (K == 1) {
          static const D es[7] = {1, 2, 4, 8, 16, 32, 50};
          t[5] = es[j / 4], t[3] = j / 2 % 2 ? 40 : 15, t[4] = al + M_PI / 2;
          if (t[3] < t[5]) t[3] = t[5];
          t[1] = px - sg * t[5] * nc, t[2] = py - sg * t[5] * ns;
        }
        if (K == 2) {
          t[3] = 8 + j / 12 * 4, t[5] = 2 + j / 4 % 3 * 2;
          D rr = t[3] + (j / 2 % 2 ? t[5] : -t[5]);
          t[1] = px - sg * rr * nc, t[2] = py - sg * rr * ns;
        }
        if (K == 3) {
          t[3] = 15, t[5] = 5 + j / 8;
          D ph = (j % 8) * t[5] / 8;
          t[1] = px + ph * nc, t[2] = py + ph * ns;
        }
        D g = score(t, SMAX);
        // keep the NREF best proposals of this kind
        int w = -1;
        for (int q = 0; q < NREF; q++)
          if (g > bg[K][q] && (w < 0 || bg[K][q] < bg[K][w])) w = q;
        if (w >= 0) {
          bg[K][w] = g;
          for (int q = 0; q < 8; q++) bo[K][w][q] = t[q];
        }
      }
    }
    D *bb = 0, bgr = GMIN;
    for (int K = 0; K < 4; K++)
      for (int q = 0; q < NREF; q++) {
        if (bg[K][q] <= 0) continue;
        D g = refine(bo[K][q]);
        if (g > bgr) bgr = g, bb = bo[K][q];
      }
    if (!bb) {
      fails++;
      continue;
    }
    fails = 0;
    addobj(bb);
    fprintf(stderr, "#%d seed (%d,%d) K%d gain %.0f  X %.2f Y %.2f R %.2f th %.3f E %.2f A %.2f Q %d  evals %ld\n", nob, px, py, (int)bb[0], bgr, bb[1], bb[2], bb[3], bb[4], bb[5], bb[6], (int)bb[7], nev);
  }
  for (i = 0; i < nob; i++) {
    D *o = ob[i];
    D th = fmod(o[4], M_PI), A = o[6];
    if (th < 0) th += M_PI;
    if (o[0] == 3 && labs(lround((o[4] - th) / M_PI)) % 2) A = -A;  // stripes flip under theta + pi
    printf("%d %.4f %.4f %.4f %.6f %.4f %.4f %d\n", (int)o[0], o[1], o[2], o[3], th, o[5], A, (int)o[7]);
  }
  fprintf(stderr, "%d objects, %ld evaluations\n", nob, nev);
}
