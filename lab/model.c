// Research version of the raster reference model with the model ideas of ideas.md as -D knobs
// (readable lab tool, not a submission). Without knobs it is candidates/000_raster_ref/src.c and
// writes the same archive.
//
// The reference: 5 hashed counters and one Gaussian input centred on an online least-squares
// (RLS) prediction from 112 causal neighbours, mixed by one logistic mixer, binary arithmetic
// coder. Knobs add further Gaussian inputs, each with its own error image and variance buckets:
//   -DKRG       idea 1: a fixed predictor computed from the generator's spectrum (simple kriging
//               weights for the same 112 neighbours, scale KS grey levels per unit, default 2.5)
//   -DNF=n      idea 2: n small RLS predictors that forget fast; predictor i uses the FN_i nearest
//               causal neighbours and forgetting FL_i (defaults below)
//   -DSEL=r     idea 3: the slow RLS learns a pixel with weight r (0..1) where the fastest small
//               predictor did better than it on the 4 nearest pixels, else with weight 1
//   -DTM=w      idea 4: template matching, best match of the 6 nearest neighbours within w
//               columns either side and TMR rows above (default 12)
//               -DTMN=n template size, -DTMK=k average of the k best matches, -DTMO match
//               the shape (mean difference removed and added back), -DTCX a counter in the
//               context (best match's value, match distance, W), -DTMH=h non-local means:
//               all candidates averaged with weights exp(-distance / h)
//   -DMIXC=1    mixer weight set chosen by which Gaussian input has the smallest local error
//   -DMIXE=1    and by the bit position (4 sets per choice)
//   -DKP=x      the slow RLS is regularised towards the kriging weights with strength x (in
//               units of its diagonal sums; the reference regularises towards 0 with 100)
//   -DKFIX=sd   the kriging input's sd fixed (no buckets)
//   -DKCX=p     a 6th counter in the context (Gaussian input p's prediction in half levels,
//               its error bucket), as the quadtree's (m, q) context
//   -DVQ=x      Gaussian variance V - x/2 for the added inputs (x = 1/6 removes the rounding
//               counted twice), floored at the sensor noise
// usage: model <image.bin >archive   (compress)   model <archive >image   (decompress)
//        -DCOST: compressing also writes cost.f32 (bits per image pixel) for lab/where_bits.py
#include <math.h>
#include <stdio.h>
#include <string.h>
#define F(i, n) for (i = 0; i < n; i++)
#define S 528  // row stride of the padded image: 8-pixel left margin, 20-pixel right margin
#define L 112  // slow RLS inputs: first 112 cells (row-major) of the 8x15 block ending at k
typedef double D;

#ifndef MLR
#define MLR .001  // mixer learning rate
#endif
#ifndef KS
#define KS 2.5
#endif
#ifndef NF
#define NF 0
#endif
#ifndef TMN
#define TMN 6  // template: the TMN nearest causal neighbours
#endif
#ifndef TMK
#define TMK 1  // average of the TMK best matches
#endif
#ifndef TMR
#define TMR 12
#endif
#ifdef KRG
#define NK 1
#else
#define NK 0
#endif
#ifdef TM
#define NT 1
#else
#define NT 0
#endif
#define NG (1 + NK + NF + NT)  // Gaussian inputs
#if defined(KCX) || defined(TCX)
#define NC 6  // counters
#else
#define NC 5
#endif
#define NI (NC + NG)  // mixer inputs
#ifdef MIXC
#define NWS0 NG
#else
#define NWS0 1
#endif
#ifdef MIXE
#define NWS (NWS0 * 4)
#else
#define NWS NWS0
#endif
#ifndef FN0
#define FN0 12
#endif
#ifndef FN1
#define FN1 6
#endif
#ifndef FN2
#define FN2 24
#endif
#ifndef FN3
#define FN3 4
#endif
#ifndef FL0
#define FL0 .99
#endif
#ifndef FL1
#define FL1 .95
#endif
#ifndef FL2
#define FL2 .995
#endif
#ifndef FL3
#define FL3 .9
#endif
static const int FN[4] = {FN0, FN1, FN2, FN3};
static const D FL[4] = {FL0, FL1, FL2, FL3};

unsigned char B[1 << 17], O[125000], I[S * S], N[1 << 22];
unsigned x, z = -1, X, M;
int h[6], H[NI + 1], d, u, i, j, k, y, c, e, r, l, v, q, b;
D P[1 << 22], E[S * S], w[NWS][NI], s[NI], o[L + 1], a[L + 1], C[L][L + 1], V[32], m, g, t, f, p;
// added Gaussian inputs: prediction, sd, error image, variance buckets, bucket
D mu[NG], gs[NG], Er[NG][S * S], Vr[NG][32];
int qr[NG], nb[24];  // nb: offsets (in o[] indices) of the causal cells, nearest first
// small fast RLS predictors
D fa[4][25], fC[4][24][25], fu[25];
// kriging weights
D kw[L];
int tmv, tms;  // best match's value and its (capped) squared distance
#ifdef COST
float Q[250000];
#endif

D st(D p) { return p = p * .9998 + .0001, log(p / (1 - p)); }
D G(int v) { return v <<= e, v < 17 ? 0 : v > 31 ? 1 : erfc((m + 16.5 - v) / g) / 2; }

// kriging weights for the 112-cell stencil from the generator's background spectrum
void kriging(void) {
  static D Sp[250000], R[15][29], A[L][L + 1];
  static const D be[3] = {3.4, 2.2, 1.0}, wt[3] = {1, .09, .0064};
  D nz[3] = {0}, fr;
  int n = 500, a1, b1, ky, kx;
  F(i, n * n) if (i) {
    ky = i / n, kx = i % n, ky -= ky > 249 ? n : 0, kx -= kx > 249 ? n : 0;
    fr = hypot(ky, kx) / n;
    F(j, 3) nz[j] += pow(fr, -be[j]);
  }
  F(i, n * n) if (i) {
    ky = i / n, kx = i % n, ky -= ky > 249 ? n : 0, kx -= kx > 249 ? n : 0;
    fr = hypot(ky, kx) / n, Sp[i] = 0;
    F(j, 3) Sp[i] += wt[j] * pow(fr, -be[j]) / nz[j] * n * n;
    Sp[i] *= KS * KS;
  }
  F(a1, 15) F(b1, 29) {
    D t2 = 0;
    F(i, n * n) if (i) {
      ky = i / n, kx = i % n, ky -= ky > 249 ? n : 0, kx -= kx > 249 ? n : 0;
      t2 += Sp[i] * cos(2 * M_PI * (ky * (a1 - 7) + kx * (b1 - 14)) / n);
    }
    R[a1][b1] = t2 / n / n;
  }
  F(i, L) {
    F(j, L) A[i][j] = R[i / 15 - j / 15 + 7][i % 15 - j % 15 + 14] + (i == j) * (.0784 + 1 / 12.);
    A[i][L] = R[i / 15 - 7 + 7][i % 15 - 7 + 14];
  }
  F(i, L) for (j = i + 1; j < L; j++) {  // Gaussian elimination, then back substitution
    D t2 = A[j][i] / A[i][i];
    for (int k2 = i; k2 <= L; k2++) A[j][k2] -= t2 * A[i][k2];
  }
  for (i = L; i--;) {
    D t2 = A[i][L];
    for (j = i + 1; j < L; j++) t2 -= A[i][j] * kw[j];
    kw[i] = t2 / A[i][i];
  }
  D sw = 0;
  F(i, L) sw += kw[i];
  F(i, L) kw[i] /= sw;  // weights sum to 1: the local mean is unknown
#ifdef KDUMP
  F(i, L) fprintf(stderr, "%.6f\n", kw[i]);
#endif
}

// error-bucket sd for Gaussian input p at padded position k
D bucket(int p) {
  qr[p] = fmin((Er[p][k - 1] + Er[p][k - S] + Er[p][k - S - 1] + Er[p][k - S + 1]) * 3, 31);
#ifdef VQ
  return sqrt(fmax(Vr[p][qr[p]] * 2 - VQ, .157));
#else
  return sqrt(Vr[p][qr[p]] * 2);
#endif
}

int main() {
  d = fread(B, 1, 1 << 17, stdin) != 125000;
  F(i, 1 << 22) P[i] = .5, I[i % (S * S)] = 8, V[i % 32] = 1, a[L] = -1;
  F(i, NG) F(j, 32) Vr[i][j] = 1;
  F(i, 4) fa[i][FN[i]] = -1;
  {  // causal cells of the block, nearest first (ties: row-major)
    int used[L] = {0}, best;
    F(i, 24) {
      best = -1;
      F(j, L) if (!used[j]) {
        int dy = j / 15 - 7, dx = j % 15 - 7, bb;
        if (best < 0) best = j;
        bb = best;
        if (dy * dy + dx * dx < (bb / 15 - 7) * (bb / 15 - 7) + (bb % 15 - 7) * (bb % 15 - 7)) best = j;
      }
      used[best] = 1, nb[i] = best;
    }
  }
#if defined(KRG) || defined(KP)
  kriging();
#endif
  F(i, 4) X = X << 8 | B[u++];
  F(r, 500) F(l, 500) {
    k = (r + 8) * S + l + 8;
    if (!l) F(j, 8) I[k - j - 1] = I[k - S];
    m = 0;
    F(e, L) m += a[e] * (o[e] = I[k - 7 * S - 7 + e / 15 * S + e % 15]);
    q = fmin((E[k - 1] + E[k - S] + E[k - S - 1] + E[k - S + 1]) * 3, 31);
    g = sqrt(V[q] * 2);
    mu[0] = m, gs[0] = g;
    int pi = 1;
#ifdef KRG
    mu[pi] = 0;
    F(i, L) mu[pi] += kw[i] * o[i];
    gs[pi] = bucket(pi), pi++;
#ifdef KFIX  // the kriging input's sd fixed at its theoretical value for background
    gs[pi - 1] = KFIX * sqrt(2);
#endif
#endif
    F(j, NF) {
      mu[pi] = 0;
      F(i, FN[j]) mu[pi] += fa[j][i] * o[nb[i]];
      gs[pi] = bucket(pi), pi++;
    }
#ifdef TM
    {  // best matches of the TMN nearest neighbours among causal positions
      D bss[TMK], ss, off[TMK], om, nlw = 0, nlv = 0;
      int bkk[TMK];
      F(i, TMK) bss[i] = 1e9, bkk[i] = k - 1, off[i] = 0;
      for (int dr = -TMR; dr <= 0; dr++)
        for (int dc = -TM; dc <= TM; dc++) {
          int kk = k + dr * S + dc;
          if (dr == 0 && dc >= 0) break;
          if (r + dr < 2 || l + dc < 2 || l + dc > 499) continue;
          ss = 0, om = 0;
          F(i, TMN) {
            D df = o[nb[i]] - I[kk + (nb[i] / 15 - 7) * S + nb[i] % 15 - 7];
            om += df;
#ifdef TMO
            ss += df * df;
#else
            ss += df * df;
#endif
          }
#ifdef TMO  // match on the shape: the mean difference is removed and added back
          ss -= om * om / TMN;
#endif
#ifdef TMH  // non-local means: every candidate weighted by exp(-distance / TMH)
          nlw += exp(-ss / TMH), nlv += exp(-ss / TMH) * I[kk];
#endif
          for (j = TMK; j-- && ss <= bss[j];)  // keep the TMK best, nearest rows win ties
            if (j + 1 < TMK) bss[j + 1] = bss[j], bkk[j + 1] = bkk[j], off[j + 1] = off[j];
          if (++j < TMK) bss[j] = ss, bkk[j] = kk, off[j] = om / TMN;
        }
      mu[pi] = 0;
      F(i, TMK) mu[pi] += I[bkk[i]]
#ifdef TMO
                          + off[i]
#endif
          ;
      mu[pi] /= TMK;
#ifdef TMH
      mu[pi] = nlw > 0 ? nlv / nlw : I[k - 1];
#endif
      tmv = I[bkk[0]], tms = fmin(bss[0], 15);
      gs[pi] = bucket(pi), pi++;
#ifdef TSHARE  // the template input uses the RLS input's variance bucket
      gs[pi - 1] = g;
#endif
    }
#endif
    int ws = 0;
#ifdef MIXC
    {
      D be2 = 1e9, le;
      F(i, NG) {
        le = i ? Er[i][k - 1] + Er[i][k - S] + Er[i][k - S - 1] + Er[i][k - S + 1]
               : E[k - 1] + E[k - S] + E[k - S - 1] + E[k - S + 1];
        if (le < be2) be2 = le, ws = i;
      }
    }
#endif
    h[0] = I[k - 1] << 8 | I[k - 2] << 4 | I[k - S + 2];
    h[1] = I[k - S] << 8 | I[k - S + 1] << 4 | I[k - 2 * S + 1];
    h[2] = I[k - 1] << 8 | I[k - S] << 4 | I[k - S - 1];
    h[3] = h[1] << 12 | h[2];
    h[4] = I[k - 3] << 8 | I[k - S + 3] << 4 | I[k - S + 1];
#ifdef KCX
    h[5] = (int)(mu[KCX] * 2 + .5) << 5 | qr[KCX];
#endif
#ifdef TCX  // counter context: the best match's value and how well it matched
    h[5] = tmv << 8 | tms << 4 | I[k - 1];
#endif
    b = (499 - l) * 500 + r;
    v = B[b / 2] >> 4 - b % 2 * 4 & 15;
    for (c = 1, e = 4; e--; c += y) {
      c *= 2;
      f = 0;
      F(j, NG) {
        m = mu[j], g = gs[j];
        t = G(c + 2);
        s[NC + j] = st((t - G(c + 1)) / (t - G(c) + 1e-30));
      }
      F(i, NC) s[i] = st(P[H[i] = (h[i] << 5 | c | i << 29) * 2654435761u >> 10]);
#ifdef NORLS  // the RLS prediction is not a mixer input (it still sets the buckets)
      s[NC] = 0;
#endif
#ifdef ONLY  // diagnostic: the mixer sees only input ONLY
      F(i, NI) if (i != ONLY) s[i] = 0;
#endif
#ifdef MIXE
      int wsel = ws * 4 + e;
#else
      int wsel = ws;
#endif
      F(i, NI) f += w[wsel][i] * s[i];
      p = 1 / (1 + exp(-f));
      M = x + (z - x) * (p * .9999 + .00005);
      y = d ? X <= M : v >> e & 1;
#ifdef COST
      Q[b] -= log2(y ? p * .9999 + .00005 : 1 - (p * .9999 + .00005));
#endif
      y ? z = M : (x = M + 1);
      for (; (x ^ z) < 1 << 24; x <<= 8, z = z << 8 | 255) d ? X = X << 8 | B[u++] : putchar(z >> 24);
      F(i, NI) w[wsel][i] += MLR * (y - p) * s[i];
      F(i, NC + 1) j = H[i], N[j] += N[j] < 80, P[j] += (y - P[j]) / (N[j] + .5);
    }
    m = mu[0];
    I[k] = c -= 16;
    O[b / 2] |= c << 4 - b % 2 * 4;
    if (l > 498) F(j, 20) I[k + j] = c;
    E[k] = fabs(c - m);
    V[q] += (E[k] * E[k] - V[q]) / 500;
    for (j = 1; j < NG; j++) Er[j][k] = fabs(c - mu[j]), Vr[j][qr[j]] += (Er[j][k] * Er[j][k] - Vr[j][qr[j]]) / 500;
    // small fast RLS predictors: update sums with forgetting, one Gauss-Seidel sweep
    F(j, NF) {
      int n = FN[j];
      F(i, n) fu[i] = o[nb[i]];
      fu[n] = c;
      F(i, n) {
        D mm = 0;
        for (int jj = 0; jj <= n; jj++) fC[j][i][jj] = fC[j][i][jj] * FL[j] + fu[i] * fu[jj], mm -= jj - i ? fC[j][i][jj] * fa[j][jj] : 0;
        fa[j][i] = mm / (fC[j][i][i] + 1);
      }
    }
    D rho = 1;
#if defined(SEL) && NF > 0
    {  // learn less where the fastest small predictor beat the slow one nearby
      int pf = 1 + NK + NF - 1;
      if (Er[pf][k - 1] + Er[pf][k - S] + Er[pf][k - S - 1] + Er[pf][k - S + 1] <
          E[k - 1] + E[k - S] + E[k - S - 1] + E[k - S + 1])
        rho = SEL;
    }
#endif
    o[L] = c;
    if (rho < 1) F(i, L + 1) o[i] *= sqrt(rho);
    F(i, L) {
      m = 0;
      F(j, L + 1) C[i][j] = C[i][j] * .9997 + o[i] * o[j], m -= j - i ? C[i][j] * a[j] : 0;
#ifdef KP
      a[i] = (m + KP * kw[i]) / (C[i][i] + 100 + KP);
#else
      a[i] = m / (C[i][i] + 100);
#endif
    }
  }
  d ? fwrite(O, 1, 125000, stdout) : putchar((x >> 24) + 1);
#ifdef COST
  d || fwrite(Q, 4, 250000, fopen("cost.f32", "wb"));
#endif
}
