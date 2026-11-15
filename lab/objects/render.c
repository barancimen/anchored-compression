// Render an object list into a layer (lab tool). Reads lines "K cx cy R th E A Q" from stdin
// (K: 0 ellipse or disc (semi-axes R >= E), 1 rectangle or line (half-sides R, E), 2 ring
// (radius R, half-width E), 3 striped disc (radius R, period E); A amplitude in grey levels;
// Q 1 = opaque: covers what is below), writes float32 grey levels (row-major) to stdout.
// -DNOCOVER treats every object as additive. -DQB=n rounds positions and sizes to 1/n px,
// angles to pi/(64 n), amplitudes to 1/16.
#include <stdio.h>
#include <math.h>
double T[250000];
float o[250000];
int main() {
  int K, Q, x, y;
  double X, Y, R, th, E, A;
  while (scanf("%d %lf %lf %lf %lf %lf %lf %d", &K, &X, &Y, &R, &th, &E, &A, &Q) == 8) {
#ifdef QB
    X = round(X * QB) / QB, Y = round(Y * QB) / QB, R = round(R * QB) / QB, E = round(E * QB) / QB;
    th = round(th / M_PI * 64 * QB) * M_PI / 64 / QB, A = round(A * 16) / 16;
#endif
    double c = cos(th), s = sin(th);
    for (y = 0; y < 500; y++)
      for (x = 0; x < 500; x++) {
        double u = (x - X) * c + (y - Y) * s, v = (y - Y) * c - (x - X) * s, g = hypot(u, v), d, m;
        d = K == 0 ? (hypot(u / R, v / E) - 1) * E : K == 1 ? fmax(fabs(u) - R, fabs(v) - E) : K == 2 ? fabs(g - R) - E : g - R;
        m = fmin(fmax(.5 - d, 0), 1);
        if (K == 3) m *= sin(2 * M_PI * u / E);
#ifdef NOCOVER
        T[y * 500 + x] += A * m;
#else
        T[y * 500 + x] = Q ? T[y * 500 + x] * (1 - m) + A * m : T[y * 500 + x] + A * m;
#endif
      }
  }
  for (x = 0; x < 250000; x++) o[x] = T[x];
  fwrite(o, 4, 250000, stdout);
}
