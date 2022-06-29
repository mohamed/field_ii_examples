#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void calc_delay(double c,
                double s,
                double * dx,
                int D,
                double * r,
                int R,
                double * out) {
  double r_x, r_y;
  for (int i = 0; i < R; i++) {
    r_x = pow(r[i] * c, 2);
    for (int j = 0; j < D; j++) {
      r_y = pow(r[i] * s - dx[j], 2);
      out[i*D + j] = r[i] + sqrt(r_x + r_y);
    }
  }
}
