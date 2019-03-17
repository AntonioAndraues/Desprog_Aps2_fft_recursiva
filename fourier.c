#include "fourier.h"
#include <math.h>
#include <stdio.h>
void nft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
  for (int k = 0; k < n; k++) {
    t[k] = 0;

    for (int j = 0; j < n; j++) {
      t[k] += s[j] * cexp(sign * 2 * PI * k * j * I / n);
    }
  }
}

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
  nft(s, t, n, -1);
}

void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
  nft(t, s, n, 1);

  for (int k = 0; k < n; k++) {
    s[k] /= n;
  }
}

void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {

  if (n == 1) {
    t[0] = s[0];
    return;
  }
  if (n==2){
    t[0]=s[0]+s[1]* cexp(-2 * PI * 0 * I / n);
    t[1]=s[0]-s[1]* cexp(-2 * PI * 0 * I / n);
  }


  int metade_de_n = n / 2;
  double complex s_p[metade_de_n];
  double complex s_i[metade_de_n];
  double complex t_p[metade_de_n];
  double complex t_i[metade_de_n];
  // for (int k = 0; k <= metade_de_n - 1; k++) {
  //   t_p[k]=0;
  //   t_i[k]=0;
  // }

  for (int j = 0; j < metade_de_n - 1; j++) {
    s_p[j] = s[2 * j];
    s_i[j] = s[2 * j + 1];
  }

  fft(s_p, t_p, metade_de_n, sign);
  fft(s_i, t_i, metade_de_n, sign);


  for (int k = 0; k <= metade_de_n - 1; k++) {
    t[k]=0;
    t[k+metade_de_n]=0;
    t[k] = t_p[k] + t_i[k] * cexp(-2 * PI * k * I / n);
    t[k + metade_de_n] = t_p[k] - t_i[k] * cexp(-2 * PI * k * I / n);
  }
}

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
  fft(s, t, n, -1);
}

void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
  fft(t, s, n, 1);

  for (int k = 0; k < n; k++) {
    s[k] /= n;
  }
}

void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
}

void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
}

void filter(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
  int center_x = width / 2;
  int center_y = height / 2;

  double variance = -2 * SIGMA * SIGMA;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int dx = center_x - (x + center_x) % width;
      int dy = center_y - (y + center_y) % height;

      double d = dx * dx + dy * dy;

      double g = exp(d / variance);

      if (flip) {
        g = 1 - g;
      }

      output[y][x] = g * input[y][x];
    }
  }
}

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
  filter(input, output, width, height, 0);
}

void filter_hp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
  filter(input, output, width, height, 1);
}
