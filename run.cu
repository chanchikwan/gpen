#include <stdio.h>
#include <unistd.h>
#include "gpen.h"

Q f0(R x, R y, R z)
{
  Q f;

  x -= 0.5;
  y -= 0.5;
  z -= 0.5;

  f.lnrho = -0.5 * (x * x + y * y + z * z) / 0.01;
  f.ux    = 0.0;
  f.uy    = 0.0;
  f.uz    = 0.0;

  return f;
}

int main(int argc, char *argv[])
{
  const char rotor[] = "-/|\\";

  const R tt = (argc > 1) ? atof(argv[1]) : 1.0;
  const Z nt = (argc > 2) ? atoi(argv[2]) : 100;
  const Z nx = (argc > 3) ? atoi(argv[3]) : 256;
  const Z ny = (argc > 4) ? atoi(argv[4]) :  nx;
  const Z nz = (argc > 5) ? atoi(argv[5]) :  ny;

  const R fo = 0.0;    /* TODO: number of floating-point operations */
  const R dt = 1.0e-3; /* TODO: compute from velocity */

  Q *f = NULL;
  Z  i = 0;
  cudaEvent_t t0, t1;
  cudaEventCreate(&t0);
  cudaEventCreate(&t1);

  printf("G-Pen: reimplementing the pencil code for GPU\n");

  f = initialize_modules(nx, ny, nz);

  initial_condition(f, f0);

  while(output(i++, f) < nt) {
    const Z ns = (Z)ceilf(tt / nt / dt);
    const R ds = tt / nt / ns;

    Z j = 0;
    float ms;

    printf("%4d: %5.2f -> %5.2f: dt ~ %.0e:       ",
           i, ds * ns * (i-1), ds * ns * i, ds);

    cudaEventRecord(t0, 0);

    while(j++ < ns) {
      printf("\b\b\b\b\b\b%c %4d", rotor[j%4], j);
      fflush(stdout);

      usleep(10000); /* TODO: time stepping */

      cudaThreadSynchronize();
    }
    cudaEventRecord(t1, 0);

    cudaEventSynchronize(t1);
    cudaEventElapsedTime(&ms, t0, t1); ms /= ns;
    printf("\b\b\b\b\b\b%.3f ms/cycle ~ %.3f GFLOPS\n",
           ms, 1.0e-6 * fo / ms);
  }

  cudaEventDestroy(t1);
  cudaEventDestroy(t0);
  return 0;
}
