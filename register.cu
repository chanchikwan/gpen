#include "gpen.h"

extern void initialize_io               (void *, const Z, const Z, const Z);
extern void initialize_initial_condition(void *, const Z, const Z, const Z);
extern void initialize_rk_2n            (void *, const Z, const Z, const Z);
extern void initialize_pde              (        const Z, const Z, const Z);

static void *f, *g, *h;

static void done(void)
{
  cudaFree(f);
  cudaFree(g);
  free(h);
}

R *initialize_modules(const Z nx, const Z ny, const Z nz)
{
  cudaError_t err;

  const Z n = nx * ny * nz;
  const Z m = n + (nx * ny + ny * nz + nz * nx) * (2 * RADIUS);

  err = cudaMalloc(&f, sizeof(R) * m * N_VAR);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMalloc(&g, sizeof(R) * n * N_VAR);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  h = malloc(sizeof(R) * n * N_VAR);
  if(!h) error("fail to allocate host memory");

  atexit(done);

  initialize_io               (h, nx, ny, nz);
  initialize_initial_condition(h, nx, ny, nz);
  initialize_rk_2n            (g, nx, ny, nz);
  initialize_pde              (   nx, ny, nz);

  return (R *)f;
}
