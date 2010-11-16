#include "gpen.h"

static void *f, *g, *h;

static void done(void)
{
  cudaFree(f);
  cudaFree(g);
  free(h);
}

void initialize_modules(const Z nx, const Z ny, const Z nz)
{
  cudaError_t err;

  const Z n = nx * ny * nz;

  err = cudaMalloc(&f, sizeof(Q) * n);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMalloc(&g, sizeof(Q) * n);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  h = malloc(sizeof(Q) * n);
  if(!h) error("fail to allocate host memory");

  atexit(done);

  /* TODO: set modules */
}
