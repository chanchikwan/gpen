#include "gpen.h"

static Q *host;
static Z  nx, ny, nz;

void initialize_initial_condition(void *h, const Z n, const Z m, const Z l)
{
  host = (Q *)h;

  nx = n;
  ny = m;
  nz = l;
}

void initial_condition(Q *f, Q (*f0)(R, R, R))
{
  cudaError_t err;

  const Z n = nx * ny * nz;

  Z i, j, k;

  for(k = 0; k < nz; ++k) {
    const R z = (R)(k) / nz;
    for(j = 0; j < ny; ++j) {
      const R y = (R)(j) / ny;
      for(i = 0; i < nx; ++i) {
        const R x = (R)(i) / nx;
        host[(k * ny + j) * nx + i] = f0(x, y, z);
      }
    }
  }

  err = cudaMemcpy(f, host, sizeof(Q) * n, cudaMemcpyHostToDevice);
  if(cudaSuccess != err) error(cudaGetErrorString(err));
}
