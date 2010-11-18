#include "gpen.h"

static R *host;
static Z  nx, ny, nz;

void initialize_initial_condition(void *h, const Z n, const Z m, const Z l)
{
  host = (R *)h;

  nx = n;
  ny = m;
  nz = l;
}

void initial_condition(R *f, Q (*f0)(R, R, R))
{
  cudaError_t err;

  const Z n = nx * ny * nz;

  R *lnrho = host + 0 * n;
  R *ux    = host + 1 * n;
  R *uy    = host + 2 * n;
  R *uz    = host + 3 * n;

  Z i, j, k;

  for(k = 0; k < nz; ++k) {
    const R z = (R)(k) / nz;
    for(j = 0; j < ny; ++j) {
      const R y = (R)(j) / ny;
      for(i = 0; i < nx; ++i) {
        const R x = (R)(i) / nx;
        const Z l = (k * ny + j) * nx + i;
        const Q f = f0(x, y, z);
        lnrho[l] = f.lnrho;
        ux   [l] = f.ux   ;
        uy   [l] = f.uy   ;
        uz   [l] = f.uz   ;
      }
    }
  }

  err = cudaMemcpy(f, host, sizeof(R) * n * N_VAR, cudaMemcpyHostToDevice);
  if(cudaSuccess != err) error(cudaGetErrorString(err));
}
