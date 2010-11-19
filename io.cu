#include <stdio.h>
#include "gpen.h"

static R *host;
static Z  nx, ny, nz;

void initialize_io(void *h, const Z n, const Z m, const Z l)
{
  host = (R *)h;

  nx = n;
  ny = m;
  nz = l;
}

Z output(Z i, const R *f)
{
  cudaError_t err;

  const Z n = nx * ny * nz;
  const Z m = n + (nx * ny + ny * nz + nz * nx) * (2 * RADIUS);

  char  name[256];
  FILE *file;

  Z h;

  for(h = 0; h < N_VAR; ++h) {
    err = cudaMemcpy(host + h * n,
                     f    + h * m + nx * ny * RADIUS,
                     sizeof(R) * n, cudaMemcpyDeviceToHost);
    if(cudaSuccess != err) error(cudaGetErrorString(err));
  }

  sprintf(name, "%04d.raw", i);
  file = fopen(name, "wb");
  fwrite(&nx,  sizeof(Z), 1, file);
  fwrite(&ny,  sizeof(Z), 1, file);
  fwrite(&nz,  sizeof(Z), 1, file);
  fwrite(host, sizeof(R), n * N_VAR, file);
  fclose(file);

  return i;
}
