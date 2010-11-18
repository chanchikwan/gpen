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

  char  name[256];
  FILE *file;

  err = cudaMemcpy(host, f, sizeof(R) * n * N_VAR, cudaMemcpyDeviceToHost);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  sprintf(name, "%04d.raw", i);
  file = fopen(name, "wb");
  fwrite(&nx,  sizeof(Z), 1, file);
  fwrite(&ny,  sizeof(Z), 1, file);
  fwrite(&nz,  sizeof(Z), 1, file);
  fwrite(host, sizeof(R), n * N_VAR, file);
  fclose(file);

  return i;
}
