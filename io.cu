#include <stdio.h>
#include "gpen.h"

static Q *host;
static Z  nx, ny, nz;

void initialize_io(void *h, const Z n, const Z m, const Z l)
{
  host = (Q *)h;

  nx = n;
  ny = m;
  nz = l;
}

Z output(Z i, const Q *f)
{
  cudaError_t err;

  const Z n = nx * ny * nz;

  char  name[256];
  FILE *file;

  err = cudaMemcpy(host, f, sizeof(Q) * n, cudaMemcpyDeviceToHost);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  sprintf(name, "%04d.raw", i);
  file = fopen(name, "wb");
  fwrite(&nx,  sizeof(Z), 1, file);
  fwrite(&ny,  sizeof(Z), 1, file);
  fwrite(&nz,  sizeof(Z), 1, file);
  fwrite(host, sizeof(Q), n, file);
  fclose(file);

  return i;
}
