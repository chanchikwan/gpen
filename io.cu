#include <stdio.h>
#include "gpen.h"

static Q *host;
static Z  n;

void initialize_io(void *h, const Z nx, const Z ny, const Z nz)
{
  host = (Q *)h;
  n    = nx * ny * nz;
}

Z output(Z i, const Q *f)
{
  cudaError_t err;

  char  name[256];
  FILE *file;

  err = cudaMemcpy(host, f, sizeof(Q) * n, cudaMemcpyDeviceToHost);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  sprintf(name, "%04d.raw", i);
  file = fopen(name, "wb");
  fwrite(host, sizeof(Q), n, file);
  fclose(file);

  return i;
}
